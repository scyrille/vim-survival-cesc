
# ----------------------------- Documentation ------------------------------#

#' This file contains functions to:
#' 1. Compute GSVA
#' 2. Get processes associated with all Hallmark gene sets
#' 3. Explore pathways distribution

# ----------------------------- Dependencies -------------------------------#

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(labelled)
library(janitor)
library(gtsummary)
library(msigdbr)   
library(DESeq2)   
library(GSVA)     

# ------------------------------ Functions ---------------------------------#

#' Compute GSVA pathway enrichment scores from gene-expression data
#'
#' This function computes sample-level pathway enrichment scores using the
#' Gene Set Variation Analysis (GSVA) framework and Hallmark gene sets from
#' the `msigdbr` package. It preprocesses gene-expression data, filters out
#' lowly expressed genes, applies variance-stabilizing transformation with
#' `DESeq2`, aggregates duplicated gene symbols, and returns pathway-level
#' GSVA scores in wide patient-level format.
#'
#' @param gene_expression A data.frame or tibble containing gene-expression
#'   data. It must include:
#'   \itemize{
#'     \item a gene identifier column named `gene_id`,
#'     \item a gene symbol column named `gene_name`,
#'     \item one column per patient/sample containing expression values.
#'   }
#'
#' @param cohort_name A character string indicating the cohort name. This value
#'   is added as a new `cohort` column in the output.
#'
#' @return A tibble in wide format with:
#' \itemize{
#'   \item one row per patient,
#'   \item one column per Hallmark pathway GSVA score,
#'   \item a `patient_id` column,
#'   \item a `cohort` column,
#'   \item cleaned column names using `janitor::clean_names()`,
#'   \item variable labels corresponding to pathway names.
#' }
#'
#' @export
compute_gsva <- function(gene_expression, 
                         cohort_name){
  
  #'To perform GSVA, we need:
  #'1. A normalized gene expression dataset, which can be provided in one of 
  #'   the following containers
  #'   - A `matrix` of expression values with genes corresponding to rows and
  #'     samples corresponding to columns
  #'2. A collection of gene sets, which can be provided in one of the 
  #'   following containers:
  #'   - A `list` object when each element names correspond to the names of 
  #'     the gene sets
  
  #'We can obtain the hallmark pathway information from the main function of 
  #'the `msigdbr` package [@Dolgalev 2020]
  hallmarks <- msigdbr::msigdbr(species = "Homo sapiens",
                                category = "H") # only hallmarks gene sets
  hallmarks_gene <- split(hallmarks$gene_symbol, 
                          hallmarks$gs_name)%>%
    lapply(function(x) unique(x))
  
  #' Metadata contains the list of patient ID.
  #' We ensure metadata and expression data are in the same order. 
  gene_expression <- gene_expression %>%
    dplyr::relocate(stringr::str_sort(names(.), numeric = TRUE))
  
  patient_id <- gene_expression %>%
    dplyr::select(-c(gene_name, gene_id))%>%
    colnames()
  metadata <- data.frame(patient_id = str_sort(patient_id, numeric = TRUE))
  
  #'We need to filter out the genes that have not been expressed or that 
  #'have low expression counts since we can not be as confident in those genes 
  #'being reliably measured. We are going to do some pre-filtering to keep only 
  #'genes with TPM>1 in >=20% of patients.
  
  gene_expression_filter <- gene_expression %>%
    # Only keep genes with TPM>1 in >=20% of patients
    dplyr::filter(rowSums(across(-c(gene_name, gene_id))>1) 
                  >= 0.2*ncol(gene_expression%>% 
                                dplyr::select(-c(gene_name, gene_id))))%>%
    # The next DESeq2 functions need the values to be converted to integers
    dplyr::mutate(across(-c(gene_name, gene_id), round))
  
  #'We will be using the `DESeq2` package for normalizing and transforming 
  #'our data, which requires us to format our data into a `DESeqDataSet` object.
  #'We turn the data frame (or matrix) into a [`DESeqDataSet` object] and 
  #'specify which variable labels our experimental groups using the 
  #'[`design` argument][@Love2014].
  #'In this chunk of code, we will not provide a specific model to the `design` 
  #'argument because we are not performing a differential expression analysis.
  
  gene_expression_matrix <- gene_expression_filter %>%
    # We keep only gene IDs 
    dplyr::select(-gene_name)%>%
    # We need to store our gene identifiers as row names
    tibble::column_to_rownames("gene_id") %>%
    # Now we can convert our object into a matrix
    as.matrix()
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = gene_expression_matrix, # Our prepped data frame with counts
    colData = metadata, # Data frame with annotation for our samples
    design = ~1 # Here we are not specifying a model
  )
  
  #'We are going to use the `vst()` function from the `DESeq2` package 
  #'to normalize and transform the data.
  dds_norm <- vst(dds)
  
  #'We need to extract the normalized counts to a matrix and make it 
  #'into a data frame so we can use with tidyverse functions later. 
  vst_df <- assay(dds_norm) %>%
    as.data.frame() %>% # Make into a data frame
    tibble::rownames_to_column("gene_id") # Make Gene IDs into their own column
  
  gene_expression_filter_norm <- 
    dplyr::inner_join(gene_expression %>% 
                        dplyr::select(gene_id, gene_name),
                      vst_df,
                      by = "gene_id")
  
  #'We will want to keep in mind that GSVA requires that data is in a 
  #'matrix with the gene identifiers as row names.
  #'In order to successfully turn our data frame into a matrix, we will 
  #'need to ensure that we do not have any duplicate gene identifiers.
  
  #'We need to calculate the gene mean for each tumor sample (or patient)
  gene_expression_filter_norm_nondup <- gene_expression_filter_norm %>%
    dplyr::group_by(gene_name)%>%
    dplyr::summarize_at(vars(-gene_id), mean)
  
  #'Now we should prep this data so GSVA can use it. 
  gene_expression_matrix <- gene_expression_filter_norm_nondup %>%
    tibble::column_to_rownames("gene_name")%>%
    # We transform data into a matrix
    as.matrix()
  
  #'GSVA fits a model and ranks genes based on their expression 
  #'level relative to the sample distribution [@Hanzelmann2013].
  #'The pathway-level score calculated is a way of asking how 
  #'genes _within_ a gene set vary as compared to genes that 
  #'are _outside_ of that gene set [@Malhotra2018].
  
  #'The idea here is that we will get pathway-level scores for each 
  #'sample that indicate if genes in a pathway vary concordantly in 
  #'one direction (over-expressed or under-expressed relative to 
  #'the overall population) [@Hanzelmann2013].
  #'This means that GSVA scores will depend on the samples included 
  #'in the dataset when you run GSVA; if you added more samples and 
  #'ran GSVA again, you would expect the scores to change [@Hanzelmann2013].
  
  #'The output is a gene set by sample matrix of GSVA scores.
  #'Let's perform GSVA using the `gsva()` function.
  
  gsvaParam <- GSVA::gsvaParam(exprData = gene_expression_matrix,
                               geneSets = hallmarks_gene,
                               kcdf = "Gaussian")
  gsva(gsvaParam)%>%
    as.data.frame()%>%
    tibble::rownames_to_column("hallmark_pathway")%>%
    tidyr::pivot_longer(-hallmark_pathway, names_to = "patient_id")%>%
    tidyr::pivot_wider(names_from = hallmark_pathway, 
                       values_from = value)%>%
    labelled::set_variable_labels(.labels = names(.))%>%
    janitor::clean_names()%>%
    dplyr::mutate(cohort = cohort_name, .before = patient_id)
}


#' Retrieve pathway-to-process annotations for DNA and RNA pathways
#'
#' This function returns a reference table linking pathway variables to
#' biological process categories for both DNA and RNA pathway features.
#'
#' The output combines:
#' \itemize{
#'   \item DNA pathway categories (manually defined),
#'   \item RNA pathway categories based on Hallmark gene sets.
#' }
#'
#' This table can be used for annotation, grouping, visualization, or
#' downstream enrichment and interpretation analyses.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item `num`: numeric identifier of the pathway,
#'   \item `variable`: pathway name,
#'   \item `process`: associated biological process category,
#'   \item `data_type`: type of pathway (`"DNA_pathways"` or `"RNA_pathways"`).
#' }
#'
#' @export
get_pathways_process <- function(){
  pathways_df <- bind_rows(
    tibble(num = 1:26,
           variable =  c("Adhesion/Migration",
                         "Apoptosis",
                         "Cell cycle",
                         "Cellular Metabolism", 
                         "Chromatin organization",
                         "Development ==> Notch",
                         "DNA repair", 
                         "Epigenetics",
                         "Genome integrity ==> p53",
                         "Hedgehog signaling pathway",
                         "Hippo signaling pathway",
                         "Immunity",
                         "JAK-STAT signaling pathway",
                         "mRNA Processing",
                         "MAPK",
                         "Myc",
                         "Oxydative stress",
                         "PI3K",
                         "Proteins Processing",
                         "Regulation Of Gene Expression",
                         "RTK/RAS",
                         "Senescence",
                         "TGF-beta",
                         "Transcription factor-regulator",
                         "Wnt Signaling Pathway",
                         "Others"),
           process = "",
           data_type = "DNA_pathways"),
    tibble(num = 1:50,
           variable =  c("HALLMARK_APICAL_JUNCTION",
                         "HALLMARK_APICAL_SURFACE",
                         "HALLMARK_PEROXISOME",
                         "HALLMARK_ADIPOGENESIS",
                         "HALLMARK_ANGIOGENESIS",
                         "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                         "HALLMARK_MYOGENESIS",
                         "HALLMARK_SPERMATOGENESIS",
                         "HALLMARK_PANCREAS_BETA_CELLS",
                         "HALLMARK_DNA_REPAIR",
                         "HALLMARK_UV_RESPONSE_DN",
                         "HALLMARK_UV_RESPONSE_UP",
                         "HALLMARK_ALLOGRAFT_REJECTION",
                         "HALLMARK_COAGULATION",
                         "HALLMARK_COMPLEMENT",
                         "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                         "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                         "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                         "HALLMARK_INFLAMMATORY_RESPONSE",
                         "HALLMARK_BILE_ACID_METABOLISM",
                         "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
                         "HALLMARK_FATTY_ACID_METABOLISM",
                         "HALLMARK_GLYCOLYSIS",
                         "HALLMARK_HEME_METABOLISM",
                         "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                         "HALLMARK_XENOBIOTIC_METABOLISM",
                         "HALLMARK_APOPTOSIS",
                         "HALLMARK_HYPOXIA",
                         "HALLMARK_PROTEIN_SECRETION",
                         "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                         "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                         "HALLMARK_E2F_TARGETS",
                         "HALLMARK_G2M_CHECKPOINT",
                         "HALLMARK_MYC_TARGETS_V1",
                         "HALLMARK_MYC_TARGETS_V2",
                         "HALLMARK_P53_PATHWAY",
                         "HALLMARK_MITOTIC_SPINDLE",
                         "HALLMARK_ANDROGEN_RESPONSE",
                         "HALLMARK_ESTROGEN_RESPONSE_EARLY",
                         "HALLMARK_ESTROGEN_RESPONSE_LATE",
                         "HALLMARK_IL2_STAT5_SIGNALING",
                         "HALLMARK_KRAS_SIGNALING_UP",
                         "HALLMARK_KRAS_SIGNALING_DN",
                         "HALLMARK_MTORC1_SIGNALING",
                         "HALLMARK_NOTCH_SIGNALING",
                         "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                         "HALLMARK_HEDGEHOG_SIGNALING",
                         "HALLMARK_TGF_BETA_SIGNALING",
                         "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                         "HALLMARK_WNT_BETA_CATENIN_SIGNALING"),
           process = c(rep("Cellular component", 3),
                       rep("Development", 6),
                       rep("DNA damage", 3),
                       rep("Immune", 7),
                       rep("Metabolic", 7),
                       rep("Pathway", 5),
                       rep("Proliferation", 6),
                       rep("Signaling", 13)),
           data_type = "RNA_pathways")
  )
  return(pathways_df)
}

#' Split a plot into predefined biological process panels
#'
#' This function takes a long-format data frame containing pathway-level
#' variables and their associated biological processes, splits the data into
#' four predefined process groups, applies a plotting function to each subset,
#' and combines the resulting plots into a 2 x 2 multi-panel layout.
#'
#' The legend is kept only in the first panel and removed from the others.
#' If a data frame of statistical annotations is provided, it is split using
#' the same predefined process groups and passed to the plotting function.
#'
#' @param df_long A long-format data frame containing at least the columns
#'   `variable` and `process`. Each row corresponds to one observation for one
#'   variable.
#' @param fun A plotting function returning a ggplot object. This function must
#'   accept at least `df_long` as an argument, and optionally `stats_df` and
#'   graphical parameters passed through `...`.
#' @param stats_df Optional data frame containing statistical annotations
#'   (for example p-values or significance symbols). It must contain at least
#'   the columns `variable` and `process` if provided.
#' @param ... Additional arguments passed to `fun()`.
#'
#' @return A combined plot object produced by `egg::ggarrange()`.
#'
#' @details
#' The data are split into four panels using the following predefined process
#' groups:
#' \itemize{
#'   \item Panel 1: `"Cellular component"`, `"DNA damage"`, `"Pathway"`
#'   \item Panel 2: `"Proliferation"`, `"Development"`
#'   \item Panel 3: `"Metabolic"`, `"Immune"`
#'   \item Panel 4: `"Signaling"`
#' }
#'
#' Each panel is then faceted by `process` using `facet_wrap()` with a fixed
#' number of rows and columns:
#' \itemize{
#'   \item Panel 1: 3 rows, 1 column
#'   \item Panel 2: 2 rows, 1 column
#'   \item Panel 3: 2 rows, 1 column
#'   \item Panel 4: 1 row, 1 column
#' }
#'
#'
#' @export
plot_process_panel <- function(df_long, fun, stats_df = NULL,
                               legend.position, ...) {
  
  df_long <- df_long %>%
    dplyr::left_join(get_pathways_process(), by = "variable")
  
  if (!is.null(stats_df)) {
    stats_df <- stats_df %>%
      dplyr::left_join(get_pathways_process(), by = "variable")
  }
  
  df_p1 <- df_long %>%
    dplyr::filter(process %in% c("Cellular component", "DNA damage", "Pathway"))
  df_p2 <- df_long %>%
    dplyr::filter(process %in% c("Proliferation", "Development"))
  df_p3 <- df_long %>%
    dplyr::filter(process %in% c("Metabolic", "Immune"))
  df_p4 <- df_long %>%
    dplyr::filter(process %in% c("Signaling"))
  
  stats_p1 <- if (!is.null(stats_df)) 
    dplyr::filter(stats_df, process %in% c("Cellular component", 
                                           "DNA damage", "Pathway")) else NULL
  stats_p2 <- if (!is.null(stats_df)) 
    dplyr::filter(stats_df, process %in% c("Proliferation", "Development")) else NULL
  stats_p3 <- if (!is.null(stats_df)) 
    dplyr::filter(stats_df, process %in% c("Metabolic", "Immune")) else NULL
  stats_p4 <- if (!is.null(stats_df)) 
    dplyr::filter(stats_df, process %in% c("Signaling")) else NULL
  
  p1 <- fun(df_long = df_p1, stats_df = stats_p1,
            legend.position, ...) +
    ggplot2::facet_wrap(~ .data[["process"]], 
                        nrow = 3, ncol = 1, scales = "free_y")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  p2 <- fun(df_long = df_p2, stats_df = stats_p2, 
            legend.position = "none", ...) +
    ggplot2::facet_wrap(~ .data[["process"]], 
                        nrow = 2, ncol = 1, scales = "free_y")+
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank())
  p3 <- fun(df_long = df_p3, stats_df = stats_p3, 
            legend.position = "none", ...) +
    ggplot2::facet_wrap(~ .data[["process"]], 
                        nrow = 2, ncol = 1, scales = "free_y")
  p4 <- fun(df_long = df_p4, stats_df = stats_p4, 
            legend.position = "none", ...) +
    ggplot2::facet_wrap(~ .data[["process"]], 
                        nrow = 1, ncol = 1, scales = "free_y")
  
  egg::ggarrange(
    p1, p2, p3, p4,
    ncol = 2, nrow = 2,
    widths = c(1, 1), heights = c(1, 1)
  )
}

#' Plot continuous variables as boxplots
#'
#' This function reshapes continuous variables selected by a common prefix into
#' long format and displays them as horizontal boxplots, optionally stratified
#' by group. Variable names can be replaced by variable labels. A process-panel
#' layout can also be applied.
#'
#' @param df A data frame.
#' @param var_prefix Character. Prefix used to select continuous variables.
#' @param with_group Logical. Whether to stratify boxplots by group.
#' @param group_var Character. Name of the grouping variable.
#' @param group_label Character. Label of the grouping variable.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param fill_values Named character vector for group colors.
#' @param ylim Numeric vector of length 2 giving y-axis limits.
#' @param legend.position Character. Legend position.
#' @param legend.direction Character. Legend direction.
#' @param process_panel Logical. Whether to split the plot into process panels.
#'
#' @return A ggplot object or a combined plot.
#' @export
plot_continuous <- function(df,
                            var_prefix = "hallmark_",
                            with_group = TRUE,
                            group_var = "cohort",
                            group_label = "", 
                            xlab = NULL,
                            ylab = NULL,
                            fill_values = c("Bio-RAIDs" = "#1b9e77",
                                            "TCGA-CESC" = "#377eb8"),
                            ylim = c(-0.9, 0.9),
                            legend.position = "top",
                            legend.direction = "horizontal",
                            process_panel = TRUE) {
  
  vars <- names(dplyr::select(df, dplyr::starts_with(var_prefix)))
  
  df_long <- df %>%
    dplyr::select(
      dplyr::any_of(if (with_group) group_var),
      dplyr::all_of(vars)
    ) %>%
    {
      lab <- labelled::var_label(.[vars], unlist = TRUE)
      
      dplyr::rename_with(
        .,
        ~ dplyr::coalesce(lab[.x], .x),
        .cols = dplyr::all_of(vars)
      )
    } %>%
    tidyr::pivot_longer(
      cols = -dplyr::any_of(if (with_group) group_var),
      names_to = "variable",
      values_to = "value"
    )
  
  .boxplot <- function(df_long,
                       stats_df = NULL,
                       with_group,
                       group_var,
                       group_label, 
                       xlab,
                       ylab,
                       fill_values,
                       ylim,
                       legend.position,
                       legend.direction) {
    
    if (with_group) {
      plot <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
          x = .data[["variable"]],
          y = .data[["value"]],
          fill = .data[[group_var]]
        )
      ) +
        ggplot2::geom_boxplot() +
        ggpubr::stat_compare_means(
          ggplot2::aes(group = .data[[group_var]]),
          method = "wilcox.test",
          label = "p.signif",
          colour = "#F8766D"
        ) +
        ggplot2::scale_fill_manual(values = fill_values)
    } else {
      plot <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
          x = .data[["variable"]],
          y = .data[["value"]]
        )
      ) +
        ggplot2::geom_boxplot()
    }
    
    plot <- plot +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = xlab,
        y = ylab,
        fill = group_label
      )
    
    if (!is.null(ylim)) {
      plot <- plot + ggplot2::ylim(ylim)
    }
    
    base_theme(
      plot,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
  
  if (process_panel) {
    plot_process_panel(
      df_long = df_long,
      fun = .boxplot,
      with_group = with_group,
      group_var = group_var,
      group_label = group_label, 
      xlab = xlab,
      ylab = ylab,
      fill_values = fill_values,
      ylim = ylim,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  } else {
    .boxplot(
      df_long = df_long,
      with_group = with_group,
      group_var = group_var,
      group_label = group_label, 
      xlab = xlab,
      ylab = ylab,
      fill_values = fill_values,
      ylim = ylim,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
}

#' Plot proportions for dichotomous variables
#'
#' This function computes the proportion of positive values for dichotomous
#' variables selected by a common prefix, optionally stratified by group,
#' and displays them as horizontal bar plots. A Fisher or chi-squared test can
#' be added for group comparison. A process-panel layout can also be applied.
#'
#' @param df A data frame.
#' @param var_prefix Character. Prefix used to select dichotomous variables.
#' @param with_group Logical. Whether to stratify proportions by group.
#' @param group_var Character. Name of the grouping variable.
#' @param xlab Character. X-axis label.
#' @param ylab Character. Y-axis label.
#' @param fill_values Named character vector for group colors.
#' @param legend.position Character. Legend position.
#' @param legend.direction Character. Legend direction.
#' @param process_panel Logical. Whether to split the plot into process panels.
#'
#' @return A ggplot object or a combined plot.
#' @export
plot_dichotomous <- function(df,
                             var_prefix = "genomic_pathway_",
                             with_group = TRUE,
                             group_var = "cohort",
                             xlab = NULL,
                             ylab = NULL,
                             fill_values = c("Bio-RAIDs" = "#1b9e77",
                                             "TCGA-CESC" = "#377eb8"),                         
                             legend.position = "top",
                             legend.direction = "horizontal",
                             process_panel = FALSE) {
  
  vars <- names(dplyr::select(df, dplyr::starts_with(var_prefix)))
  
  df_raw <- df %>%
    dplyr::select(
      dplyr::any_of(if (with_group) group_var),
      dplyr::all_of(vars)
    ) %>%
    {
      lab <- labelled::var_label(.[vars], unlist = TRUE)
      
      dplyr::rename_with(
        .,
        ~ dplyr::coalesce(lab[.x], .x),
        .cols = dplyr::all_of(vars)
      )
    }
  
  value_vars <- setdiff(names(df_raw), if (with_group) group_var else character(0))
  
  df_long <- df_raw %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(if (with_group) group_var))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(value_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(value_vars),
      names_to = "variable",
      values_to = "value"
    )
  
  stats_df <- NULL
  
  if (with_group) {
    df_test <- df_raw %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(value_vars),
        names_to = "variable",
        values_to = "value"
      )
    
    stats_df <- df_test %>%
      dplyr::group_by(.data[["variable"]]) %>%
      dplyr::summarise(
        pval = stats::fisher.test(table(.data[[group_var]], 
                                        .data[["value"]]))$p.value,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        pval_stars = dplyr::case_when(
          .data[["pval"]] <= 0.001 ~ "***",
          .data[["pval"]] <= 0.01 ~ "**",
          .data[["pval"]] <= 0.05 ~ "*",
          .data[["pval"]] <= 0.1 ~ ".",
          TRUE ~ "ns"
        )
      )
  }
  
  .barplot <- function(df_long,
                       stats_df = NULL,
                       with_group,
                       group_var,
                       xlab,
                       ylab,
                       fill_values,
                       legend.position,
                       legend.direction) {
    
    if (with_group) {
      plot <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
          x = stats::reorder(.data[["variable"]], .data[["value"]]),
          y = .data[["value"]],
          fill = .data[[group_var]]
        )
      ) +
        ggplot2::geom_col(
          position = ggplot2::position_dodge(width = 0.8),
          colour = NA,
          linewidth = 0.2
        ) +
        ggplot2::scale_fill_manual(values = fill_values)
    } else {
      plot <- ggplot2::ggplot(
        df_long,
        ggplot2::aes(
          x = stats::reorder(.data[["variable"]], .data[["value"]]),
          y = .data[["value"]]
        )
      ) +
        ggplot2::geom_col(
          colour = NA,
          linewidth = 0.2
        )
    }
    
    plot <- plot +
      ggplot2::labs(x = xlab, y = ylab, fill = "") +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent)
    
    if (!is.null(stats_df)) {
      ymax <- max(df_long$value, na.rm = TRUE)
      
      plot <- plot +
        ggplot2::geom_text(
          data = stats_df,
          ggplot2::aes(
            x = .data[["variable"]],
            y = 1.05 * ymax,
            label = .data[["pval_stars"]]
          ),
          inherit.aes = FALSE,
          size = 6,
          colour = "#F8766D"
        ) +
        ggplot2::expand_limits(y = 0.6 * ymax)
    }
    
    base_theme(
      plot,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
  
  if (process_panel) {
    plot_process_panel(
      df_long = dplyr::left_join(df_long, get_pathways_process(), 
                                 by = "variable"),
      fun = .barplot,
      stats_df = if (!is.null(stats_df)) 
        dplyr::left_join(stats_df, get_pathways_process(), 
                         by = "variable") else NULL,
      with_group = with_group,
      group_var = group_var,
      xlab = xlab,
      ylab = ylab,
      fill_values = fill_values,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  } else {
    .barplot(
      df_long = df_long,
      stats_df = stats_df,
      with_group = with_group,
      group_var = group_var,
      xlab = xlab,
      ylab = ylab,
      fill_values = fill_values,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
}

#' Plot correlation heatmap between continuous pathway scores
#'
#' @param df A data.frame containing pathway variables.
#' @param var_prefix A character string (or multiple prefixes separated by "|")
#'   used to select pathway columns.
#'
#' @return A ggplot object containing a corrplot-based correlation heatmap.
#' 
#' @export
plot_cor_continuous <- function(df, var_prefix = "hallmark_") {

  data <- df %>%
    dplyr::select(dplyr::starts_with(unlist(strsplit(var_prefix, split = "[|]"))))

  colnames(data) <- var_label(data, unlist = TRUE)

  plot <- ggplotify::as.ggplot(~corrplot::corrplot(
    cor(data),
    method = "color",
    type = "lower",
    tl.cex = 0.5,
    tl.col = "grey30",
    diag = FALSE
  ))

  plot
}

#' Pairwise Fisher's exact tests with heatmap visualization
#'
#' Performs pairwise Fisher's exact tests between binary/categorical variables
#' selected from a data frame based on a prefix pattern. For each variable pair,
#' the function computes the raw p-value and odds ratio, and summarizes results
#' as symmetric matrices. It also generates heatmaps to visualize pairwise
#' associations.
#'
#' @param df A data.frame containing the variables of interest.
#' @param var_prefix A character string (or multiple prefixes separated by "|")
#'   used to select variables via `dplyr::starts_with()`.
#'
#' @return A list with:
#'   \item{mat_pval}{A symmetric matrix of raw p-values from Fisher's exact tests.}
#'   \item{heatmap_pval}{A ggplot2 heatmap of pairwise p-values.}
#'   \item{mat_or}{A symmetric matrix of odds ratios.}
#'   \item{heatmap_or}{A ggplot2 heatmap of odds ratios.}
#'
#' @export
plot_cor_dichotomous <- function(df, var_prefix = "genomic_pathway_") {
  
  data <- vars <- df %>% 
    dplyr::select(starts_with(unlist(strsplit(var_prefix, split = "[|]"))))
    
  colnames(data) <- var_label(data, unlist = T)
  vars <- colnames(data)

  combn_pairs <- combn(vars, 2, simplify = FALSE)
  
  results <- purrr::map_dfr(combn_pairs, function(pair) {
    
    x <- data[[pair[1]]]
    y <- data[[pair[2]]]
    
    tab <- table(x, y)
    
    test <- fisher.test(tab)
    
    tibble::tibble(
      var1 = pair[1],
      var2 = pair[2],
      p.value = test$p.value,
      odds.ratio = unname(test$estimate)
    )
  })
  
  # Create symmetric matrix of raw p-values
  mat_pval <- matrix(NA, length(vars), length(vars),
                     dimnames = list(vars, vars))
  
  for(i in seq_len(nrow(results))) {
    mat_pval[results$var1[i], results$var2[i]] <- results$p.value[i]
    mat_pval[results$var2[i], results$var1[i]] <- results$p.value[i]
  }
  diag(mat_pval) <- 0
  
  # Heatmap of raw p-values
  heatmap_pval <- mat_pval %>%
    as.table()%>%
    as.data.frame()%>%
    mutate(Freq = ifelse(Var1==Var2, NA, Freq))%>%
    {
      ggplot(., aes(Var1, Var2, fill = Freq)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C", direction = -1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = "Raw p-value", x ="", y = "")
    }
  
  # Create symmetric matrix of raw odds ratio
  mat_or <- matrix(NA, length(vars), length(vars),
                   dimnames = list(vars, vars))
  
  for(i in seq_len(nrow(results))) {
    mat_or[results$var1[i], results$var2[i]] <- results$odds.ratio[i]
    mat_or[results$var2[i], results$var1[i]] <- results$odds.ratio[i]
  }
  diag(mat_or) <- 0
  
  # Heatmap of raw p-values
  heatmap_or <- mat_or %>%
    as.table()%>%
    as.data.frame()%>%
    mutate(Freq = ifelse(Var1==Var2, NA, Freq))%>%
    {
      ggplot(., aes(Var1, Var2, fill = Freq)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C", direction = -1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(fill = "Raw odds ratios", x ="", y = "")
    }
  
  list(heatmap_pval = heatmap_pval, 
       heatmap_or = heatmap_or)
}
