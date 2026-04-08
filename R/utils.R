
# ----------------------------- Documentation ------------------------------#

#' This file centralizes utility functions that:
#' - Provide helpers for repetitive tasks related to environment setup or 
#'   configuration
#' - Standardize data reading and saving 
#' - Perform common data manipulations
#' - Standardize plotting or summary behaviours 

# ----------------------------- Dependencies -------------------------------#

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(TCGAbiolinks)
library(ggVennDiagram)

# Data import and export helpers ------------------------------------------

#' Import TCGA-CESC data (clinical, DNA, RNA, and external annotations)
#'
#' This function imports TCGA cervical cancer (TCGA-CESC) data, including
#' clinical, DNA, and RNA data retrieved from the GDC portal, and augments
#' them with external annotations such as HPV status and clinical endpoints.
#'
#' Data are loaded from a specified directory, with an option to trigger
#' download if files are not already available locally.
#'
#' @param raw_dir A character string specifying the directory containing
#'   raw TCGA data files. Default is `"data/raw/tcga"`.
#'
#' @param download Logical; if `TRUE`, data are downloaded from the GDC portal
#'   using `fetch_tcga_raw()`. If `FALSE`, existing local files are used.
#'   Default is `FALSE`.
#'
#' @return A named list containing:
#' \itemize{
#'   \item TCGA data imported via `fetch_tcga_raw()` (typically including
#'   clinical, DNA, and RNA datasets),
#'   \item `hpv`: a data.frame containing HPV status annotations,
#'   \item `endpoints`: a data.frame containing clinical endpoint data
#'   (e.g., progression-free survival).
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Calls `fetch_tcga_raw()` to retrieve TCGA-CESC data from the GDC
#'   portal (or load them locally).
#'   \item Imports HPV status from an external study
#'   (Panayiotou et al., 2020).
#'   \item Imports clinical endpoint data (including PFS) from another study
#'   (Liu et al., 2019), and filters for TCGA-CESC samples.
#'   \item Combines all datasets into a single named list.
#' }
#'
#' External files expected in `raw_dir`:
#' \itemize{
#'   \item `"ppat.1008468.s009.xlsx"`: HPV status annotations,
#'   \item `"NIHMS978596-supplement-1.xlsx"`: clinical endpoints.
#' }
#'
#' @export
import_tcga <- function(raw_dir = "data/raw/tcga",
                        download = FALSE){
  
  # Download and import TCGA-CESC data
  #'We download and import clinical, DNA and RNA data from the GDC portal. 
  tcga_gdc <- fetch_tcga_raw(
    directory = here::here(raw_dir), 
    download = download
  )
  
  # Import HPV status
  #'We import HPV status from the `@Panayiotou` et. al. study (2020)
  tcga_hpv <- readxl::read_xlsx(here::here(raw_dir, "ppat.1008468.s009.xlsx"))
  
  # Import clinical endpoints 
  #'We import PFS data from the `@Liu` et. al. study (2019)
  tcga_endpoints <- readxl::read_xlsx(
    here::here(raw_dir, "NIHMS978596-supplement-1.xlsx"),
    sheet = "ExtraEndpoints"
  ) %>%
    dplyr::filter(type == "CESC")%>%
    dplyr::select(-c(1,3))
  
  append(tcga_gdc, 
         list(hpv = tcga_hpv,
              endpoints = tcga_endpoints))
}

# Import Bio-RAIDs data

#' This function imports raw Bio-RAIDs data from multiple sources, including
#' clinical, genomic, transcriptomic, and proteomic datasets, and returns them
#' as a named list.
#'
#' The imported objects include raw tables related to patient/sample
#' correspondence, clinical annotations, HPV status, necrosis evaluation,
#' genomic profiling (WES, shallow WGS, mutational signatures, HRD, TMB/MSI),
#' RNA-seq, and RPPA proteomics.
#'
#' @param raw_dir A character string specifying the root directory containing
#'   the Bio-RAIDs raw data. Default is `"data/raw/raids"`.
#'
#' @return A named list of raw imported objects. Object names correspond to all
#'   variables ending with `"_raw"` created inside the function, including:
#' \itemize{
#'   \item clinical data tables,
#'   \item patient/sample correspondence tables,
#'   \item HPV and necrosis data,
#'   \item genomic data tables (HRD, TMB/MSI, WES mutations, mutational
#'   signatures),
#'   \item transcriptomic RNA-seq tables,
#'   \item proteomic RPPA tables.
#' }
#'
#' @details
#' The function imports data from several subdirectories within `raw_dir`:
#'
#' \itemize{
#'   \item `clinical_data`: patient/sample correspondence, clinical
#'   annotations, HPV status, necrosis evaluation, and additional clinical
#'   objects loaded from an `.Rdata` file,
#'   \item `genomic_data`: sequencing sample identifiers, HRD results, TMB/MSI
#'   results, WES mutation data, and mutational signature data,
#'   \item `transcriptomic_data`: RNA-seq expression tables,
#'   \item `proteomic_data`: RPPA protein expression data and RPPA-related
#'   patient annotations.
#' }
#'
#' At the end of the function, all objects whose names end with `"_raw"` are
#' collected and returned using `mget()`.
#'
#' @export
import_raids <- function(raw_dir = "data/raw/raids"){
  
  # Table of correspondance between patient ID and sample ID 
  clin_dir <- file.path(raw_dir, "clinical_data")
  patient_sample_raw <- readxl::read_xlsx(
    here::here(clin_dir, "RAIDs all patients ID analyzed.xlsx"),
    na = "NA"
    )
  
  #  Clinical data 
  clin_raw <- readxl::read_xlsx(
    here::here(clin_dir, "Listing_clinique_complet_RAIDS_23042024.xlsx"), 
    na = "NA"
    )

  hpv_raw <- read.csv(here::here(clin_dir, "RAIDs_HPV_377_patients.csv"))
  necrosis_raw <- readxl::read_xlsx(
    here::here(clin_dir, 
              "Evaluation RAIDS project LNS_01042021_fichier propre_01042021_LLC.xlsx")
    )
  load(here::here(clin_dir, "BIORAID_20210108.Rdata"))

  
  # Genomic data (WES, sWGS) 
  dna_dir <- here::here(raw_dir, "genomic_data")
  patient_sample_seq_raw <- read.csv(here::here(dna_dir, "patientID.csv"))
  hrd_raw <- read.csv(
    here::here(dna_dir, "RAIDS_shallowHRD_results_on_302_patients.csv")
    )
  tmb_msi_raw <- read.csv(
    here::here(dna_dir, "RAIDs_WES_majorSBS_TMB_MSI.csv")
    )
  mut_wes_raw <- read.csv(
    here::here(dna_dir, "RAIDs_WES_DecisionAlgo_OncoKBgenes_Pathways.csv")
    )
  signatures_subgroups_raw <- read.csv(
    here::here(
      dna_dir,
      "mutational_signatures_SBS_by_subgroups_without_artifact_and_unknown_signatures.csv")
    )
  signatures_raw <- read.csv(
    here::here(
      dna_dir, 
      "mutational_signatures_SBS_without_artifact_and_unknown_signatures.csv")
    )
  
  # Transcriptomic data (RNA-Seq) 
  rna_dir <- here::here(raw_dir, "transcriptomic_data")
  rna_seq_coding_raw <- read.csv(
    here::here(rna_dir, "RAIDs_RNAseq_tpm_270_patients_19594_coding_genes.csv")
    )
  rna_seq_gene_raw <- read.csv(
    here::here(rna_dir, "tableannot.csv")
    )
  
  # Proteomic data (RPPA) 
  prot_dir <- here::here(raw_dir, "proteomic_data")
  rppa_raw <- 
    readxl::read_xls(
      here::here(prot_dir, 
                "SuperCurve_Corrected_Normalized_data_RAIDs-col_utérus.xls")
    )
  rppa_patient_raw <- 
    readxl::read_xlsx(
      here::here(prot_dir, "clusterParys_vs_clusterIsabel.xlsx")
    )
  
  mget(ls(pattern = "_raw"))
}

#' Read processed cohort data from RDS file
#'
#' This function loads processed data for a given cohort from an `.RDS` file
#' stored in a standardized directory structure.
#'
#' @param cohort A character string specifying the cohort name. This should
#'   match the subdirectory name and file prefix (e.g., `"tcga"`, `"raids"`).
#'
#' @param processed_dir A character string specifying the root directory where
#'   processed data are stored. Default is `"data/processed"`.
#'
#' @return An R object read from the corresponding `.RDS` file, typically a
#'   list containing processed datasets ready for analysis.
#'
#' @export
read_processed <- function(cohort, processed_dir = "data/processed"){
  readRDS(here::here(processed_dir, cohort, paste0(cohort, "_processed_data.RDS")))
}

#' Integrate multi-omics data into gene-level and pathway-level datasets
#'
#' This function concatenates clinical, DNA, and RNA data into integrated
#' patient-level datasets for downstream analyses. It builds gene-level and
#' pathway-level datasets by joining available omics layers in different
#' combinations: clinical + DNA, clinical + RNA, and clinical + DNA + RNA.
#'
#' Optional preprocessing can be applied to remove non-informative DNA and RNA
#' features before returning the integrated datasets.
#'
#' @param id_col A character vector specifying the identifier columns used to
#'   join datasets. Default is `c("patient_id", "cohort")`.
#'
#' @param clin A data.frame or tibble containing clinical variables and the
#'   identifier columns.
#'
#' @param dna An optional data.frame or tibble containing DNA-derived variables
#'   and the identifier columns. Default is `NULL`.
#'
#' @param rna_gene An optional data.frame or tibble containing gene-level RNA
#'   variables and the identifier columns. Default is `NULL`.
#'
#' @param gsva_clin_rna An optional data.frame or tibble containing pathway-level
#'   RNA variables (e.g., GSVA scores) and the identifier columns. Default is
#'   `NULL`.
#'
#' @param gsva_clin_dna_rna An optional data.frame or tibble containing
#'   pathway-level RNA variables to be combined with clinical and DNA pathway
#'   data. Default is `NULL`.
#'
#' @param keep_only_complete Logical; if `TRUE`, datasets are joined using
#'   `dplyr::inner_join()` so that only patients present in all joined tables
#'   are retained. If `FALSE`, `dplyr::full_join()` is used instead. Default is
#'   `TRUE`.
#'
#' @param pre_filter Logical; if `TRUE`, applies feature filtering to the
#'   integrated gene-level datasets:
#'   \itemize{
#'     \item removes non-altered DNA gene features,
#'     \item keeps only the main mutational signature groups,
#'     \item removes lowly expressed RNA genes.
#'   }
#'   Default is `TRUE`.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item `clin`: the input clinical dataset,
#'   \item `gene`: a list of gene-level integrated datasets:
#'   \itemize{
#'     \item `clin_dna`,
#'     \item `clin_rna`,
#'     \item `clin_dna_rna`
#'   }
#'   \item `pathway`: a list of pathway-level integrated datasets:
#'   \itemize{
#'     \item `clin_dna`,
#'     \item `clin_rna`,
#'     \item `clin_dna_rna`
#'   }
#' }
#' 
#' @export
integrate_data <- function(
    id_col = c("patient_id","cohort"),
    clin,
    dna = NULL,
    rna_gene = NULL,
    gsva_clin_rna = NULL,
    gsva_clin_dna_rna = NULL, 
    keep_only_complete = TRUE,
    pre_filter = TRUE
){
  
  .join <- function(x, y){
    if (is.null(y)) return(x)
    if (keep_only_complete) {
      dplyr::inner_join(x, y, by = id_col)
    } else {
      dplyr::full_join(x, y, by = id_col)
    }
  }
  
  .integrate <- function(clin_df, dna_df, rna_df){
    out <- clin_df %>% 
      # Keep only complete cases (necrosis for Bio-RAIDs only)
      drop_na(any_of(c("age","figo","hpv_negative","necrosis"))) 
    out <- .join(out, dna_df)
    out <- .join(out, rna_df)
    out
  }
  
  .select_prefix <- function(df, prefix){
    df %>%
      dplyr::select(all_of(id_col), 
                    dplyr::matches(paste0("^(", prefix, ")")))
  }
  
  dna_gene <- .select_prefix(
    dna, 
    "altered_|cnv_|snv_|fusion_|hrd|tmb|msi|sig_group_"
    )
  rna_gene <- .select_prefix(
    rna_gene, 
    "rna_seq_")
  dna_pathway <- .select_prefix(
    dna, 
    "genomic_pathway_|hrd|tmb"
    )
  gsva_clin_rna <- .select_prefix(
    gsva_clin_rna, 
    "hallmark_"
    )
  gsva_clin_dna_rna <- .select_prefix(
    gsva_clin_dna_rna, 
    "hallmark_"
    )
    
  # Gene-level aggregation 
  gene_clin_dna     <- .integrate(clin, dna_gene, NULL)
  gene_clin_rna     <- .integrate(clin, NULL, rna_gene)
  gene_clin_dna_rna <- .integrate(clin, dna_gene, rna_gene)
  
  if (isTRUE(pre_filter)){
    gene_clin_dna <- gene_clin_dna %>%
      filter_out_non_altered_genes()%>%
      filter_main_sig_group()
    
    gene_clin_rna <- gene_clin_rna %>%
      filter_out_low_expressed_genes()
    
    gene_clin_dna_rna <- gene_clin_dna_rna %>%
      filter_out_non_altered_genes()%>%
      filter_main_sig_group()%>%
      filter_out_low_expressed_genes()
  }
  
  # Pathway-level aggregation
  pathway_clin_dna     <- .integrate(clin, dna_pathway, NULL)
  pathway_clin_rna     <- .integrate(clin, NULL, gsva_clin_rna)
  pathway_clin_dna_rna <- .integrate(clin, dna_pathway, gsva_clin_dna_rna)
  
  # Output 
  list(
    clin = clin, 
    gene = list(
      clin_dna     = gene_clin_dna,
      clin_rna     = gene_clin_rna,
      clin_dna_rna = gene_clin_dna_rna
      ),
    pathway = list(
      clin_dna     = pathway_clin_dna,
      clin_rna     = pathway_clin_rna,
      clin_dna_rna = pathway_clin_dna_rna
    )
  )
}

#' Combine integrated datasets from Bio-RAIDs and TCGA
#'
#' This function combines harmonized datasets from the Bio-RAIDs and TCGA
#' cohorts into a single multi-cohort object. Clinical, gene-level, and
#' pathway-level datasets are merged while preserving common variable labels.
#'
#' The function assumes that both inputs follow the same nested structure,
#' with components `clin`, `gene`, and `pathway`.
#'
#' @param raids A named list containing integrated Bio-RAIDs datasets. It is
#'   expected to include:
#'   \itemize{
#'     \item `clin`: clinical dataset,
#'     \item `gene`: a list of gene-level integrated datasets,
#'     \item `pathway`: a list of pathway-level integrated datasets.
#'   }
#'
#' @param tcga A named list containing integrated TCGA datasets with the same
#'   structure as `raids`.
#'
#' @return A named list with three elements:
#' \itemize{
#'   \item `clin`: the combined clinical dataset,
#'   \item `gene`: a list of combined gene-level datasets,
#'   \item `pathway`: a list of combined pathway-level datasets.
#' }
#'
#' Clinical datasets are combined directly. Gene-level and pathway-level
#' datasets are combined pairwise using `purrr::map2()`, assuming that the
#' corresponding elements in `raids$gene` and `tcga$gene` (and similarly for
#' `pathway`) are aligned in the same order and represent the same data type.
#' 
#' @export
combine_data <- function(raids, tcga){
  
  .combine <- function(x, y){
    common_cols <- intersect(names(x), names(y))
    labs <- var_label(x)[common_cols]
    out <- full_join(x, y, by = common_cols)
    var_label(out) <- labs 
    out
  }
  
  clin <- .combine(raids$clin, tcga$clin)
  
  gene <- purrr::map2(
    raids$gene, tcga$gene,
    ~ .combine(.x, .y)
    )
  
  pathway <- purrr::map2(
    raids$pathway, tcga$pathway,
    ~ .combine(.x, .y)
  )
  
  list(clin = clin, 
       gene = gene, 
       pathway = pathway)
}

# Data handling helpers ---------------------------------------------------

#' Split a data frame by grouping variables and assign names to groups
#'
#' This function splits a data.frame or tibble into a list of tibbles based on
#' one or more grouping variables, and assigns informative names to each list
#' element corresponding to the grouping levels.
#'
#' @param df A data.frame or tibble to split.
#'
#' @param ... One or more unquoted column names used as grouping variables.
#'   These variables define how the data will be split.
#'
#' @return A named list of tibbles, where each element corresponds to a group
#'   defined by the unique combinations of the grouping variables. List names
#'   are constructed by concatenating grouping values with `"_"`.
#'
#'@export
group_split_custom <- function(df, ...) {
  vars <- enquos(...)
  df %>%
    group_split(!!!vars, .keep = TRUE) %>%
    set_names(
      df %>%
        distinct(!!!vars) %>%
        mutate(name = apply(across(everything()), 1, paste, collapse = "_")) %>%
        pull(name)
    )
}

#' Retrieve variable labels from a dataset
#'
#' This function extracts variable labels from a data.frame or tibble using
#' the `labelled` package. If a variable does not have a label, its name is
#' used as a fallback.
#'
#' @param X A data.frame or tibble containing variables with optional labels.
#'
#' @return A tibble with two columns:
#' \itemize{
#'   \item `variable`: variable names,
#'   \item `label`: corresponding variable labels (or variable names if no
#'   label is defined).
#' }
.get_var_labels <- function(X) {
  require(labelled)
  labs <- vapply(names(X), function(v) {
    lab <- labelled::var_label(X[[v]])
    if (is.null(lab) || is.na(lab) || identical(lab, "")) v else as.character(lab)
  }, character(1))
  tibble::tibble(variable = names(X), label = unname(labs))
}

# Plotting helpers --------------------------------------------------------

#' Apply a consistent minimal theme for publication-ready plots
#'
#' This function standardizes the appearance of ggplot objects using a
#' clean minimal theme with large text sizes, improved facet display,
#' and customizable legend positioning.
#'
#' @param plot A ggplot object.
#' @param legend.position Character. Position of the legend
#'   (e.g., "right", "bottom", "none").
#' @param legend.direction Character. Layout of the legend
#'   (e.g., "vertical", "horizontal").
#' @param base_size Numeric. Base font size. Default is 20.
#'
#' @return A ggplot object with the applied theme.
#'
#' @export
base_theme <- function(plot,
                       legend.position = "right",
                       legend.direction = "vertical",
                       base_size = 20) {
  
  plot +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size),
      
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      
      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size
      ),
      strip.placement = "outside",
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        colour = NA
      ),
      
      plot.margin = ggplot2::margin(8, 8, 8, 8),
      panel.spacing = grid::unit(0.8, "lines"),
      
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      
      panel.border = ggplot2::element_rect(
        colour = "grey80",
        fill = NA,
        linewidth = 0.4
      )
    )
}

#' Save a plot in both PDF and PNG formats
#'
#' This function saves a plot object to both PDF and PNG files in a specified
#' directory. It ensures that the output directory exists and applies
#' consistent dimensions and resolution across formats.
#'
#' @param plot A plot object (typically a `ggplot2` object) to be saved.
#'
#' @param path A character string specifying the directory where the files
#'   will be saved. The directory is created if it does not exist.
#'
#' @param filename A character string specifying the base name of the output
#'   files (without extension).
#'
#' @param width Numeric value specifying the plot width.
#'
#' @param height Numeric value specifying the plot height.
#'
#' @param units A character string specifying the units for the PNG output.
#'   Default is `"px"`.
#'
#' @param res Numeric value specifying the resolution (in DPI) for the PNG
#'   output. Default is `300`.
#'
#' @param newpage Logical; passed to `print()` when saving the PDF. Controls
#'   whether a new page is started. Default is `TRUE`.
#'
#' @param ... Additional arguments passed to the graphics devices (`pdf()`
#'   and `png()`).
#'
#' @return Invisibly returns `NULL`. Files are written to disk.
#'
#' @export
save_plot <- function(plot, path, filename, width = 8, height = 6,
                      units = "px", res = 300, newpage = T, ...) {
  
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  # Define output file paths
  png_file <- file.path(path, paste0(filename, ".png"))
  pdf_file <- file.path(path, paste0(filename, ".pdf"))
  
  # Save PDF
  pdf(pdf_file, width = width, height = height, ...)
  print(plot, newpage = newpage)
  dev.off()
  
  # Save PNG
  png(png_file, width = width*res, height = height*res, unit = units, 
      res = res, ...)
  print(plot)
  dev.off()
  
  message("Saved as:\n- ", png_file, "\n- ", pdf_file)
}

#' Plot a Venn diagram of patient/sample availability
#'
#' This function creates a Venn diagram showing the overlap of patients across
#' different data modalities or sequencing data types. The sets are defined
#' from non-missing identifiers in selected columns of the input data frame.
#'
#' Two display modes are supported:
#' \itemize{
#'   \item `"sequencing"`: overlap between clinical, WES, WGS, and RNA-seq data,
#'   \item `"omics"`: overlap between clinical, genomic, and transcriptomic data.
#' }
#'
#' @param df A data.frame or tibble containing a `patient_id` column and the
#'   identifier columns required for the selected `type`.
#'
#' @param type A character string specifying the type of Venn diagram to draw.
#'   Supported values are:
#'   \itemize{
#'     \item `"sequencing"` for clinical, WES, WGS, and RNA-seq availability,
#'     \item `"omics"` for clinical, genomic, and transcriptomic availability.
#'   }
#'
#' @return A `ggplot2` object representing the Venn diagram.
#'
#' @details
#' The function first defines the relevant groups according to `type`:
#'
#' \itemize{
#'   \item for `"sequencing"`, the sets are based on the columns `clin`,
#'   `wes_id`, `wgs_id`, and `rna_id`,
#'   \item for `"omics"`, the sets are based on the columns `clin`, `dna`,
#'   and `rna`.
#' }
#'
#' For each group, patients with non-missing values in the corresponding column
#' are selected, and unique `patient_id` values are used as set elements.
#'
#' The Venn diagram is then generated with `ggVennDiagram::ggVennDiagram()`,
#' using type-specific set sizes and color scales.
#'
#' @export
plot_venn_diagram_samples <- function(df, type) {
  
  groups <- switch(
    type,
    sequencing = list(
      "Clinical\n data" = "clin",
      "WES data"        = "wes_id",
      "WGS data"        = "wgs_id",
      "RNA-Seq\n data"  = "rna_id"
    ),
    omics = list(
      "Clinical\n data"         = "clin",
      "Genomic\n data"          = "dna",
      "\n Transcriptomic data"  = "rna"
    )
  )
  
  sets <- purrr::map(
    groups,
    ~ df %>% tidyr::drop_na(dplyr::all_of(.x)) %>% 
      dplyr::pull(patient_id) %>% 
      unique()
  )
  
  p <- ggVennDiagram::ggVennDiagram(
    sets,
    label_alpha = 0,
    set_size    = if (type == "sequencing") 4.2 else 5.2,
    category.names = names(groups),
    set_color  = if (type == "omics") c("#0AB329",
                                        "#FF5733",
                                        "#0771A2") else "black"
  ) +
    (if (type == "omics")
      ggplot2::scale_fill_gradient(low = "#F4FAFE", 
                                   high = "#4981BF")
     else
       ggplot2::scale_fill_gradient(low = "grey90",  
                                    high = "blue")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_cartesian(clip = "off")
  
  p
}
