
#'@description Descriptive statistics across cohorts 

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))
source(here::here("R","surv.R"))

# Load data ----
raids <- read_processed("raids")
tcga <- read_processed("tcga")
combined <- read_processed("combined")
combined_sets <- list(
  clin           = combined$clin, 
  dna            = combined$clin %>% drop_na(dna),
  rna            = combined$clin %>% drop_na(rna),
  clin_dna_rna   = combined$clin %>% drop_na(clin, dna, rna)
)

## Venn diagrams ----

venn_diagram_sequencing <- combined$clin %>%
  group_split_custom(cohort2)%>%
  map(~plot_venn_diagram_samples(.x, type = "sequencing"))

venn_diagram_sequencing %>%
  iwalk(~save_plot(.x, here::here("outputs","figures"),
                   paste0(.y, "_venn_diagram_sequencing"), 7, 5))

cowplot::plot_grid(venn_diagram_sequencing$raids,
                   venn_diagram_sequencing$tcga,
                   labels = c("Bio-RAIDs", "TCGA-CESC"),
                   ncol = 2, align = "hv", hjust = -3, vjust = 1)%>%
  save_plot(here::here("outputs","figures"),
            "compare_venn_diagram_sequencing", 13, 5)

venn_diagram_omics <-  combined$clin %>%
  group_split_custom(cohort2)%>%
  map(~plot_venn_diagram_samples(.x, type = "omics"))

venn_diagram_omics %>%
  iwalk(~save_plot(.x, here::here("outputs","figures"),
                   paste0(.y, "_venn_diagram_omics"), 6, 5))

cowplot::plot_grid(venn_diagram_omics$raids,
                   venn_diagram_omics$tcga,
                   labels = c("Bio-RAIDs", "TCGA-CESC"),
                   ncol = 2, align = "hv", hjust = -2, vjust = 1)%>%
  save_plot(here::here("outputs","figures"),
            "compare_venn_diagram_omics", 10, 5)

## Patients characteristics ----

tbl_compare_pat_char <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, age_c_f, figo_c_f, histo_f, hpv_negative_f)%>%
  tbl_summary(by = cohort, digits = list(everything()~0), 
              type = list(c(hpv_negative_f) ~ "categorical"))%>%
  add_p()%>%
  bold_p()%>%
  modify_header(label="")%>%
  bold_labels()%>%
  modify_footnote(everything()~NA)

tbl_compare_pat_char_omics_sets <- list(
  combined$clin,
  combined$clin %>% drop_na(dna),
  combined$clin %>% drop_na(rna), 
  combined$clin %>% drop_na(clin, dna, rna)
)%>%
  map(~.x %>% 
        dplyr::select(cohort, age_c_f, figo_c_f, histo_f, hpv_negative_f)%>%
        tbl_summary(by = cohort, digits = list(everything()~0), 
                    type = list(c(hpv_negative_f) ~ "categorical"))%>%
        add_p()%>%
        bold_p())%>%
  tbl_merge(tab_spanner = c("**Whole cohort**",
                            "**Genomic data set**",
                            "**Transcriptomic data set**",
                            "**Analysis set**"))%>%
  modify_header(label="")%>%
  bold_labels()%>%
  modify_footnote(everything()~NA)

tbl_combined_pat_char <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, age_c_f, histo_f, hpv_negative_f, figo_c_f) %>%
  tbl_strata(
    strata = figo_c_f,
    .header = "**FIGO stage {strata}**<br>N = {n}", 
    .tbl_fun = ~ .x %>%
      tbl_summary(
        by = cohort,
        type = list(hpv_negative_f ~ "categorical")
      ) %>%
      add_overall() %>%
      add_p() %>%
      bold_p() %>%
      modify_header(label = "") %>%
      bold_labels() %>%
      modify_footnote(everything() ~ NA)
  )

## Pathway-level frequencies of DNA alterations ----

plot_dichotomous(combined$pathway$clin_dna_rna,
                 var_prefix = "genomic_pathway_",
                 with_group = TRUE,
                 group_var = "cohort",
                 xlab = "",
                 ylab = "Pathway alteration frequency",
                 fill_values = c("Bio-RAIDs" = "#1b9e77",
                                 "TCGA-CESC" = "#377eb8"),
                 legend.position = "top",
                 legend.direction = "horizontal",
                 process_panel = F)%>%
  save_plot(here::here("outputs","figures"), 
            "compare_barplot_DNA_pathways", 12, 12)

tbl_compare_DNA_pathways <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, starts_with("genomic_pathway"))%>%
  tbl_summary(by = cohort, digits = everything()~c(0,0))%>%
  add_p(test = everything()~"fisher.test")%>%
  modify_header(label = "")%>%
  bold_labels()

tbl_combined_DNA_pathways <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, figo_c_f, starts_with("genomic_pathway"))%>%
  tbl_strata(
    strata = figo_c_f,
    .header = c("**FIGO stage {strata}**<br>N = {n}"), 
    .tbl_fun = ~.x %>%
      tbl_summary(by = cohort, digits = everything()~c(0,0))%>%
      add_p(test = everything()~"fisher.test")%>%
      modify_header(label = "")%>%
      bold_labels()%>%
      add_overall()
    )

## RNA-based pathway activity scores ----

plot_continuous(combined$pathway$clin_dna_rna,
                var_prefix = "hallmark_",
                with_group = TRUE,
                group_var = "cohort",
                xlab = "",
                ylab = "Pathway-activity score",
                fill_values = c("Bio-RAIDs" = "#1b9e77",
                                "TCGA-CESC" = "#377eb8"),
                legend.position = "top",
                legend.direction = "horizontal",
                process_panel = T)%>%
  save_plot(here::here("outputs","figures"), 
            "compare_boxplot_RNA_pathways", 24, 20)

tbl_compare_RNA_pathways <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, starts_with("hallmark_"))%>%
  tbl_summary(by = cohort, digits = everything()~c(0,0))%>%
  add_p(test = everything()~"wilcox.test")%>%
  modify_header(label = "")%>%
  bold_labels()

tbl_combined_RNA_pathways <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort, figo_c_f, starts_with("hallmark_"))%>%
  tbl_strata(
    strata = figo_c_f,
    .header = c("**FIGO stage {strata}**<br>N = {n}"), 
    .tbl_fun = ~.x %>%
      tbl_summary(by = cohort, digits = everything()~c(0,0))%>%
      add_p(test = everything()~"wilcox.test")%>%
      modify_header(label = "")%>%
      bold_labels()%>%
      add_overall()
    )

## Association between clinical variables and DNA-based pathways ----

clin_dna_fisher_heatmap <- combined$pathway$clin_dna_rna %>%
  dplyr::select(cohort2, age, figo, hpv_negative, 
                starts_with("genomic_pathway_"))%>%
  group_split_custom(cohort2)%>%
  purrr::map(~plot_cor_dichotomous(
    .x, 
    var_prefix = "age|figo|hpv_negative|genomic_pathway_"
    ))

purrr::iwalk(
  list(
    pval = "heatmap_pval",
    or   = "heatmap_or"
  ),
  function(slot, type) {
    
    list(
      raids = clin_dna_fisher_heatmap$raids[[slot]],
      tcga  = clin_dna_fisher_heatmap$tcga[[slot]]
    ) %>%
      purrr::iwalk(~save_plot(
        .x,
        here::here("outputs", "figures"),
        paste0(.y, "_plot_cor_", type, "_clin_DNA_pathways"),
        10, 8
      ))
  }
)

## Association between clinical variables RNA-based pathways ----

params <- tibble::tibble(
  var = c("age_c_f", "figo_c_f", "hpv_negative_f"),
  var_label = c("Age","FIGO stage","HPV status"),
  fill = list(
    c("\u226450" = "#377EB8", ">50" = "#E41A1C"),
    c("I/II" = "#377EB8", "III/IV" = "#E41A1C"),
    c("Negative" = "#377EB8", "Positive" = "#E41A1C")
  )
)

purrr::iwalk(
  list(raids = raids$pathway$clin_dna_rna,
       tcga  = tcga$pathway$clin_dna_rna),
  function(dat, cohort_name) {
    purrr::pwalk(
      params,
      function(var, var_label, fill) {
        plot_continuous(
          dat,
          var_prefix = "hallmark_",
          with_group = TRUE,
          group_var = var,
          group_label = var_label,
          xlab = "",
          ylab = "Pathway-activity score",
          fill_values = fill,
          legend.position = "top",
          legend.direction = "horizontal",
          process_panel = TRUE
        ) %>%
          save_plot(
            here::here("outputs", "figures"),
            paste0(cohort_name, "_cor_boxplot_RNA_pathways_", var),
            24, 20
          )
      }
    )
  }
)

## Association between RNA-based pathways ----

combined$pathway$clin_dna_rna %>%
  group_split_custom(cohort2)%>%
  purrr::map(~.x %>%
               dplyr::select(starts_with("hallmark")) %>%
               plot_cor_continuous())%>%
  iwalk(~save_plot(.x, here::here("outputs","figures"),
                   paste0(.y, "_plot_cor_RNA_pathways"), 10, 10))

# Survival curves ---------------------------------------------------------

## Progression-free survival in the different subpopulations 
combined_sets %>%
  plot_surv(formula      = Surv(time, event) ~ cohort,
            data         = ., 
            pval         = T, 
            # legend.title = rep("Cohort", 4), 
            legend.title = rep("", 4), 
            palette      = c("#1b9e77","#377eb8"))%>%
  iwalk(~save_plot(., here::here("outputs","figures"),
                   paste0("compare_surv_plot_", .y),
                   7.3, 5.8, newpage = F))

## Progression-free survival on subsets defined by variables 
combined$pathway$clin_dna_rna %>%
  group_split_custom(cohort2)%>%
  map(~.x %>% 
        plot_surv(formula      = list(Surv(time, event) ~ figo_c_f,
                                      Surv(time, event) ~ age_c_f,
                                      Surv(time, event) ~ hpv_negative_f),
                  data         = ., 
                  pval         = T, 
                  legend.title = c("FIGO stage","Age in years","HPV status"),
                  palette      = c("#377EB8","#E41A1C"))%>%
        set_names(gsub("data","surv_plot_clin_dna_rna", names(.))))%>%
  purrr::list_flatten()%>%
  iwalk(~save_plot(., here::here("outputs","figures"), .y,
                   6.3, 6, newpage = F))

# Survival tables ---------------------------------------------------------

## Progression-free survival in the different subpopulations 
tbl_compare_surv <- combined_sets %>%
  map(~survfit(Surv(time, event) ~ cohort, data = .x)%>%
        tbl_survfit(times        = c(12,24), 
                    label_header = "**{time} months-PFS rate (95% CI)**", 
                    statistic    = "{estimate} [{conf.low}-{conf.high}]"))%>%
  tbl_stack(group_header = c("Clinical data sets",
                             "Genomic data sets",
                             "Transcriptomic data sets",
                             "Analysis sets"))%>%
  modify_table_body(~.x %>% dplyr::filter(label !="cohort"))%>%
  modify_header(label = "")

# Cox proportional hazards models -----------------------------------------

tbl_compare_cox_univ <- tbl_cox_strata(
  data   = combined$pathway$clin_dna_rna,
  time   = "time",
  event  = "event",
  strata = "cohort",
  covars = c("age_c_f", "figo_c_f", "hpv_negative_f"),
  model  = "univ"
)

tbl_compare_cox_multi <- tbl_cox_strata(
  data   = combined$pathway$clin_dna_rna,
  time   = "time",
  event  = "event",
  strata = "cohort",
  covars = c("age_c_f", "figo_c_f", "hpv_negative_f"),
  model  = "multi"
)

# Global test -------------------------------------------------------------

tbl_compare_global_test <- 
  list(
    raids = list(gene    = raids$gene$clin_dna_rna,
                 pathway = raids$pathway$clin_dna_rna),
    tcga = list(gene     = tcga$gene$clin_dna_rna,
                pathway  = tcga$pathway$clin_dna_rna)) %>%
  purrr::imap(~compute_global_test(
    data_gene            = .x$gene, 
    data_pathway         = .x$pathway,
    time                 = "time",
    event                = "event",
    clin                 = c("age","figo","hpv_negative"), 
    dna_pattern          = "altered_",
    dna_pathway_pattern  = "genomic_pathway", 
    rna_pattern          = "rna_seq_", 
    rna_pathway_pattern  = "hallmark_"
  )
  )%>%
  imap(~{.x %>% rename_at(vars(-data_type), 
                          function(x) paste(x, .y, sep = "_"))})%>%
  purrr::reduce(inner_join, by = "data_type")%>%
  gt::gt()%>%
  gt::tab_spanner(label = gt::md("**Bio-RAIDs**"), columns = matches("Bio-RAIDs"))%>%
  gt::tab_spanner(label = gt::md("**TCGA-CESC**"), columns = matches("TCGA"))%>%
  gt::cols_label(data_type = gt::md("**Data type**"), 
                 starts_with("x_cov")~ gt::md("**N.<br>features**"),
                 starts_with("statistic")~ gt::md("**Statistic**"),
                 starts_with("expected")~ gt::md("**Expected**"),
                 starts_with("std_dev")~ gt::md("**Standard<br>deviation**"),
                 starts_with("p_value")~ gt::md("**p-value**"))%>%
  gt::tab_style(style = list(gt::cell_text(weight = "bold")),
                locations = gt::cells_column_labels())%>%
  gt::tab_style_body(style = gt::cell_text(weight = "bold"),
                     columns =  starts_with("p"),
                     fn = function(x) x < 0.05)%>%
  gt::fmt_number(decimals = 2, columns = -c(matches(c("p_val","x_cov"))))%>%
  gt::fmt(columns = starts_with("p_val"), 
          fns = function(x) ifelse(
            x < 0.05, formatC(x, format = "f", digits = 3),
            ifelse(x < 0.001, "<0.001",
                   formatC(x, format = "f", digits = 2))))%>%
  gt::tab_source_note(
    source_note = html(paste(
      "Clinical: clinical features",
      "DNA: gene-level DNA alterations features",
      "RNA: gene-level expression features",
      "DNA_pathways: DNA-based pathways features",
      "RNA_pathways: RNA-based pathways features",
      sep = "<br>"))
    )

# ---------------------------- Save all tables --------------------------#

for (i in setdiff(ls(pattern = "tbl_"), lsf.str())){
  saveRDS(get(i), here::here("outputs","tables", paste0(i,".rds")))
}
