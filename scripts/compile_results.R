
#'@description Compile results across cohorts 

library(here)

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))
source(here::here("R","cox_lasso.R"))
source(here::here("R","cox_boosting.R"))
source(here::here("R","surv.R"))
source(here::here("R","vimp.R"))

cohorts <-  c("raids","tcga")
cohorts_name <- c("Bio-RAIDs","TCGA-CESC")

pathways_process <- get_pathways_process()
labels <- setNames(as.list(pathways_process$label), 
                   pathways_process$variable_name)
dna_labels <- labels %>% subset(startsWith(names(.), "genomic_pathway"))
rna_labels <- labels %>% subset(startsWith(names(.), "hallmark"))
dna_process_labels <- pathways_process %>%
  dplyr::filter(data_type =='DNA_pathways')%>%
  dplyr::select(process)%>%
  distinct()%>%
  pull()
rna_process_labels <- pathways_process %>%
  dplyr::filter(data_type =='RNA_pathways')%>%
  dplyr::select(process)%>%
  distinct()%>%
  pull()

combined <- read_processed("combined")$pathway$clin_dna_rna %>%
  set_variable_labels(.labels = labels)
raids <- read_processed("raids")
tcga <- read_processed("tcga")
 
# Cox lasso ---------------------------------------------------------------

cox_lasso_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_lasso_fit.rds"))))%>%
  set_names(cohorts)#%>%
  # purrr::list_flatten()

## Plot selected coefficients 
cox_lasso_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lasso_coef"), 9, 6))

cox_lasso_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x, data_type = T))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lasso_coef_omics"), 9, 6))
               
## Table of selected coefficients
tbl_cox_lasso <- cox_lasso_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# Cox group lasso ---------------------------------------------------------

cox_group_lasso_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_group_lasso_fit.rds"))))%>%
  set_names(cohorts)#%>%
# purrr::list_flatten()

## Table of selected coefficients
tbl_cox_group_lasso <- cox_group_lasso_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# Cox priority lasso ------------------------------------------------------

cox_priority_lasso_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_priority_lasso_fit.rds"))))%>%
  set_names(cohorts)#%>%
  # purrr::list_flatten()

## Table of selected coefficients
tbl_cox_priority_lasso <- cox_priority_lasso_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# With unpenalized clinical features
cox_priority_lasso_unpen_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_priority_lasso_unpen_fit.rds"))))%>%
  set_names(cohorts)#%>%
# purrr::list_flatten()

## Table of selected coefficients
tbl_cox_priority_lasso_unpen <- cox_priority_lasso_unpen_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# Likelihood-based boosting -----------------------------------------------

## Gene-level ----

# Bio-RAIDs only (with necrosis, HRD, TMB, MSI, mutationals signatures)
cox_lboost_gene_fit <- readRDS(
  here::here("outputs/results/raids_cox_lboost_necrosis_fit.rds")
)

# Plot selected coefficients
cox_lboost_gene_fit %>%
  plot_coef()%>%
  save_plot(here::here("outputs","figures"), 
            paste0("raids_gene_clin_dna_rna_cox_lboost_coef"), 9, 6)

cox_lboost_gene_fit %>%
  plot_coef(data_type = T)%>%
  save_plot(here::here("outputs","figures"),
            paste0("raids_gene_clin_dna_rna_cox_lboost_coef_omics"), 9, 6)

# Plot coefficients paths
cox_lboost_gene_fit %>%
  plot_coef_paths()%>%
  save_plot(here::here("outputs","figures"), 
            "raids_gene_clin_dna_rna_cox_lboost_coef_paths", 6, 4)

# Plot coefficients updates at each boosting step until optimal step
cox_lboost_gene_fit %>%
  (\(fit_obj) 0:fit_obj$opt_step %>%
     purrr::walk(
       ~ plot_coef_updates(fit_obj, step_show = .x) %>%
         save_plot(
           here::here("outputs","figures"),
           paste0("raids_gene_clin_dna_rna_cox_lboost_coef_updates_", .x),
           9, 6)
     )
  )()

# Table with number selected variables
tbl_cox_lboost_gene_nselect <- cox_lboost_gene_fit %>%
  tbl_n_select_coef()

## Table of selected coefficients
tbl_cox_lboost_gene_coef <- cox_lboost_gene_fit %>%
  tbl_select_coef()

## Cox regression with selected variables
clin <- c("age_c_f", "figo_c_f", "hpv_negative_f","necrosis_f")
var_select <- cox_lboost_gene_fit$var_select %>%
  subset(!. %in% c("age","figo","hpv_negative","necrosis"))
raids_gene_cox_lboost <- raids$gene$clin_dna_rna %>%
  dplyr::select(time, event, all_of(clin),  all_of(var_select))%>%
  dplyr::mutate(across(starts_with("rna_seq"), 
                       ~as.numeric(scale(.x))))%>%
  labelled::set_variable_labels(
    .labels =  var_label(raids$gene$clin_dna_rna[,var_select])%>%
      purrr::map(~gsub("RNA-Seq - ", "", .x))
  )

tbl_cox_univ_multi_gene_cox_lboost_select <- list(
  tbl_cox(
    data   = raids_gene_cox_lboost, 
    time   = "time",
    event  = "event",
    covars = c(clin, var_select),
    model  = "univ"
  )%>%
    gtsummary::add_q(method = "BH") %>%
    gtsummary::bold_p(q = TRUE),
  tbl_cox(
    data   = raids_gene_cox_lboost, 
    time   = "time",
    event  = "event",
    covars = c(clin, var_select),
    model  = "multi"
  )%>%
    gtsummary::add_q(method = "BH") %>%
    gtsummary::bold_p(q = TRUE) %>%
    gtsummary::modify_column_hide(c(stat_n, stat_nevent))
)%>%
  gtsummary::tbl_merge(tab_spanner = c("**Univariate Cox regression**",
                                       "**Multivariate Cox regression**"))%>%
  gtsummary::modify_spanning_header(c(stat_n_1, stat_nevent_1)~NA)

## Pathway-level ---- 
cox_lboost_pathway_fit <- cohorts %>%
  purrr::map(~readRDS(here::here(
    "outputs","results",
    paste0(.x,"_cox_lboost_fit.rds"))))%>%
  set_names(cohorts)%>%
  purrr::list_flatten()

## Plot selected coefficients
height <- set_names(c(4,9,4,4,5,5), 
                    names(cox_lboost_pathway_fit))
cox_lboost_pathway_fit %>%
  purrr::map(~plot_coef(.x))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lboost_coef"), 9, height[.y]))

cox_lboost_pathway_fit %>%
  purrr::map(~plot_coef(.x, data_type = T))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lboost_coef_omics"), 9, height[.y]))

# Plot coefficients paths
cox_lboost_pathway_fit %>%
  purrr::map(~plot_coef_paths(.x))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lboost_coef_paths"), 6, 4))

## Table of number of selected variables
tbl_cox_lboost_pathway_nselect <- cox_lboost_pathway_fit %>%
  purrr::map(~tbl_n_select_coef(.x))

## Table of of selected coefficients 
tbl_cox_lboost_pathway_coef <- cox_lboost_pathway_fit %>%
  purrr::map(~tbl_select_coef(.x))

## Cox regression with selected variables
clin <- c("age_c_f", "figo_c_f", "hpv_negative_f")
var_select <- cox_lboost_pathway_fit[c("raids_pathway_clin_dna_rna",
                                       "tcga_pathway_clin_dna_rna")]%>%
  purrr::map(~.x$var_select %>%
               subset(!. %in% c("age","figo","hpv_negative")))
tbl_cox_univ_multi_pathway_cox_lboost_select <- list(
  raids$pathway$clin_dna_rna, 
  tcga$pathway$clin_dna_rna
  ) %>%
  set_names(c("raids_pathway_clin_dna_rna",
              "tcga_pathway_clin_dna_rna"))%>%
  purrr::imap(
    ~.x %>%
      dplyr::select(time, event, all_of(clin), all_of(var_select[[.y]]))%>%
      dplyr::mutate(across(starts_with("hallmark"), 
                           ~as.numeric(scale(.x))))%>%
      labelled::set_variable_labels(
        .labels =  var_label(.x[,var_select[[.y]]])
      )
  )%>% 
  purrr::imap(
    ~list(
      tbl_cox(
        data   = .x, 
        time   = "time",
        event  = "event",
        covars = c(clin, var_select[[.y]]),
        model  = "univ"
      )%>%
        gtsummary::add_q(method = "BH") %>%
        gtsummary::bold_p(q = TRUE),
      tbl_cox(
        data   = .x, 
        time   = "time",
        event  = "event",
        covars = c(clin, var_select[[.y]]),
        model  = "multi"
      )%>% 
        gtsummary::add_q(method = "BH") %>%
        gtsummary::bold_p(q = TRUE) %>%
        gtsummary::modify_column_hide(c(stat_n, stat_nevent))
    )%>%
    gtsummary::tbl_merge(tab_spanner = c("**Univariate Cox regression**",
                                         "**Multivariate Cox regression**"))%>%
    gtsummary::modify_spanning_header(c(stat_n_1, stat_nevent_1)~NA)
)

# Model-based boosting ----------------------------------------------------

# cox_mboost_fit <- cohorts %>%
#   purrr::map(~readRDS(here::here("outputs","results",
#                                  paste0(.x,"_cox_mboost_fit.rds"))))%>%
#   set_names(cohorts)#%>%
# purrr::list_flatten()

## Plot selected coefficients 
# cox_mboost_fit %>%
#   purrr::list_flatten()%>%
#   purrr::map(~plot_coef(.x))%>%
#   purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
#                           paste0(.y, "_cox_mboost_coef"), 9, 6))
# 
# cox_mboost_fit %>%
#   purrr::list_flatten()%>%
#   purrr::map(~plot_coef(.x, data_type = T))%>%
#   purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
#                           paste0(.y, "_cox_mboost_coef_omics"), 9, 6))

## Table of selected coefficients
# tbl_cox_mboost <- cox_mboost_fit %>%
#   purrr::map(
#     ~purrr::map(.x, tbl_n_select_coef)%>%
#       .custom_list_tbl_n_select_coef())


# survML: VIM relative to full model (Bio-RAIDs only) ---------------------

# Load results 
vimp_survML_full_est <- readRDS(
  here::here("outputs","results","raids_vimp_survML_full_est.rds")
  )%>%
  purrr::map(~.x %>% dplyr::mutate(data_type = tolower(data_type)))

# Plots of VIMs
vimp_survML_full_est_split <- vimp_survML_full_est %>%
  purrr::map(~.x %>% group_split_custom(landmark_time))%>%
  purrr::list_flatten()

width <- set_names(rep(9, 9), 
                   names(vimp_survML_full_est_split))
height <- set_names(c(rep(18, 3), rep(20, 3), rep(25, 3)), 
                    names(vimp_survML_full_est_split))

vimp_survML_full_est_split %>%
  purrr::imap(~ plot_vimp_est(
    .x,
    ylab   = "",
    type   = "dotplot",
    process_panel = F
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("raids_dotplot_surv_vimp_full_est_", .y),
                           width[.y], height[.y]))

# Summary table of top-10 ranked pathways
tbl_compare_vimp_survML_full_top10 <- vimp_survML_full_est %>%
  purrr::map(
    ~tbl_top10_vimp(
      vims    = .x,
      compare = FALSE,
      n_top   = 10
    )
)

# survML: VIM relative to base model --------------------------------------

## Parallel analysis across cohorts ----

### Main analysis ----

# Load all results 
vimp_survML_base_est <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_vimp_survML_base_est.rds"))))%>%
  set_names(cohorts_name)%>%
  purrr::imap(~.x %>% mutate(cohort = .y))%>%
  bind_rows()%>%
  dplyr::mutate(data_type = tolower(data_type))

# Plots of VIMs
vimp_survML_base_est_split <- vimp_survML_base_est %>%
  group_split_custom(data_type, landmark_time)
width <- set_names(c(20,9,23,9), names(vimp_survML_base_est_split))
height <- set_names(c(15,5,20,5), names(vimp_survML_base_est_split))
process_panel <- set_names(c(T,F,T,F), names(vimp_survML_base_est_split))

vimp_survML_base_est_split %>%
  purrr::imap(~ plot_vimp_est(
    .x,
    ylab   = "",
    type   = "barplot",
    process_panel = process_panel[.y], 
    fill_by = "cohort",
    fill_label = ""
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_barplot_surv_vimp_base_est_", .y,
                                  "_panel"),
                           width[.y], height[.y]))

# for (i in seq_along(vimp_survML_base_est_split)){
#   jpeg(
#     filename = paste0(
#       here::here("docs","articles","briefings_in_bioinformatics"),
#       "/Figure_",1+i,"_revised.jpg"),
#     width = width[i]*300, height = height[i]*300, units = "px", 
#     res = 300, quality = 100
#   )
#   plot_vimp_est(
#     vimp_survML_base_est_split[[i]],
#     ylab   = "",
#     type   = "barplot",
#     process_panel = process_panel[i], 
#     fill_by = "cohort",
#     fill_label = ""
#   )
#   dev.off()
# }

# vimp_survML_base_est_split %>%
#   purrr::imap(~ plot_vimp_est(
#     .x,
#     ylab   = "",
#     type   = "barplot",
#     process_panel = process_panel[.y], 
#     fill_by = "cohort",
#     fill_label = ""
#   )) %>%
#   purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
#                            paste0("compare_barplot_surv_vimp_base_est_", .y,
#                                   "_panel"),
#                            width[.y], height[.y]))
# 
# list.files(here::here("outputs/figures"), 
#            pattern = "^compare_barplot_surv_vimp_base_est_.*panel\\.pdf$", 
#            full.names = T)%>%
#   sort()%>%
#   purrr::iwalk(~{
#     new_name <- file.path(here::here("docs/articles", 
#                                      "computers_in_biology_and_medicine"), 
#                           paste0("Figure_", .y+1, ".pdf"))
#     file.copy(.x, new_name, overwrite = TRUE)
#   })

vimp_survML_base_est_split %>%
  purrr::imap(~ plot_vimp_est(
    .x,
    ylab   = "",
    process_panel = process_panel[.y], 
    type   = "dotplot",
    fill_by = "cohort",
    fill_label = ""
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_dotplot_surv_vimp_base_est_", .y, 
                                  "_panel"),
                           width[.y], height[.y]))

# Summary table of all pathways  
tbl_compare_vimp_survML_base_all <- tbl_vimp(
  vims    = vimp_survML_base_est,
  compare = TRUE
)

# Summary table of top-10 ranked pathways
tbl_compare_vimp_survML_base_top10 <- tbl_top10_vimp(
  vims    = vimp_survML_base_est,
  compare = TRUE,
  n_top   = 10
)

# Univariate Cox regression of top-10 ranked pathways
tbl_cox_univ_survML_base_top10 <- tbl_cox_univ_top10_vimp(
  vims    = vimp_survML_base_est,
  data    = combined, 
  compare = TRUE
)
  
# Cross-cohort comparison 
tbl_vimp_survML_base_top_overlap <- c(5,10,20) %>%
  purrr::map(~tbl_top_overlap_vimp(vims = vimp_survML_base_est, k = .x))%>%
  set_names(paste0("top-", c(5,10,20)))

# Scatter plot with Spearman correlation 
scatter_plot_vimp(vims = vimp_survML_base_est)%>%
  save_plot(here::here("outputs","figures"),
            "compare_scatter_plot_vimp_survML_base_24", 6, 5)

###  Sensitivity analyses ----

####  Adjustment on alteration burden ----

# Load all results 
vimp_survML_base_sens1_adj_est <- cohorts %>%
  purrr::map(~readRDS(here::here(
    "outputs","results",
    paste0(.x,"_vimp_survML_base_sens1_adj_est.rds"))))%>%
  set_names(cohorts_name)%>%
  purrr::imap(~.x %>% mutate(cohort = .y))%>%
  bind_rows()%>%
  dplyr::mutate(data_type = tolower(data_type))

# Summary table of top-10 ranked pathways
tbl_compare_vimp_survML_base_sens1_adj_top10 <- tbl_top10_vimp(
  vims    = vimp_survML_base_sens1_adj_est,
  compare = TRUE,
  n_top   = 10
)

# Cross-cohort comparison 
tbl_vimp_survML_base_sens1_adj_top10_overlap <- tbl_top_overlap_vimp(
  vims = vimp_survML_base_sens1_adj_est, 
  k = 10
)

#### Alternative DNA pathway definitions ----

# At least two altered genes 

# Load all results 
vimp_survML_base_sens2_two_genes_est <- cohorts %>%
  purrr::map(~readRDS(here::here(
    "outputs","results",
    paste0(.x,"_vimp_survML_base_sens2_two_genes_est.rds"))))%>%
  set_names(cohorts_name)%>%
  purrr::imap(~.x %>% mutate(cohort = .y))%>%
  bind_rows()%>%
  dplyr::mutate(data_type = tolower(data_type))

# Summary table of top-10 ranked pathways
tbl_compare_vimp_survML_base_sens2_two_genes_top10 <- tbl_top10_vimp(
  vims    = vimp_survML_base_sens2_two_genes_est,
  compare = TRUE,
  n_top   = 10
)

# Cross-cohort comparison 
tbl_vimp_survML_base_sens2_two_genes_top10_overlap <- tbl_top_overlap_vimp(
  vims = vimp_survML_base_sens2_two_genes_est, 
  k = 10
)

## Pooled analysis ----

# Load all results
combined_vimp_survML_base_est <- readRDS(
  here::here("outputs/results/combined_vimp_survML_base_est.rds"))%>%
  purrr::imap(~.x %>% dplyr::mutate(data_type = tolower(data_type)))%>%
  purrr::set_names(janitor::make_clean_names(names(.)))
  
# Plots of VIMs
combined_vimp_survML_base_est_split <- combined_vimp_survML_base_est %>%
  purrr::map(~.x %>% group_split_custom(data_type, landmark_time))%>%
  purrr::list_flatten()

width <- set_names(rep(c(23,9,23,9),2), 
                   names(combined_vimp_survML_base_est_split))
height <- set_names(rep(c(15,5,20,5),2), 
                    names(combined_vimp_survML_base_est_split))
process_panel <- set_names(rep(c(T,F,T,F),2), 
                           names(combined_vimp_survML_base_est_split))

combined_vimp_survML_base_est_split %>%
  purrr::imap(
    ~ plot_vimp_est(
      .x,
      ylab   = "",
      process_panel = process_panel[.y],
      type   = "barplot"
      )
    )%>%
  purrr::iwalk(
    ~ save_plot(
      .x, here::here("outputs","figures"),
      janitor::make_clean_names(
        paste0("combined_barplot_surv_vimp_base_est_figo_", .y, "_panel")
        ),
      width[.y], height[.y]
      )
    )

list.files(here::here("outputs/figures"), 
           pattern = "^combined_barplot_surv_vimp_base_est_figo_i_ii.*\\.pdf$", 
           full.names = T)%>%
  sort()%>%
  purrr::iwalk(~{
    new_name <- file.path(here::here("docs/articles", 
                                     "computers_in_biology_and_medicine"), 
                          paste0("Supplementary_Figure_S", .y+3, ".pdf"))
    file.copy(.x, new_name, overwrite = TRUE)
  })

combined_vimp_survML_base_est_split %>%
  purrr::imap(
    ~ plot_vimp_est(
      .x,
      ylab   = "",
      process_panel = process_panel[.y],
      type   = "dotplot"
    )
  )%>%
  purrr::iwalk(
    ~ save_plot(
      .x, here::here("outputs","figures"),
      janitor::make_clean_names(
        paste0("combined_dotplot_surv_vimp_base_est_figo_", .y, "_panel")
      ),
      width[.y], height[.y]
    )
  )

# Summary table of all pathways  
tbl_combined_vimp_survML_base_all <- combined_vimp_survML_base_est %>%
  purrr::map(~tbl_vimp(
    vims    = .x,
    compare = FALSE
))%>%
  purrr::list_flatten()

# Summary table of top-10 ranked pathways
tbl_combined_vimp_survML_base_top10 <- combined_vimp_survML_base_est %>%
  purrr::map(~tbl_top10_vimp(
    vims    = .x,
    compare = FALSE,
    n_top   = 10
))%>%
  purrr::list_flatten()

# Univariate Cox regression of top-10 ranked pathways
combined_figo <- combined %>%
  group_split_custom(figo_c_f)%>%
  purrr::set_names(janitor::make_clean_names(names(.)))

tbl_combined_cox_univ_survML_base_top10 <- combined_vimp_survML_base_est %>%
  purrr::imap(~tbl_cox_univ_top10_vimp(
    vims    = .x,
    data    = combined_figo[[.y]], 
    compare = FALSE
))%>%
  purrr::list_flatten()

# ranger: permutation importance ------------------------------------------

# Load all results 
vimp_ranger_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_vimp_ranger_fit.rds"))))%>%
  set_names(cohorts_name)

vimp_ranger_est <- vimp_ranger_fit %>%
  purrr::imap(~data.frame(.x$variable.importance)%>%
                tibble::rownames_to_column()%>%
                set_names(c("variable_name","est"))%>%
                dplyr::mutate(cohort = .y,
                              landmark_time = 0)%>%
                dplyr::left_join(pathways_process %>% 
                                   dplyr::select(data_type, label, variable_name), 
                                 by = c("variable_name"))%>%
                dplyr::mutate(data_type = tolower(data_type),
                              variable = label))%>%
  bind_rows()

# Summary table of all pathways  
tbl_vimp_ranger_est_all <- tbl_vimp(
  vims    = vimp_ranger_est,
  compare = TRUE
)

# Summary table of top-10 ranked pathways
tbl_vimp_ranger_top10 <- tbl_top10_vimp(
  vims    = vimp_ranger_est,
  compare = TRUE,
  n_top   = 10,
  digits  = 4
)

# survex: variable importance with explainable machine learning -----------

# Load all results 
vimp_survex_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_vimp_survex_fit.rds"))))%>%
  set_names(cohorts_name)

# VIMs at all oberved times 
vimp_survex_fit <- vimp_survex_fit %>%
  purrr::map(function(x) {
    x[["model_parts_pathways"]]$result <- x[["model_parts_pathways"]]$result %>%
      dplyr::select(-`_full_model_`)
    var_label(x[["model_parts_pathways"]]$result) <- c(
      labels,
      `_times_`="_times_",
      `_baseline_`="_baseline_",                       
      `_permutation_`="_permutation_",
      label = "label"
    )
    names(x[["model_parts_pathways"]]$result) <- 
      as.character(var_label(x[["model_parts_pathways"]]$result))
    x
  })

## DNA pathways 
cowplot::plot_grid(plotlist = vimp_survex_fit %>%
  purrr::imap(function(x, y) {
    x[["model_parts_pathways"]]$result <- x[["model_parts_pathways"]]$result %>%
      dplyr::select(-all_of(as.character(rna_labels)))
    plot(x[["model_parts_pathways"]], 
         title = y, 
         subtitle = "", 
         max_vars = 10)+
      scale_y_continuous(limits = c(-0.075,0.075))+
      scale_x_continuous(limits = c(0,40))+
      labs(x = "Time in months")
  }), nrow = 2)%>%
  save_plot(here::here("outputs","figures"), "compare_survex_vimp_est_dna_pathways",
            18, 9)

## RNA pathways 
cowplot::plot_grid(plotlist = vimp_survex_fit %>%
                     purrr::imap(function(x, y) {
                       x[["model_parts_pathways"]]$result <- x[["model_parts_pathways"]]$result %>%
                         dplyr::select(-all_of(as.character(dna_labels)))
                       plot(x[["model_parts_pathways"]], 
                            title = y, 
                            subtitle = "", 
                            max_vars = 10)+
                         scale_y_continuous(limits = c(-0.12,0.12))+
                         scale_x_continuous(limits = c(0,40))+
                         labs(x = "Time in months")
                     }), nrow = 2)%>%
  save_plot(here::here("outputs","figures"), "compare_survex_vimp_est_rna_pathways",
            18, 9)

## DNA processes 
cowplot::plot_grid(
  plotlist = vimp_survex_fit %>%
    purrr::imap(function(x, y) {
      x[["model_parts_pathways_groups"]]$result <- 
        x[["model_parts_pathways_groups"]]$result %>%
        dplyr::select(-`_full_model_`,-matches("(RNA)"))
      names(x[["model_parts_pathways_groups"]]$result) <- 
        gsub(' \\(DNA\\)',"", names(x[["model_parts_pathways_groups"]]$result))
      plot(x[["model_parts_pathways_groups"]], 
           title = y,    
           subtitle = "", 
           max_vars = 10)+
        scale_y_continuous(limits = c(-0.075,0.075))+
        scale_x_continuous(limits = c(0,40))+
        labs(x = "Time in months")
      }), nrow = 2)%>%
  save_plot(here::here("outputs","figures"), "compare_survex_vimp_est_dna_processes",
            18, 9)

## RNA processes 
cowplot::plot_grid(
  plotlist = vimp_survex_fit %>%
    purrr::imap(function(x, y) {
      x[["model_parts_pathways_groups"]]$result <- 
        x[["model_parts_pathways_groups"]]$result %>%
        dplyr::select(-`_full_model_`,-matches("(DNA)"))
      names(x[["model_parts_pathways_groups"]]$result) <- 
        gsub(' \\(RNA\\)',"", names(x[["model_parts_pathways_groups"]]$result))
      plot(x[["model_parts_pathways_groups"]], 
           title = y,    
           subtitle = "", 
           max_vars = 10)+
        scale_y_continuous(limits = c(-0.16,0.16))+
        scale_x_continuous(limits = c(0,40))+
        labs(x = "Time in months")
    }), nrow = 2)%>%
  save_plot(here::here("outputs","figures"), "compare_survex_vimp_est_rna_processes",
            18, 9)

# Get all VIMs  
vimp_survex_est <- vimp_survex_fit %>%
  purrr::imap(~get_vimp_est(.x, method = "survex")%>%
                mutate(cohort = .y))%>%
  bind_rows()%>%
  dplyr::mutate(data_type = tolower(data_type))

# Plots of VIMs
vimp_survex_est_split <- vimp_survex_est %>%
  group_split_custom(data_type, landmark_time)
width <- set_names(c(20,9,23,9), names(vimp_survex_est_split))
height <- set_names(c(15,5,20,5), names(vimp_survex_est_split))
process_panel <- set_names(c(T,F,T,F), names(vimp_survex_est_split))
limits <- set_names(list(c(0, 0.004),
                         c(0, 0.007), 
                         c(0, 0.03),
                         c(0, 0.12)), names(vimp_survex_est_split))
vimp_survex_est_split %>%
  purrr::imap(~ plot_vimp_est(
    .x,
    xlab   = "Permutation VIM",
    ylab   = "",
    type   = "barplot",
    method = "survex", 
    process_panel = process_panel[.y], 
    fill_by = "cohort",
    fill_label = "",
    limits = limits[[.y]],
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_barplot_survex_vimp_base_est_", .y,
                                  "_panel"),
                           width[.y], height[.y]))

vimp_survex_est_split %>%
  purrr::imap(~ plot_vimp_est(
    .x,
    ylab   = "",
    process_panel = process_panel[.y], 
    type   = "dotplot",
    method = "survex",
    fill_by = "cohort",
    fill_label = "",
    limits = limits[[.y]]
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_dotplot_survex_vimp_base_est_", .y, 
                                  "_panel"),
                           width[.y], height[.y]))

# Summary table of all pathways  
tbl_vimp_survex_est_all <- tbl_vimp(
  vims    = vimp_survex_est,
  compare = TRUE
)

# Summary table of top-10 ranked pathways
tbl_vimp_survex_top10 <- tbl_top10_vimp(
  vims    = vimp_survex_est,
  compare = TRUE,
  n_top   = 10
)

# Univariate Cox regression of top-10 ranked pathways
tbl_cox_univ_survex_top10 <- tbl_cox_univ_top10_vimp(
  vims    = vimp_survex_est,
  data    = combined, 
  compare = TRUE
)

# Cross-cohort comparison 
tbl_vimp_survex_top_overlap <- c(5,10,20) %>%
  purrr::map(~tbl_top_overlap_vimp(vimp_survex_est, k = .x))%>%
  set_names(paste0("top-", c(5,10,20)))

# Scatter plot with Spearman correlation 
scatter_plot_vimp(vims = vimp_survex_est,
                  annotate_x = 0, 
                  annotate_overall_y = 0.05,
                  annotate_data_type_y = c(0.045, 0.04))%>%
  save_plot(here::here("outputs","figures"),
            "compare_scatter_plot_vimp_survex_24", 5, 4.5)

# ---------------------------- Save all tables --------------------------#

for (i in setdiff(ls(pattern = "tbl_"), lsf.str())){
  saveRDS(get(i), here::here("outputs","tables", paste0(i,".rds")))
}
