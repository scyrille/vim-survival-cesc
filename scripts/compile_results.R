
#'@description Compile results across cohorts 

library(here)

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))
source(here::here("R","cox_lasso.R"))
source(here::here("R","cox_boosting.R"))
source(here::here("R","vimp.R"))

cohorts <-  c("raids","tcga")
cohorts_name <- c("Bio-RAIDs","TCGA-CESC")

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

# Cox sparse group lasso --------------------------------------------------

cox_sparse_group_lasso_fit <- cohorts %>%
  purrr::map(~readRDS(
    here::here("outputs","results",
               paste0(.x,"_cox_sparse_group_lasso_fit.rds"))))#%>%
  # purrr::list_flatten()

## Table of selected coefficients
tbl_cox_sparse_group_lasso <- cox_sparse_group_lasso_fit %>%
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

## Bio-RAIDs only (with necrosis, HRD, TMB, MSI, mutationals signatures)
cox_lboost_necrosis_fit <- readRDS(
  here::here("outputs/results/raids_cox_lboost_necrosis_fit.rds")
)

# Plot selected coefficients
cox_lboost_necrosis_fit %>%
  purrr::map(~plot_coef(.x))%>%
  purrr::iwalk(
    ~save_plot(.x, here::here("outputs","figures"), 
               paste0("raids_cox_lboost_necrosis_", .y, "_coef"), 9, 6))

cox_lboost_necrosis_fit %>%
  purrr::map(~plot_coef(.x, data_type = T))%>%
  purrr::iwalk(
    ~save_plot(.x, here::here("outputs","figures"),
               paste0("raids_cox_lboost_necrosis_", .y ,"_coef_omics"), 9, 6))

# Plot coefficients paths
cox_lboost_necrosis_fit$clin_dna_rna %>%
  plot_coef_paths()%>%
  save_plot(here::here("outputs","figures"), 
            "raids_cox_lboost_necrosis_clin_dna_rna_coef_paths", 10, 6)

# Plot coefficients updates at each boosting step until optimal step
cox_lboost_necrosis_fit$clin_dna_rna %>%
  (\(fit_obj) 0:fit_obj$opt_step %>%
     purrr::walk(
       ~ plot_coef_updates(fit_obj, step_show = .x) %>%
         save_plot(
           here::here("outputs","figures"),
           paste0("raids_cox_lboost_necrosis_clin_dna_rna_coef_updates_", .x),
           9, 6)
     )
  )()

# Table of selected coefficients
tbl_cox_lboost_necrosis <- cox_lboost_necrosis_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

## Across both cohorts 
cox_lboost_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_lboost_fit.rds"))))%>%
  set_names(cohorts)#%>%
# purrr::list_flatten()

## Plot selected coefficients 
cox_lboost_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lboost_coef"), 9, 6))

cox_lboost_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x, data_type = T))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_lboost_coef_omics"), 9, 6))

## Table of selected coefficients
tbl_cox_lboost <- cox_lboost_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# Model-based boosting ----------------------------------------------------

cox_mboost_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_cox_mboost_fit.rds"))))%>%
  set_names(cohorts)#%>%
# purrr::list_flatten()

## Plot selected coefficients 
cox_mboost_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_mboost_coef"), 9, 6))

cox_mboost_fit %>%
  purrr::list_flatten()%>%
  purrr::map(~plot_coef(.x, data_type = T))%>%
  purrr::iwalk(~save_plot(.x, here::here("outputs","figures"), 
                          paste0(.y, "_cox_mboost_coef_omics"), 9, 6))

## Table of selected coefficients
tbl_cox_mboost <- cox_mboost_fit %>%
  purrr::map(
    ~purrr::map(.x, tbl_n_select_coef)%>%
      .custom_list_tbl_n_select_coef())

# VIM relative to all features --------------------------------------------

## Bio-RAIDs only 


# VIM relative to clinical features ---------------------------------------

# Load all results 
vimp_survML_base_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_vimp_survML_base_fit.rds"))))%>%
  set_names(cohorts_name)

# VIMs  
vimp_survML_base_est <- vimp_survML_base_fit %>%
  purrr::imap(~get_vimp_survML_est(.x)%>%
                mutate(cohort = .y))

# Plots of VIMs
vimp_survML_base_est_bind <- vimp_survML_base_est %>%
  bind_rows()%>%
  group_split_custom(data_type, landmark_time)
# width <- set_names(rep(10, 9), names(vimp_survML_base_est))
# height <- set_names(c(rep(10,3), rep(20,3), rep(6,3)),
#                     names(vimp_survML_base_est))
width <- set_names(c(12,23,9), names(vimp_survML_base_est_bind))
height <- set_names(c(12,20,5), names(vimp_survML_base_est_bind))
# ylab <- set_names(names(vimp_survML_base_est_bind) %>%
#   stringr::str_replace_all(c("DNA_"="DNA-based ",
#                              "RNA_"="RNA-based "))%>%
#   gsub(paste0("_",c(12,24,36), collapse = "|"), "", .), 
#   names(vimp_survML_base_est_bind))

vimp_survML_base_est_bind %>%
  purrr::imap(~ plot_vimp_survML_est(
    .x,
    # xlab = "", 
    ylab   = "",
    type   = "barplot",
    fill_by = "cohort",
    fill_label = ""
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_barplot_surv_vimp_base_est_", .y),
                           width[.y], height[.y]))

vimp_survML_base_est_bind %>%
  purrr::imap(~ plot_vimp_survML_est(
    .x,
    ylab   = "",
    type   = "dotplot",
    fill_by = "cohort",
    fill_label = ""
  )) %>%
  purrr::iwalk(~ save_plot(.x, here::here("outputs","figures"),
                           paste0("compare_dotplot_surv_vimp_base_est_", .y),
                           width[.y], height[.y]))

# Tables of VIMs
tbl_vimp_survML_base <- vimp_survML_base_est %>%
  purrr::map(~.x %>% 
        group_split_custom(data_type, landmark_time)%>%
        purrr::map(~tbl_vimp_survML(.x)))

# Scatter plot with Spearman correlation 
vimp_survML_base_est %>%
  bind_rows()%>%
  scatter_plot_vimp()%>%
  save_plot(here::here("outputs","figures"),
            "compare_scatter_plot_vimp_survML_base_24", 6, 5)

# Top-10 ranked pathways 
tbl_vimp_survML_base_top10 <-  vimp_survML_base_est %>%
  map(~.x %>%
        # dplyr::filter(!grepl("RNA_processes", data_type))%>%
        group_split_custom(data_type, landmark_time)%>%
        map(~get_top10_vimp_survML(.x)))

# Overlap of top-10 ranked pathways 
tbl_vimp_survML_base_top10_overlap <- vimp_survML_base_est %>%
  get_overlap_top10_vimp_survML()


# ---------------------------- Save all tables --------------------------#

for (i in setdiff(ls(pattern = "tbl_"), lsf.str())){
  saveRDS(get(i), here::here("outputs","tables", paste0(i,".rds")))
}
