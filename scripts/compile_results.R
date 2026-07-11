
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

## Bio-RAIDs only (with necrosis, HRD, TMB, MSI, mutationals signatures)
cox_lboost_necrosis_fit <- readRDS(
  here::here("outputs/results/raids_cox_lboost_necrosis_fit.rds")
)

# Plot selected coefficients
cox_lboost_necrosis_fit %>%
  plot_coef()%>%
  save_plot(here::here("outputs","figures"), 
            paste0("raids_cox_lboost_necrosis_clin_dna_rna_coef"), 9, 6)

cox_lboost_necrosis_fit %>%
  plot_coef(data_type = T)%>%
  save_plot(here::here("outputs","figures"),
            paste0("raids_cox_lboost_necrosis_clin_dna_rna_coef_omics"), 9, 6)

# Plot coefficients paths
cox_lboost_necrosis_fit %>%
  plot_coef_paths()%>%
  save_plot(here::here("outputs","figures"), 
            "raids_cox_lboost_necrosis_clin_dna_rna_coef_paths", 10, 6)

# Plot coefficients updates at each boosting step until optimal step
cox_lboost_necrosis_fit %>%
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
  tbl_n_select_coef()

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


# survML: VIM relative to base model ----------------------------------------------

## Parallel analysis across cohorts ----

# Load all results 
vimp_survML_base_fit <- cohorts %>%
  purrr::map(~readRDS(here::here("outputs","results",
                                 paste0(.x,"_vimp_survML_base_fit.rds"))))%>%
  set_names(cohorts_name)

# Get all VIMs  
vimp_survML_base_est <- vimp_survML_base_fit %>%
  purrr::imap(~get_vimp_est(.x)%>%
                mutate(cohort = .y))%>%
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

list.files(here::here("outputs/figures"), 
           pattern = "^compare_barplot_surv_vimp_base_est_.*panel\\.pdf$", 
           full.names = T)%>%
  sort()%>%
  purrr::iwalk(~{
    new_name <- file.path(here::here("docs/articles", 
                                     "computers_in_biology_and_medicine"), 
                          paste0("Figure_", .y+1, ".pdf"))
    file.copy(.x, new_name, overwrite = TRUE)
  })

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
tbl_vimp_survML_base_all <- tbl_vimp(
  vims    = vimp_survML_base_est,
  compare = TRUE
)

# Summary table of top-10 ranked pathways
tbl_vimp_survML_base_top10 <- tbl_top10_vimp(
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
tbl_vimp_survML_top_overlap <- c(5,10,20) %>%
  purrr::map(~tbl_top_overlap_vimp(vims = vimp_survML_base_est, k = .x))%>%
  set_names(paste0("top-", c(5,10,20)))

# Scatter plot with Spearman correlation 
scatter_plot_vimp(vims = vimp_survML_base_est)%>%
  save_plot(here::here("outputs","figures"),
            "compare_scatter_plot_vimp_survML_base_24", 6, 5)

## Pooled analysis ----

# Load all results
combined_vimp_survML_base_fit <- readRDS(
  here::here("outputs/results/combined_vimp_survML_base_fit.rds"))

# VIMs  
combined_vimp_survML_base_est <- combined_vimp_survML_base_fit %>%
  purrr::imap(~get_vimp_est(.x)%>%
                dplyr::mutate(data_type = tolower(data_type)))%>%
  purrr::set_names(janitor::make_clean_names(names(.)))
  
# Plots of VIMs
combined_vimp_survML_base_est_split <- combined_vimp_survML_base_est %>%
  purrr::map(~.x %>%
               group_split_custom(data_type, landmark_time))%>%
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

tbl_cox_univ_survML_base_top10 <- combined_vimp_survML_base_est %>%
  purrr::imap(~tbl_cox_univ_top10_vimp(
    vims    = .x,
    data    = combined_figo[[.y]], 
    compare = FALSE
))%>%
  purrr::list_flatten()

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
      scale_y_continuous(limits = c(-0.2,0.2))+
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
                         scale_y_continuous(limits = c(-0.2,0.2))+
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
        scale_y_continuous(limits = c(-0.2,0.2))+
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
        scale_y_continuous(limits = c(-0.2,0.2))+
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
    limits = c(0, 0.15),
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
    limits = c(0, 0.15)
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
