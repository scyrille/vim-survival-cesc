
#'@description Compute variable importance estimates with `survML`

library(here)

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))
source(here::here("R","vimp.R"))

# Load data 
combined <- read_processed(cohort = "combined")

# VIM relative to clinical features 

#' Estimating variable importance relative to base model
#' 
#' We consider the importance of each molecular feature relative to a 
#' baseline set of clinical features. The molecular feature of interest 
#' is added to a baseline set of clinical features, with importance 
#' measured by the resulting gain in predictiveness. 

SL.library <- make_SL_library()

combined_input_survML_base_fit <- combined$pathway$clin_dna_rna %>%
  dplyr::filter(figo==1)%>%
  make_input_vimp_survML_base(
    var_clin = c("age","hpv_negative"),
    dna_prefix = "genomic_pathway_",
    rna_prefix = "hallmark_"
  )

start <- Sys.time()
combined_vimp_survML_base_fit <- compute_vimp_survML_base(
    time           = combined_input_survML_base_fit$time,
    event          = combined_input_survML_base_fit$event,
    X              = combined_input_survML_base_fit$X,
    base_features  = combined_input_survML_base_fit$base_features, 
    feature_groups = combined_input_survML_base_fit$feature_groups, 
    SL.library     = SL.library,
    seed           = 123 
    )
end <- Sys.time()
combined_vimp_survML_base_fit_runtime <- 
  as.numeric(difftime(end, start, units = "mins"))

# VIMs  
combined_vimp_survML_base_est <- combined_vimp_survML_base_fit %>%
  get_vimp_survML_est()

# Plots of VIMs
combined_vimp_survML_base_est_split <- combined_vimp_survML_base_est %>%
  group_split_custom(data_type, landmark_time)

width <- set_names(c(12,23,9), names(combined_vimp_survML_base_est_split))
height <- set_names(c(12,20,5), names(combined_vimp_survML_base_est_split))

combined_vimp_survML_base_est_split %>%
  imap(~ plot_vimp_survML_est(
    .x,
    ylab   = "",
    type   = "barplot"
  )) %>%
  iwalk(~ save_plot(.x, here::here("outputs","figures"),
                    paste0("combined_barplot_surv_vimp_base_est_", .y),
                    width[.y], height[.y]))

combined_vimp_survML_base_est_split %>%
  imap(~ plot_vimp_survML_est(
    .x,
    ylab   = "",
    type   = "dotplot"
  )) %>%
  iwalk(~ save_plot(.x, here::here("outputs","figures"),
                    paste0("combined_dotplot_surv_vimp_base_est_", .y),
                    width[.y], height[.y]))

# Tables of VIMs
tbl_combined_vimp_survML_base <- combined_vimp_survML_base_est_split %>%
  map(~tbl_vimp_survML(.x))

# Top-10 ranked pathways 
tbl_combined_vimp_survML_base_top10 <-  combined_vimp_survML_base_est %>%
  # dplyr::filter(!grepl("RNA_processes",data_type))%>%
  group_split_custom(data_type, landmark_time)%>%
  map(~get_top10_vimp_survML(.x))


# ---------------------------- Save all results -------------------------#

for (i in setdiff(ls(pattern = "_fit|_runtime"), ls(pattern = "input"))){
  saveRDS(get(i), here::here("outputs","results", paste0(i,".rds")))
}

# ---------------------------- Save all tables --------------------------#

for (i in setdiff(ls(pattern = "tbl_"), lsf.str())){
  saveRDS(get(i), here::here("outputs","tables", paste0(i,".rds")))
}
