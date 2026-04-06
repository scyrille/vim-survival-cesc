
#'@description Compute variable importance estimates with `survML`

library(here)

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))
source(here::here("R","vimp.R"))

# Load data 
raids <- read_processed(cohort = "raids")

# VIM relative to all features --------------------------------------------

#' Estimating variable importance relative to all features 

#' We consider the importance of each variable relative to a full 
#' feature vector composed of all clinical features and pathways.
#' We will use the exclusion variable importance method to quantify 
#' importance of variables. This procedure involves iteratively 
#' removing each feature from the data set, retraining the model without 
#' that feature and evaluating its impact on model performance. 

input_survML_full_fit <- raids$pathway %>%
  purrr::map(~.x %>% 
               make_input_vimp_survML_full(
                 var_clin = c("age","hpv_negative","figo","necrosis"),
                 dna_prefix = "genomic_pathway_|hrd|tmb_high",
                 rna_prefix = "hallmark_"
                 )
             )

# library(furrr)
# library(future)
# message("Number of parallel workers: ", future::nbrOfWorkers())
# future::plan(multisession, workers = 3)
# message("Number of parallel workers: ", future::nbrOfWorkers())

start <- Sys.time()
vimp_survML_full_fit <- input_survML_full_fit %>%
  # furrr::future_map(
  purrr::map(
    ~ compute_vimp_survML_full(
        time           = .x$time,
        event          = .x$event,
        X              = .x$X,
        feature_groups = .x$feature_groups,
        seed           = 123
        )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("survML", "SuperLearner", "glmnet", "ranger", "xgboost", 
    #                "withr", "dplyr", "stringr", "purrr")
    # )
  )
end <- Sys.time()
vimp_survML_full_fit_runtime <- as.numeric(difftime(end, start, units = "mins"))

# plan(sequential)

# VIM relative to clinical features ---------------------------------------

#' Estimating variable importance relative to base model
#' 
#' We consider the importance of each molecular feature relative to a 
#' baseline set of clinical features. The molecular feature of interest 
#' is added to a baseline set of clinical features, with importance 
#' measured by the resulting gain in predictiveness. 

SL.library <- make_SL_library()

input_survML_base_fit <- raids$pathway$clin_dna_rna %>%
  make_input_vimp_survML_base(
    var_clin = c("age","hpv_negative","figo"),
    dna_prefix = "genomic_pathway_",
    rna_prefix = "hallmark_"
  )

start <- Sys.time()
vimp_survML_base_fit <- compute_vimp_survML_base(
    time           = input_survML_base_fit$time,
    event          = input_survML_base_fit$event,
    X              = input_survML_base_fit$X,
    base_features  = input_survML_base_fit$base_features, 
    feature_groups = input_survML_base_fit$feature_groups, 
    SL.library     = SL.library,
    seed           = 123 
    )
end <- Sys.time()
vimp_survML_base_fit_runtime <- 
  as.numeric(difftime(end, start, units = "mins"))

# ---------------------------- Save all results -------------------------#

for (i in setdiff(ls(pattern = "_fit|_runtime"), ls(pattern = "input"))){
  saveRDS(get(i), here::here("outputs","results", paste0("raids_",i,".rds")))
}
