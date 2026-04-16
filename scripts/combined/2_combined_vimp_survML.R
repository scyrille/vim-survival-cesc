
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

SL.library <- make_SL_library() %>% 
  subset(.!="ridge_1") # remove ridge 

combined_input_survML_base_fit <- combined$pathway$clin_dna_rna %>%
  group_split_custom(figo_c_f)%>%
  purrr::map(~make_input_vimp_survML_base(
    data       = .x, 
    var_clin   = c("age","hpv_negative"),
    dna_prefix = "genomic_pathway_",
    rna_prefix = "hallmark_"
    )
  )

start <- Sys.time()
combined_vimp_survML_base_fit <- combined_input_survML_base_fit %>%
  purrr::map(~compute_vimp_survML_base(
    time           = .x$time,
    event          = .x$event,
    X              = .x$X,
    base_features  = .x$base_features, 
    feature_groups = .x$feature_groups, 
    SL.library     = SL.library,
    seed           = 123 
    )
  )
end <- Sys.time()
combined_vimp_survML_base_fit_runtime <- 
  as.numeric(difftime(end, start, units = "mins"))

# ---------------------------- Save all results -------------------------#

for (i in setdiff(ls(pattern = "_fit|_runtime"), ls(pattern = "input"))){
  saveRDS(get(i), here::here("outputs","results", paste0(i,".rds")))
}
