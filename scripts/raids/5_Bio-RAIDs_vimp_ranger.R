
#'@description Compute variable importance estimates with `ranger` and `survex`

library(here)
library(survex)
library(survival)
library(ranger)

source(here::here("R","utils.R"))
source(here::here("R","features_pathways.R"))

# Load data 
raids <- read_processed(cohort = "raids")

X_pathways <- raids$pathway$clin_dna_rna %>%
  dplyr::select(starts_with("genomic_pathway"),
                starts_with("hallmark_"))%>%
  colnames()

raids_vimp <- raids$pathway$clin_dna_rna %>%
  dplyr::select(time, event, all_of(X_pathways))

# ranger ------------------------------------------------------------------

vimp_ranger_fit <- ranger::ranger(Surv(time, event)~., raids_vimp,
                                  importance = "permutation", 
                                  splitrule = "logrank",
                                  seed = 123)

saveRDS(vimp_ranger_fit, 
        here::here("outputs","results", "raids_vimp_ranger_fit.rds"))

# survex ------------------------------------------------------------------

horizon_time <- 24 
times <- sort(unique(c(0, raids_vimp$time, horizon_time)))
# times <- times[times <= max(horizon_time)]

# Pathways grouped by processes
group_DNA_pathway_features <- get_pathways_process() %>%
  dplyr::filter(data_type == "DNA_pathways") %>%
  dplyr::select(process, variable_name) %>%
  dplyr::transmute(process = paste0(process, " (DNA)"), variable_name) %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(variable_name = list(variable_name), .groups = "drop") %>%
  tibble::deframe()

group_RNA_pathway_features <- get_pathways_process() %>%
  dplyr::filter(data_type == "RNA_pathways") %>%
  dplyr::transmute(process = paste0(process, " (RNA)"), variable_name) %>%
  dplyr::group_by(process) %>%
  dplyr::summarise(variable_name = list(variable_name), .groups = "drop") %>%
  tibble::deframe()

# Fit random survival forests
ranger_fit <- ranger::ranger(Surv(time, event)~., raids_vimp,
                             importance = "none", 
                             splitrule = "logrank",
                             seed = 123)

# Model performance
ranger_exp <- survex::explain(ranger_fit, data = raids_vimp,
                              y = Surv(raids_vimp$time, raids_vimp$event),
                              times = times, seed = 123)
ranger_perf <- survex::model_performance(
  ranger_exp,
  metrics = c(`C/D AUC` = cd_auc),
  times = horizon_time
)
ranger_perf$result

# Variable importance of all pathways
ranger_vimp_pathways <-  survex::model_parts(
  ranger_exp, 
  B = 10,
  type = "difference",
  variables = X_pathways,
  loss_function = survex::loss_one_minus_cd_auc
)

# Variable importance of all pathway groups 
ranger_vimp_pathways_groups <-  survex::model_parts(
  ranger_exp, 
  B = 10,
  type = "difference",
  variable_groups = c(group_DNA_pathway_features,
                      group_RNA_pathway_features),
  loss_function = survex::loss_one_minus_cd_auc
)

vimp_ranger_survex_fit <- list(
  fit                         = ranger_fit, 
  exp                         = ranger_exp, 
  perf                        = ranger_perf, 
  model_parts_pathways        = ranger_vimp_pathways,
  model_parts_pathways_groups = ranger_vimp_pathways_groups
)

saveRDS(vimp_ranger_survex_fit, 
        here::here("outputs","results", "raids_vimp_survex_fit.rds"))

rm(list = ls(pattern = "fit"))
