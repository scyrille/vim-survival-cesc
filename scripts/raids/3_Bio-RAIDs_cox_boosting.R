
#'@description Cox regression with model-based and likelihood-based boosting 

library(here)

source(here::here("R","utils.R"))
source(here::here("R","cox_boosting.R"))

# library(furrr)
# library(future)
# message("Number of parallel workers: ", future::nbrOfWorkers())
# future::plan(multisession, workers = 3)
# message("Number of parallel workers: ", future::nbrOfWorkers())

# Load data 
raids <- read_processed(cohort = "raids")%>%
  purrr::list_flatten()%>%
  { .[setdiff(names(.), "clin")]}

clin <- c("age","figo","hpv_negative")
dna_pattern <- "altered_|genomic_pathway_"
rna_pattern <- "rna_seq_|hallmark_"

X_names <- raids %>% 
  purrr::map(~.x %>% 
               dplyr::select(all_of(clin), 
                             tidyselect::matches(dna_pattern),
                             tidyselect::matches(rna_pattern)) %>% 
               names())

# Model-based boosting ----------------------------------------------------

cox_mboost_fit <- raids %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_mboost(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      seed  = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("mboost", "dplyr", "tidyselect", "labelled")
    # )
  )

# Likelihood-based boosting -----------------------------------------------

cox_lboost_fit <- raids %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_lboost(
      time      = .x$time,
      event     = .x$event,
      X         = .x %>% dplyr::select(all_of(X_names[[.y]])),
      mandatory = c("age", "figo", "hpv_negative"),
      seed      = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   globals  = TRUE,
    #   packages = c("CoxBoost", "dplyr", "tidyselect", "labelled")
    # )
  )

## Model with necrosis, HRD, TMB, MSI, mutational signatures (for Bio-RAIDs only)
cox_lboost_necrosis_fit <- fit_cox_lboost(
  time      = raids$gene_clin_dna_rna$time,
  event     = raids$gene_clin_dna_rna$event,
  X         = raids$gene_clin_dna_rna %>% dplyr::select(
    age, figo, hpv_negative, necrosis,
    dplyr::any_of(c("hrd", "tmb_high", "msi_high")),
    tidyselect::starts_with("altered_"),
    tidyselect::starts_with("rna_seq_"),
    tidyselect::starts_with("sig_group_")),
  mandatory = c("age", "figo", "hpv_negative", "necrosis"),
  seed      = 123
)

# plan(sequential)

# ---------------------------- Save all results -------------------------#

for (i in ls(pattern = "_fit")){
  saveRDS(get(i), here::here("outputs","results", paste0("raids_",i,".rds")))
}
