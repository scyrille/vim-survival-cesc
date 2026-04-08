
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
tcga <- read_processed(cohort = "tcga")%>%
  purrr::list_flatten()%>%
  { .[setdiff(names(.), c("clin",
                          "gene_clin_dna",
                          "gene_clin_rna",
                          "gene_clin_dna_rna"
                          ))]}

clin <- c("age","figo","hpv_negative")
dna_pattern <- "genomic_pathway_"
rna_pattern <- "hallmark_"

X_names <- tcga %>% 
  purrr::map(~.x %>% 
               dplyr::select(all_of(clin), 
                             tidyselect::matches(dna_pattern),
                             tidyselect::matches(rna_pattern)) %>% 
               names())

# Model-based boosting ----------------------------------------------------

cox_mboost_fit <- tcga %>%
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

cox_lboost_fit <- tcga %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_lboost(
      time      = .x$time,
      event     = .x$event,
      X         = .x %>% dplyr::select(all_of(X_names[[.y]])),
      mandatory = c("age", "figo", "hpv_negative"),
      stepno    = 30,    
      seed      = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   globals  = TRUE,
    #   packages = c("CoxBoost", "dplyr", "tidyselect", "labelled")
    # )
  )

# plan(sequential)

# ---------------------------- Save all results -------------------------#

for (i in ls(pattern = "_fit")){
  saveRDS(get(i), here::here("outputs","results", paste0("tcga_",i,".rds")))
}
