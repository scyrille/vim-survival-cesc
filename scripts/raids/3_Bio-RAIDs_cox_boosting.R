
#'@description Cox regression with model-based and likelihood-based boosting 

library(here)

source(here::here("R","utils.R"))
source(here::here("R","cox_boosting.R"))

# Load data 
raids <- read_processed(cohort = "raids")%>%
  purrr::list_flatten()%>% 
  # Pathway and gene-level analyses
  { .[setdiff(names(.), c("clin",
                          "gene_clin_dna",
                          "gene_clin_rna"))]}

clin <- c("age","figo","hpv_negative")
dna_pattern <- "altered_|genomic_pathway_"
rna_pattern <- "rna_seq_|hallmark_"

X_names <- raids %>%
  { .[setdiff(names(.), "gene_clin_dna_rna")]}%>%
  purrr::map(~.x %>% 
               dplyr::select(all_of(clin), 
                             tidyselect::matches(dna_pattern),
                             tidyselect::matches(rna_pattern),
                             -dplyr::ends_with(c("two_genes","count","prop",
                                                 "alteration_burden"))) %>% 
               names())

# Model-based boosting ----------------------------------------------------

## Pathway-level
cox_mboost_fit <- raids %>%
  { .[setdiff(names(.), "gene_clin_dna_rna")]}%>%
  purrr::imap(
    ~ fit_cox_mboost(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      seed  = 123
    )
  )

# Likelihood-based boosting -----------------------------------------------

## Gene-level 
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

## Pathway-level
cox_lboost_fit <- raids %>%
  { .[setdiff(names(.), "gene_clin_dna_rna")]}%>%
  purrr::imap(
    ~ fit_cox_lboost(
      time      = .x$time,
      event     = .x$event,
      X         = .x %>% dplyr::select(all_of(X_names[[.y]])),
      mandatory = c("age", "figo", "hpv_negative"),
      seed      = 123, 
      stepno    = 30
    )
  )

# ---------------------------- Save all results -------------------------#

for (i in ls(pattern = "_fit")){
  saveRDS(get(i), here::here("outputs","results", paste0("raids_",i,".rds")))
}
