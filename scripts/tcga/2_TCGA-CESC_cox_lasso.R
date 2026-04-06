
#'@description Cox regression with lasso penalties and likelihood-based boosting 

library(here)

source(here::here("R","utils.R"))
source(here::here("R","cox_lasso.R"))

# library(furrr)
# library(future)
# message("Number of parallel workers: ", future::nbrOfWorkers())
# future::plan(multisession, workers = 3)
# message("Number of parallel workers: ", future::nbrOfWorkers())

# Load data 
tcga <- read_processed(cohort = "tcga")%>%
  purrr::list_flatten()%>%
  { .[setdiff(names(.), "clin")]}

# Feature groups and block structures -------------------------------------

#' Variables were grouped into clinical, genomic, and transcriptomic 
#' modalities based on naming conventions. Corresponding group labels and 
#' numeric indices were defined for use in group-based penalized models. 
#' Block structures were then generated for different combinations of data 
#' types, including all possible block orderings, to support priority 
#' lasso modeling.

clin <- c("age","figo","hpv_negative")
dna_pattern <- "altered_|genomic_pathway_"
rna_pattern <- "rna_seq_|hallmark_"

X_names <- tcga %>% 
  purrr::map(~.x %>% 
               dplyr::select(all_of(clin), 
                             tidyselect::matches(dna_pattern),
                             tidyselect::matches(rna_pattern)) %>% 
               names())

# Define clinical, DNA and RNA groups 
group <- purrr::imap(
  X_names,
  ~ {as.factor(dplyr::case_when(
    .x %in% clin ~ "Clinical",
    stringr::str_starts(.x, dna_pattern) ~ "DNA",
    stringr::str_starts(.x, rna_pattern)  ~ "RNA",
    TRUE ~ NA_character_))
  })

index <- purrr::map(
  X_names,
  ~ {
    has_dna <- any(stringr::str_starts(.x, dna_pattern))
    
    dplyr::case_when(
      .x %in% clin ~ 1,
      stringr::str_starts(.x, dna_pattern) ~ 2,
      stringr::str_starts(.x, rna_pattern) ~ if (has_dna) 3 else 2,
      TRUE ~ NA_real_
    )
  }
)

blocks_list <- purrr::map(
  X_names,
  ~ {
    index_list <- purrr::compact(list(
      index_clin = which(.x %in% clin),
      index_dna  = which(stringr::str_starts(.x, dna_pattern)),
      index_rna  = which(stringr::str_starts(.x, rna_pattern))
    ))
    
    do.call(.make_blocks_priority_lasso, unname(index_list))
  }
)

blocks_unpen_list <- purrr::map(
  X_names,
  ~ {
    index_clin <- which(.x %in% clin)
    index_dna  <- which(stringr::str_starts(.x, dna_pattern))
    index_rna  <- which(stringr::str_starts(.x, rna_pattern))
    
    has_dna <- length(index_dna) > 0
    has_rna <- length(index_rna) > 0
    
    if (has_dna & has_rna) {
      list(
        list(bp1 = index_clin, bp2 = index_dna, bp3 = index_rna),
        list(bp1 = index_clin, bp2 = index_rna, bp3 = index_dna)
      )
    } else if (has_dna) {
      list(list(bp1 = index_clin, bp2 = index_dna))
    } else if (has_rna) {
      list(list(bp1 = index_clin, bp2 = index_rna))
    }
  }
) 

# Lasso -------------------------------------------------------------------

cox_lasso_fit <- tcga %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_lasso(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      seed  = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("glmnet", "dplyr", "tidyselect", "labelled")
    # )
  )

# Group lasso -------------------------------------------------------------

cox_group_lasso_fit <- tcga %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_group_lasso(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      group = group[[.y]], 
      seed  = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("grpreg", "dplyr", "tidyselect", "labelled")
    # )
  )

# Sparse group lasso ------------------------------------------------------

# cox_sparse_group_lasso_fit <- tcga %>%
#   # furrr::future_imap(
#   purrr::imap(
#     ~ fit_cox_sparse_group_lasso(
#       time  = .x$time,
#       event = .x$event,
#       X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
#       index = index[[.y]],
#       seed  = 123
#     )#,
#     # .options = furrr::furrr_options(
#     #   seed     = TRUE,
#     #   packages = c("SGL", "dplyr", "tidyselect", "labelled")
#     # )
#   )

# Priority lasso ----------------------------------------------------------

cox_priority_lasso_fit <- tcga %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_priority_lasso(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      blocks_list = blocks_list[[.y]], 
      block1_penalization = TRUE, 
      seed  = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("prioritylasso", "dplyr", "tidyselect", "labelled")
    # )
  )

cox_priority_lasso_unpen_fit <- tcga %>%
  # furrr::future_imap(
  purrr::imap(
    ~ fit_cox_priority_lasso(
      time  = .x$time,
      event = .x$event,
      X     = .x %>% dplyr::select(all_of(X_names[[.y]])),
      blocks_list = blocks_unpen_list[[.y]], 
      block1_penalization = FALSE, 
      seed  = 123
    )#,
    # .options = furrr::furrr_options(
    #   seed     = TRUE,
    #   packages = c("prioritylasso", "dplyr", "tidyselect", "labelled")
    # )
  )

# plan(sequential)

# ---------------------------- Save all results -------------------------#

for (i in ls(pattern = "_fit")){
  saveRDS(get(i), here::here("outputs","results", paste0("tcga_",i,".rds")),
          compress = FALSE)
}
