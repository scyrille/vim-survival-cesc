
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
  purrr::map(
    ~.x %>% 
      make_input_vimp_survML_full(
        var_clin = c("age","hpv_negative","figo","necrosis","hrd","tmb_high"),
        dna_prefix = "genomic_pathway_",
        dna_suffix = "any_gene", 
        rna_prefix = "hallmark_"
      )
    )

vimp_survML_full_fit <- input_survML_full_fit %>%
  purrr::map(
    ~ compute_vimp_survML_full(
        time           = .x$time,
        event          = .x$event,
        X              = .x$X,
        feature_groups = .x$feature_groups
        )
    )

vimp_survML_full_est <- vimp_survML_full_fit %>%
  purrr::map(
    ~get_vimp_est(
      fit           = .x, 
      landmark_time = c(12,24,36), 
      method        = "survML"
  )
)

# VIM relative to clinical features ---------------------------------------

## Main analysis ----

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
    dna_suffix = "any_gene", 
    rna_prefix = "hallmark_"
  )

vimp_survML_base_fit <- compute_vimp_survML_base(
    time           = input_survML_base_fit$time,
    event          = input_survML_base_fit$event,
    X              = input_survML_base_fit$X,
    base_features  = input_survML_base_fit$base_features, 
    feature_groups = input_survML_base_fit$feature_groups, 
    SL.library     = SL.library
    )

vimp_survML_base_est <- get_vimp_est(
    fit = vimp_survML_base_fit, 
    landmark_time = 24, 
    method = "survML"
)

## Sensitivity analyses ----

SL.library <- make_SL_library() 

### Adjustment on genomic alteration burden ----
input_survML_base_sens1_adj_fit <- raids$pathway$clin_dna_rna %>%
  make_input_vimp_survML_base(
    var_clin = c("age","hpv_negative","figo",
                 "genomic_pathway_alteration_burden"),
    dna_prefix = "genomic_pathway_",
    dna_suffix = "any_gene", 
    rna_prefix = "hallmark"
  )

vimp_survML_base_sens1_adj_fit <- compute_vimp_survML_base(
  time           = input_survML_base_sens1_adj_fit$time,
  event          = input_survML_base_sens1_adj_fit$event,
  X              = input_survML_base_sens1_adj_fit$X,
  base_features  = input_survML_base_sens1_adj_fit$base_features, 
  feature_groups = input_survML_base_sens1_adj_fit$feature_groups, 
  SL.library     = SL.library
)

vimp_survML_base_sens1_adj_est <- get_vimp_est(
  fit = vimp_survML_base_sens1_adj_fit, 
  landmark_time = 24, 
  method = "survML"
)

### Alternative DNA pathway definition ----

# At least two altered genes 

# Remove non-altered pathway features 
dna_non_altered_pathways <- raids$pathway$clin_dna_rna %>%
  dplyr::select(dplyr::starts_with("genomic_pathway") &
                  ends_with("two_genes"))%>%
  dplyr::summarise(
    dplyr::across(everything(),
                  ~ sum(. == 1, na.rm = TRUE) / dplyr::n())
  ) %>%
  dplyr::select(where(~ .x == 0))%>%
  colnames()

input_survML_base_sens2_two_genes_fit <- raids$pathway$clin_dna_rna %>%
  dplyr::select(-starts_with(dna_non_altered_pathways))%>%
  make_input_vimp_survML_base(
    var_clin   = c("age", "hpv_negative", "figo"),
    dna_prefix = "genomic_pathway",
    dna_suffix = "two_genes",
    rna_prefix = "hallmark"
  )%>%
  purrr::modify_in(
    "feature_groups",
    \(x) {
      remove <- endsWith(names(x), "(DNA process)") |
        endsWith(names(x), "(RNA process)") |
        endsWith(names(x), "(RNA pathway)")
      
      x[!remove]
    }
  )

vimp_survML_base_sens2_two_genes_fit <- compute_vimp_survML_base(
  time           = input_survML_base_sens2_two_genes_fit$time,
  event          = input_survML_base_sens2_two_genes_fit$event,
  X              = input_survML_base_sens2_two_genes_fit$X,
  base_features  = input_survML_base_sens2_two_genes_fit$base_features, 
  feature_groups = input_survML_base_sens2_two_genes_fit$feature_groups, 
  SL.library     = SL.library
)

vimp_survML_base_sens2_two_genes_est <- get_vimp_est(
  fit = vimp_survML_base_sens2_two_genes_fit, 
  landmark_time = 24, 
  method = "survML"
)

# ---------------------------- Save all results -------------------------#

for (i in setdiff(ls(pattern = "_est|_fit"), ls(pattern = "input"))){
  saveRDS(get(i), here::here("outputs","results", paste0("raids_",i,".rds")))
}
