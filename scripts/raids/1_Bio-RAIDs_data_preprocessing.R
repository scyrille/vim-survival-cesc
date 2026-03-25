
#'@description Data preprocessing

library(here)

source(here::here("R","utils.R"))
source(here::here("R","features_clin.R"))
source(here::here("R","features_genes.R"))
source(here::here("R","features_pathways.R"))

# Import data -------------------------------------------------------------

raids_raw <- import_raids()
saveRDS(raids_raw, here::here("data/raw/raids/raids_raw_data.RDS"))

raids_temp <- raids_raw %>%
  purrr::imap(function(x, y) {
    if (y != "signatures_subgroups_raw") {
      janitor::clean_names(x)
    } else {
      x
    }
  })

# Prepare clinical data  --------------------------------------------------

raids_clin <- raids_temp %>%
  {
    clin <- purrr::pluck(., "clin_raw") %>%
      dplyr::select(patient_id:hpv_status_ngs) %>%
      dplyr::distinct(patient_id, .keep_all = TRUE)
    
    noninfo <- setdiff(names(clin), 
                       c("patient_id", "hpv_status", "hpv_status_ngs"))
    
    hpv <- purrr::pluck(., "hpv_raw")
    
    necrosis <- purrr::pluck(., "necrosis_raw")%>%
      group_by(patient_id)%>%
      dplyr::filter(row_number()==1)
    
    samples <- purrr::pluck(., "patient_sample_seq_raw")
    
    clin %>%
      dplyr::filter(!if_all(dplyr::all_of(noninfo), is.na))%>%
      left_join(hpv, by = "patient_id")%>%
      left_join(necrosis, by = "patient_id")%>%
      left_join(samples, by = "patient_id")
  }%>%
  dplyr::select(patient_id, wes_tumor, s_wgs_tumor, rn_aseq_tumor, 
                delai_pfs, statut_pfs, 
                age, figo_18_4_classes, hpv, type_histo_4_classes, 
                necrosis_in_percent_of_total_area)

raids_clin_stand <- standardize_clin(
  df             = raids_clin,
  cohort_name    = "Bio-RAIDs",
  id_col         = "patient_id", 
  wes_id_col     = "wes_tumor", 
  wgs_id_col     = "s_wgs_tumor", 
  rna_id_col     = "rn_aseq_tumor",
  time_col       = "delai_pfs", 
  event_col      = "statut_pfs", 
  age_col        = "age", 
  figo_col       = "figo_18_4_classes", 
  hpv_col        = "hpv", 
  histo_col      = "type_histo_4_classes",
  necrosis_col   = "necrosis_in_percent_of_total_area")

# Prepare molecular data --------------------------------------------------

## Genomic data ----

### Gene-level and pathway-level data 

raids_dna <- raids_temp %>%
  pivot_wider_dna(mut_key             = "mut_wes_raw",
                  cohort_name         = "Bio-RAIDs",
                  gene_col            = "gene",
                  id_col              = "patient_id",
                  altered_col         = "altered_gene",
                  cnv_col             = "cnv",
                  snv_col             = "snv", 
                  fusion_col          = "fusion",
                  pathway_col         = "pathway",
                  altered_pathway_col = "altered_pathway")

## Transcriptomic data ----

### Gene-level data 

raids_rna_prep_gene <- raids_raw$rna_seq_coding_raw %>%
  dplyr::rename(., gene_id = X)%>%
  rename_with(
    ~raids_clin$patient_id[match(.x, raids_clin$rn_aseq_tumor)], 
    .cols = intersect(raids_clin$rn_aseq_tumor, names(.)))%>%
  dplyr::inner_join(raids_raw$rna_seq_gene_raw %>% 
                      dplyr::select(gene_name, gene_id), 
                    by = "gene_id")%>% # retrieving genes name
  group_by(gene_name)%>%
  mutate(gene_name_old = gene_name, 
         gene_name = ifelse(grepl("PAR_Y", gene_id), 
                            paste0(gene_name_old, "_PAR_Y"), 
                            gene_name_old))%>%
  mutate(gene_name = ifelse(row_number()==2, 
                            paste0(gene_name, "_readthrough"), 
                            gene_name))%>%
  dplyr::select(-c(gene_id, gene_name_old))%>%
  pivot_longer(-gene_name, 
               names_to = "patient_id", 
               values_to = "tpm")

raids_rna_gene <- raids_rna_prep_gene %>%
  pivot_wider_rna(cohort_name = "Bio-RAIDs",
                  gene_col    = "gene_name",
                  id_col      = "patient_id",
                  rna_tpm_col = "tpm")

### Pathway-level data 

raids_rna_prep_gsva <- left_join(raids_raw$rna_seq_coding_raw %>%
                                   dplyr::rename(., gene_id = X)%>%
                                   relocate(sort(names(.))), 
                                 raids_raw$rna_seq_gene_raw %>%
                                   dplyr::select(gene_id, gene_name), 
                                 by = "gene_id")%>%
  dplyr::relocate(gene_id, gene_name)%>%
  rename_with(
    ~raids_clin$patient_id[match(.x, raids_clin$rn_aseq_tumor)], 
    .cols = intersect(raids_clin$rn_aseq_tumor, names(.)))

#' Gene Set Variation Analysis (GSVA)
#'The idea here is that we will get pathway-level scores for each sample that 
#'indicate if genes in a pathway vary concordantly in one direction 
#'(over-expressed or under-expressed relative to the overall population)
#'[@Hänzelmann 2013a]. This means that GSVA scores will depend on the 
#'samples included in the dataset when you run GSVA; if you added more samples 
#'and ran GSVA again, you would expect the scores to change [@Hänzelmann 2013a].

raids_gsva <- list(
  clin_rna  = c("clin", "rna"),
  clin_dna_rna = c("clin", "dna", "rna")
) %>%
  purrr::imap(\(req, nm)
              raids_rna_prep_gsva %>%
                dplyr::select(gene_name, gene_id,
                              dplyr::any_of(raids_clin_stand %>% 
                                              drop_na(dplyr::all_of(req)) %>% 
                                              dplyr::pull(patient_id))
                ) %>%
                compute_gsva(cohort_name = "Bio-RAIDs")
  )

# Data integration  -------------------------------------------------------

raids <- integrate_data(
  id_col             = c("patient_id","cohort"),
  clin               = raids_clin_stand,
  dna                = raids_dna %>% drop_na(),
  rna_gene           = raids_rna_gene,
  gsva_clin_rna      = raids_gsva$clin_rna,
  gsva_clin_dna_rna  = raids_gsva$clin_dna_rna,
  keep_only_complete = TRUE, 
  pre_filter         = TRUE
)

saveRDS(raids, here::here("data/processed/raids/raids_processed_data.RDS"))
