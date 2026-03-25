
#'@description Data preprocessing

library(here)

source(here::here("R","utils.R"))
source(here::here("R","tcga_cesc.R"))
source(here::here("R","features_clin.R"))
source(here::here("R","features_genes.R"))
source(here::here("R","features_pathways.R"))

# Import data -------------------------------------------------------------

## Download TCGA-CESC data (optional) 
# tcga_raw <- import_tcga(raw_dir = "data/raw/tcga", download = TRUE)
# saveRDS(tcga_raw, here::here("data/raw/tcga/tcga_raw_data.RDS"))

# Prepare clinical data  --------------------------------------------------

tcga_clin <- tcga_raw %>%
  map(~.x %>% janitor::clean_names()) %>%
  {
    patient <- .$clin$patient %>%
      dplyr::select(bcr_patient_barcode, 
                    histological_type,
                    age_at_initial_pathologic_diagnosis)
    
    stage <- .$clin$stage_event %>%
      dplyr::select(bcr_patient_barcode, clinical_stage)
    
    hpv <- .$hpv %>%
      dplyr::mutate(bcr_patient_barcode = substr(sample, 1, 12)) %>%
      dplyr::select(bcr_patient_barcode, hpv_status)
    
    endpoints <- .$endpoints %>%
      dplyr::mutate(pfs_time_month = pfs_time / 30.417)
    
    samples_wes <- .$SNVs %>%
      dplyr::select(tumor_sample_barcode)%>%
      dplyr::mutate(bcr_patient_barcode = 
                      substring(tumor_sample_barcode, 1, 12))%>%
      dplyr::group_by(bcr_patient_barcode)%>%
      dplyr::filter(row_number()==1)%>%
      dplyr::rename(wes_id = tumor_sample_barcode)
    
    samples_wgs <- .$CNVs %>%
      dplyr::select(gdc_aliquot_id)%>%
      dplyr::mutate(bcr_patient_barcode = substring(gdc_aliquot_id, 1, 12))%>%
      group_by(bcr_patient_barcode)%>%
      dplyr::filter(row_number()==1)%>%
      dplyr::rename(wgs_id = gdc_aliquot_id)
    
    samples_rna <- tcga_raw$RNA %>%
      dplyr::select(starts_with("tpm_unstranded_")) %>%
      colnames() %>%
      gsub("^tpm_unstranded_","", .)%>%
      tibble()%>%
      set_names(nm = "rna_id")%>%
      dplyr::mutate(bcr_patient_barcode = substring(rna_id, 1, 12))
    
    list(patient, stage, hpv, endpoints, 
         samples_wes, samples_wgs, samples_rna)%>%
      purrr::reduce(left_join, by = "bcr_patient_barcode")
  } 

tcga_clin_stand <- standardize_clin(
  df             = tcga_clin,
  cohort_name    = "TCGA-CESC",
  id_col         = "bcr_patient_barcode",
  wes_id_col     = "wes_id",
  wgs_id_col     = "wgs_id",
  rna_id_col     = "rna_id",
  time_col       = "pfs_time_month", 
  event_col      = "pfs", 
  age_col        = "age_at_initial_pathologic_diagnosis", 
  figo_col       = "clinical_stage", 
  hpv_col        = "hpv_status", 
  histo_col      = "histological_type")

# Prepare molecular data --------------------------------------------------

## Genomic data ----

### Gene-level and pathway-level data 

tcga_snv <- process_tcga_snv(tcga_raw$SNVs)
tcga_cnv <- process_tcga_cnv(tcga_raw$CNVs)$gene_status
tcga_dna_processed <- process_tcga_dna(tcga_snv, tcga_cnv)

tcga_dna <- tcga_dna_processed %>%
  pivot_wider_dna(mut_key             = "dna_pathway_long",
                  cohort_name         = "TCGA-CESC",
                  gene_col            = "gene_name",
                  id_col              = "bcr_patient_barcode",
                  altered_col         = "altered_gene",
                  cnv_col             = "cnv",
                  snv_col             = "snv", 
                  fusion_col          = NULL,
                  pathway_col         = "pathway",
                  altered_pathway_col = "altered_pathway")

## Transcriptomic data ----

### Gene-level data 

tcga_rna_prep_gene  <- prep_tcga_rna_gene(tcga_raw$RNA)

tcga_rna_gene <- tcga_rna_prep_gene %>%
  pivot_wider_rna(cohort_name = "TCGA-CESC", 
                  gene_col    = "gene_name",
                  id_col      = "bcr_patient_barcode",
                  rna_tpm_col = "tpm")

### Pathway-level data 

#' Gene Set Variation Analysis (GSVA)
#'The idea here is that we will get pathway-level scores for each sample that 
#'indicate if genes in a pathway vary concordantly in one direction 
#'(over-expressed or under-expressed relative to the overall population)
#'[@Hänzelmann 2013a]. This means that GSVA scores will depend on the 
#'samples included in the dataset when you run GSVA; if you added more samples 
#'and ran GSVA again, you would expect the scores to change [@Hänzelmann 2013a].

tcga_rna_prep_gsva  <- prep_tcga_rna(tcga_raw$RNA)

tcga_gsva <- list(
  clin_rna  = c("clin", "rna"),
  clin_dna_rna = c("clin", "dna", "rna")
) %>%
  purrr::imap(\(req, nm)
              tcga_rna_prep_gsva %>%
                dplyr::select(gene_name, gene_id,
                              dplyr::any_of(tcga_clin_stand %>% 
                                              drop_na(dplyr::all_of(req)) %>% 
                                              dplyr::pull(patient_id))
                ) %>%
                compute_gsva(cohort_name = "TCGA-CESC")
  )

# Data integration  -------------------------------------------------------

tcga <- integrate_data(
  id_col             = c("patient_id","cohort"),
  clin               = tcga_clin_stand,
  dna                = tcga_dna,
  rna_gene           = tcga_rna_gene,
  gsva_clin_rna      = tcga_gsva$clin_rna,
  gsva_clin_dna_rna  = tcga_gsva$clin_dna_rna, 
  keep_only_complete = TRUE, 
  pre_filter         = TRUE
)

saveRDS(tcga, here::here("data/processed/tcga/tcga_processed_data.RDS"))
