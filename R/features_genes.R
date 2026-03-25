
# ----------------------------- Documentation ------------------------------#

#'This file contains functions to:
#'- converts molecular data from long format to wide format
#'- filter gene-level features (DNA and RNA) and mutational signature 
#'variables. 
#'The goal is to retain only informative and reliable molecular features 
#'for association and survival analyses.

# ----------------------------- Dependencies -------------------------------#

library(dplyr)
library(tidyr)
library(janitor)
library(labelled)

# ------------------------------ Functions ---------------------------------#

#' Pivot DNA table wider (one row per patient)
#'
#' This function converts DNA-level molecular data from long format to wide
#' format, with one row per patient and one column per genomic feature.
#' Depending on the cohort, it builds patient-level indicators for altered
#' genes, copy-number variations (CNV), single-nucleotide variants (SNV),
#' fusion events, mutational signature subgroups, homologous recombination
#' deficiency (HRD), tumor mutational burden (TMB), microsatellite instability
#' (MSI), and altered DNA pathways.
#'
#' For the `"Bio-RAIDs"` cohort, additional features are derived from
#' `hrd_raw`, `tmb_msi_raw`, `signatures_subgroups_raw`, and
#' `patient_sample_seq_raw`. For the `"TCGA-CESC"` cohort, only gene-level and
#' pathway-level DNA features are generated from the table indicated by
#' `mut_key`.
#'
#' @param df A named list containing one or more raw DNA-related data frames.
#'   It must contain at least the element referenced by `mut_key`. For
#'   `"Bio-RAIDs"`, it must also contain `hrd_raw`, `tmb_msi_raw`,
#'   `signatures_subgroups_raw`, and `patient_sample_seq_raw`.
#'
#' @param mut_key A character string indicating the name of the element in `df`
#'   containing the mutation-level table to reshape.
#'
#' @param cohort_name A character string indicating the cohort name. Supported
#'   values in the current implementation are `"Bio-RAIDs"` and `"TCGA-CESC"`.
#'   The selected cohort determines which additional derived variables are
#'   created.
#'
#' @param gene_col A character string specifying the column containing gene
#'   identifiers in the mutation table. Default is `"gene"`.
#'
#' @param id_col A character string specifying the column containing patient
#'   identifiers. Default is `"patient_id"`. This column is renamed to
#'   `"patient_id"` in the final output.
#'
#' @param altered_col A character string specifying the column containing the
#'   binary indicator of whether a gene is altered. Default is
#'   `"altered_gene"`.
#'
#' @param cnv_col A character string specifying the column containing the
#'   binary indicator of copy-number variation events. Default is `"cnv"`.
#'
#' @param snv_col A character string specifying the column containing the
#'   binary indicator of single-nucleotide variant events. Default is `"snv"`.
#'
#' @param fusion_col A character string specifying the column containing the
#'   binary indicator of fusion events. Default is `"fusion"`. This argument is
#'   only used for the `"Bio-RAIDs"` cohort.
#'
#' @param pathway_col A character string specifying the column containing DNA
#'   pathway names. Default is `"pathway"`.
#'
#' @param altered_pathway_col A character string specifying the column
#'   containing the binary indicator of pathway-level alteration. Default is
#'   `"altered_pathway"`.
#'
#' @return A tibble in wide format with:
#' \itemize{
#'   \item one row per patient,
#'   \item a `cohort` column,
#'   \item a standardized `patient_id` column,
#'   \item gene-level DNA features prefixed with `"altered_"`, `"cnv_"`,
#'   and `"snv_"`,
#'   \item for `"Bio-RAIDs"`, additional features prefixed with `"fusion_"`,
#'   `"sig_group_"`, as well as `hrd`, `tmb_high`, and `msi_high`,
#'   \item pathway-level DNA features prefixed with `"genomic_pathway_"`,
#'   \item cleaned variable names using `janitor::clean_names()`,
#'   \item variable labels for metadata and molecular features.
#' }
#'
#'@export
pivot_wider_dna <- function(df, 
                            mut_key, 
                            cohort_name, 
                            gene_col = "gene",
                            id_col = "patient_id",
                            altered_col = "altered_gene",
                            cnv_col = "cnv",
                            snv_col = "snv", 
                            fusion_col = "fusion",
                            pathway_col = "pathway",
                            altered_pathway_col = "altered_pathway"){
  
  if (cohort_name=="Bio-RAIDs"){
    # HRD 
    hrd <- df[["hrd_raw"]]%>%
      janitor::clean_names()%>%
      dplyr::mutate(hrd = ifelse(shallow_hrd_status=="HRD", 1, 0))%>%
      labelled::set_variable_labels(hrd = "Homologous recombination deficiency")%>%
      dplyr::select(dplyr::all_of(id_col), hrd)
    
    # TMB / MSI 
    tmb_msi <- df[["tmb_msi_raw"]] %>%
      janitor::clean_names()%>%
      dplyr::mutate(tmb_high = ifelse(tmb_status == "TMB high", 1, 0),
                    msi_high = ifelse(msi_status == "MSI high", 1, 0))%>%
      labelled::set_variable_labels(
        tmb_high = "Tumor mutational burden high", 
        msi_high = "Microsatellite instability high")%>%
      dplyr::select(all_of(id_col), tmb_high, msi_high)
  }
  
  mut <- df[[mut_key]] %>%
    janitor::clean_names()
  
  # Altered genes 
  altered <- mut %>%
    dplyr::select(dplyr::all_of(c(id_col, gene_col, altered_col))) %>%
    dplyr::distinct()%>% 
    tidyr::pivot_wider(names_from = dplyr::all_of(gene_col),
                       values_from = dplyr::all_of(altered_col)) %>%
    labelled::set_variable_labels(.labels = paste("DNA alteration -", 
                                                  names(.)))%>%
    janitor::clean_names() %>%
    dplyr::rename_with(~ paste0("altered_", .x), -dplyr::all_of(id_col)) 
  
  # CNV
  cnv <- mut %>%
    dplyr::select(dplyr::all_of(c(id_col, gene_col, cnv_col))) %>%
    dplyr::distinct()%>% 
    tidyr::pivot_wider(names_from = dplyr::all_of(gene_col),
                       values_from = dplyr::all_of(cnv_col)) %>%
    labelled::set_variable_labels(.labels = paste("CNV -", names(.)))%>%
    janitor::clean_names() %>%
    dplyr::rename_with(~ paste0("cnv_", .x), -dplyr::all_of(id_col))
  
  # SNV
  snv <- mut %>%
    dplyr::select(dplyr::all_of(c(id_col, gene_col, snv_col))) %>%
    dplyr::distinct()%>% 
    tidyr::pivot_wider(names_from = dplyr::all_of(gene_col),
                       values_from = dplyr::all_of(snv_col)) %>%
    labelled::set_variable_labels(.labels = paste("SNV -", names(.)))%>%
    janitor::clean_names() %>%
    dplyr::rename_with(~ paste0("snv_", .x), -dplyr::all_of(id_col)) 
  
  if (cohort_name=="Bio-RAIDs"){
    # Fusion
    fusion <- mut %>%
      dplyr::select(dplyr::all_of(c(id_col, gene_col, fusion_col))) %>%
      dplyr::distinct()%>% 
      tidyr::pivot_wider(names_from = dplyr::all_of(gene_col),
                         values_from = dplyr::all_of(fusion_col)) %>%
      labelled::set_variable_labels(.labels = paste("Fusion -", names(.)))%>%
      janitor::clean_names() %>%
      dplyr::rename_with(~ paste0("fusion_", .x), -dplyr::all_of(id_col))
    
    # Mutation signature subgroups 
    sig <- df[["signatures_subgroups_raw"]] %>% 
      dplyr::rename(., signature = X)%>%
      tidyr::pivot_longer(-c(signature), names_to = "wes_tumor")%>%
      dplyr::inner_join(df[["patient_sample_seq_raw"]] %>% 
                          janitor::clean_names()%>%
                          dplyr::select(all_of(id_col), wes_tumor), 
                        by = "wes_tumor")%>%
      tidyr::pivot_wider(names_from = signature, values_from = value)%>%
      dplyr::select(-wes_tumor)%>%
      dpyr::mutate(across(-all_of(id_col), ~ifelse(.x!=0, 1, 0)))%>%
      labelled::set_variable_labels(
        .labels = paste("Mutational signature subgroup -", 
                        names(.)))%>%
      janitor::clean_names()%>%
      dplyr::rename_with(~ paste0("sig_group_", .x), -dplyr::all_of(id_col))
  }
  
  # DNA pathways 
  pathways_sort <- c(sort(unique(mut$pathway))%>% subset(.!="Others "),
                     "Others ")
  dna_pathway <- mut %>%
    dplyr::select(dplyr::all_of(c(id_col, pathway_col, altered_pathway_col))) %>%
    dplyr::distinct()%>% 
    tidyr::pivot_wider(names_from = dplyr::all_of(pathway_col),
                       values_from = dplyr::all_of(altered_pathway_col)) %>%
    labelled::set_variable_labels(.labels = trimws(names(.)))%>%
    dplyr::relocate(all_of(id_col), all_of(pathways_sort))%>%
    janitor::clean_names() %>%
    dplyr::rename_with(~ paste0("genomic_pathway_", .x), 
                       -dplyr::all_of(id_col)) 
  
  if (cohort_name=="Bio-RAIDs"){
    list_df <- list(hrd, tmb_msi, altered, cnv, snv, fusion, sig, dna_pathway)
  } else if (cohort_name=="TCGA-CESC"){
    list_df <- list(altered, cnv, snv, dna_pathway)
  }
  
  list_df %>%
    purrr::reduce(full_join, by = id_col)%>%
    dplyr::rename(patient_id = all_of(id_col))%>% 
    dplyr::mutate(cohort = cohort_name)%>%
    dplyr::relocate(cohort, patient_id)%>%
    labelled::set_variable_labels(patient_id = "Patient ID",
                                  cohort = "Cohort")
}

#' Pivot RNA table wider (one row per patient)
#'
#' @param df A data.frame or tibble containing RNA expression data in long format.
#'   Must include at least a patient identifier column, a gene column, 
#'   and an expression value column.
#' 
#' @param cohort_name A character string indicating the cohort name.
#'   This value will be added as a new column (`cohort`) for all observations.
#' 
#' @param gene_col A character string specifying the column name containing 
#'   gene identifiers. Default is `"gene_name"`. These values will become 
#'   column names in the output.
#' 
#' @param id_col A character string specifying the column name containing 
#'   patient identifiers. Default is `"patient_id"`. This column will be 
#'   renamed to `"patient_id"` in the output.
#' 
#' @param rna_tpm_col A character string specifying the column containing 
#'   RNA expression values (e.g., TPM). Default is `"tpm"`. These values 
#'   will populate the wide-format table.
#'
#' @return A tibble in wide format with:
#' \itemize{
#'   \item One row per patient
#'   \item One column per gene (prefixed with `"rna_seq_"`)
#'   \item A `cohort` column
#'   \item A standardized `patient_id` column
#'   \item Cleaned column names (via `janitor::clean_names()`)
#'   \item Variable labels for RNA features and metadata
#' }
#' @export
pivot_wider_rna <- function(df,
                            cohort_name, 
                            gene_col = "gene_name",
                            id_col = "patient_id",
                            rna_tpm_col = "tpm"){
  df %>%
    tidyr::pivot_wider(names_from = dplyr::all_of(gene_col), 
                       values_from = dplyr::all_of(rna_tpm_col))%>%
    labelled::set_variable_labels(.labels = paste("RNA-Seq -", names(.)))%>%
    janitor::clean_names()%>%
    dplyr::rename_with(~paste0("rna_seq_", .x), -dplyr::all_of(id_col))%>%
    dplyr::rename(patient_id = all_of(id_col))%>%
    dplyr::mutate(cohort = cohort_name)%>%
    dplyr::relocate(cohort, patient_id)%>%
    labelled::set_variable_labels(
      patient_id = "Patient ID",
      cohort = "Cohort")
}

#' Remove non-altered gene features
#'
#' This function filters out gene-level DNA features that show no alterations
#' across all samples. A gene is considered non-altered if the proportion of
#' altered samples (value equal to 1) is zero.
#'
#' The function operates on columns with specified prefixes corresponding to
#' different types of genomic alterations (e.g., `"altered_"`, `"cnv_"`,
#' `"snv_"`, `"fusion_"`).
#'
#' @param df A data.frame or tibble containing gene-level DNA features.
#'   Columns are expected to include binary indicators (0/1) of genomic
#'   alterations for each gene.
#'
#' @param prefix A character vector of prefixes used to identify gene-level
#'   DNA feature columns. Default is
#'   `c("altered_", "cnv_", "snv_", "fusion_")`.
#'
#' @return A data.frame with the same structure as the input, but with gene
#'   features removed if they have no alterations (i.e., all values are 0 or NA).
#'   
#' @export
filter_out_non_altered_genes <- function(df, 
                                         prefix = c("altered_", 
                                                    "cnv_", 
                                                    "snv_", 
                                                    "fusion_")){
  non_altered_genes <- df %>%
    dplyr::select(dplyr::starts_with(prefix)) %>%
    dplyr::summarise(
      dplyr::across(everything(),
                    ~ sum(. == 1, na.rm = TRUE) / dplyr::n())) %>%
    dplyr::select(where(~ .x == 0)) %>%
    names()
  
  data %>%
    dplyr::select(-dplyr::all_of(non_altered_genes))
}

#' Filter out non-expressed or lowly expressed genes
#'
#' This function removes gene-level RNA features with low expression across
#' samples. A gene is retained only if its expression is at least 1 TPM in
#' at least 20\% of samples.
#'
#' Columns are identified using the specified prefix, which is expected to
#' correspond to RNA expression variables.
#'
#' @param df A data.frame or tibble containing gene-level RNA features.
#'   RNA variables are expected to be numeric expression values, typically TPM.
#'
#' @param prefix A character string used to identify RNA feature columns.
#'   Default is `"rna_seq_"`.
#'
#' @return A data.frame with the same structure as the input, but with
#'   lowly expressed genes removed.
#'
#' @export
filter_out_low_expressed_genes <- function(df, prefix = "rna_seq_"){
  
  n_patients <- nrow(df)
  low_expressed_genes <- data %>%
    dplyr::select(starts_with(prefix))%>%
    # count number of patients with expressed genes
    dplyr::summarise(across(everything(),~sum(.x >= 1)))%>% 
    # select genes expressed in less than 20% of patients
    dplyr::select_if(.<0.2*n_patients)%>% 
    names()
  
  data %>% 
    dplyr::select(-all_of(low_expressed_genes))
}

#' Keep the main mutational signature groups
#'
#' This function retains the most frequent mutational signature groups
#' (default: top 6) and removes less frequent ones from the dataset.
#' It is designed for datasets containing binary indicators of mutational
#' signature subgroup membership (e.g., variables prefixed with `"sig_group_"`).
#'
#' @param df A data.frame or tibble containing mutational signature group
#'   features. These features are expected to be binary indicators (0/1).
#'
#' @return A data.frame with only the most frequent mutational signature
#'   groups retained (top 6 by proportion of altered samples).
#'
#' @export
filter_main_sig_group <- function(df){
  
  var_non_main_sig_group <- df %>%
    dplyr::select(starts_with("sig_group_")) %>%
    # mutational signature subgroups rates in the sample
    dplyr::summarise(across(everything(), ~sum(.==1) / n()))%>% 
    as.vector()%>%
    unlist()%>%
    sort(decreasing = TRUE)%>%
    .[-c(1:6)]%>%
    names()
  
  df %>%
    dplyr::select(-all_of(var_non_main_sig_group)) 
}
