
# ----------------------------- Documentation ------------------------------#

#'This file contains functions to standardize and harmonize clinical data.

# ----------------------------- Dependencies -------------------------------#

library(dplyr)
library(labelled)

# ------------------------------ Functions ---------------------------------#

#' Standardize clinical data across cohorts
#'
#' This function harmonizes clinical variables coming from heterogeneous cohort-
#' specific data frames into a common analysis-ready format. It extracts selected
#' columns, renames them to a standardized nomenclature, recodes several
#' clinical variables, creates binary and factor versions of key predictors, and
#' defines availability indicators for clinical, DNA, and RNA data.
#'
#' The function is designed for cohorts such as `"Bio-RAIDs"` and `"TCGA-CESC"`.
#'
#' @param df A data frame containing raw clinical and sample-level information.
#' @param cohort_name A character string indicating the cohort name. Used to
#' populate the standardized `cohort` column and to derive `cohort2`.
#' Supported values include `"Bio-RAIDs"` and `"TCGA-CESC"`.
#' @param id_col Character string. Name of the patient identifier column.
#' @param wes_id_col Character string. Name of the WES sample identifier column.
#' @param wgs_id_col Character string. Name of the WGS sample identifier column.
#' @param rna_id_col Character string. Name of the RNA-seq sample identifier column.
#' @param time_col Character string. Name of the follow-up time column.
#' @param event_col Character string. Name of the event indicator column.
#' @param age_col Character string. Name of the age column.
#' @param figo_col Character string. Name of the FIGO stage column.
#' @param hpv_col Character string. Name of the HPV status column.
#' @param histo_col Character string. Name of the histology column.
#' @param necrosis_col Optional character string. Name of the necrosis percentage
#' column. Default is `NULL`. If `NULL`, necrosis-derived variables are not kept
#' in the output.
#'
#' @details
#' The function performs the following transformations:
#'
#' \itemize{
#'   \item Standardizes selected columns into a common format.
#'   \item Converts identifiers to character and key quantitative variables to numeric.
#'   \item Creates a short cohort label:
#'   \itemize{
#'     \item `"Bio-RAIDs"` -> `"raids"`
#'     \item `"TCGA-CESC"` -> `"tcga"`
#'   }
#'   \item Dichotomizes age into:
#'   \itemize{
#'     \item `age_c_f`: factor with levels `"\u226450"` and `">50"`
#'     \item `age`: binary indicator (`0` for `\u226450`, `1` for `>50`)
#'   }
#'   \item Recodes FIGO stage into:
#'   \itemize{
#'     \item `figo_c_f`: factor with levels `"I/II"` and `"III/IV"`
#'     \item `figo`: binary indicator (`0` for stages I/II, `1` for stages III/IV)
#'   }
#'   \item Recodes HPV status into:
#'   \itemize{
#'     \item `hpv_negative`: binary indicator (`1` for negative, `0` otherwise)
#'     \item `hpv_negative_f`: factor with levels `"Positive"` and `"Negative"`
#'   }
#'   \item If `necrosis_col` is provided:
#'   \itemize{
#'     \item `necrosis`: binary indicator (`0` for no necrosis, `1` otherwise)
#'     \item `necrosis_f`: factor with levels `"No"` and `"Yes"`
#'   }
#'   \item Recodes histology into:
#'   \itemize{
#'     \item `"Squamous cell carcinoma"`
#'     \item `"Adenocarcinoma"`
#'     \item `"Other"`
#'   }
#'   \item Creates data availability indicators:
#'   \itemize{
#'     \item `clin`: clinical data availability
#'     \item `dna`: DNA data availability
#'     \item `rna`: RNA data availability
#'   }
#' }
#'
#' Clinical data availability is cohort-specific:
#' \itemize{
#'   \item For `"Bio-RAIDs"`, `clin = 1` if `age`, `figo`, `hpv_negative`, and
#'   `necrosis` are all non-missing.
#'   \item For other cohorts, `clin = 1` if `age`, `figo`, and `hpv_negative`
#'   are all non-missing.
#' }
#'
#' DNA availability is defined as non-empty WES and WGS identifiers.
#' RNA availability is defined as a non-empty RNA-seq identifier.
#'
#' @return A tibble with standardized columns. The output includes:
#' \describe{
#'   \item{cohort}{Original cohort name.}
#'   \item{cohort2}{Short cohort label.}
#'   \item{patient_id}{Standardized patient identifier.}
#'   \item{wes_id, wgs_id, rna_id}{Sequencing identifiers.}
#'   \item{time, event}{Survival outcome variables.}
#'   \item{clin, dna, rna}{Availability indicators.}
#'   \item{age_n, age_c_f, age}{Age variables.}
#'   \item{figo, figo_c_f}{Recoded FIGO variables.}
#'   \item{hpv, hpv_negative, hpv_negative_f}{HPV variables.}
#'   \item{necrosis, necrosis_f}{Necrosis variables, only if `necrosis_col` is provided.}
#'   \item{histo, histo_f}{Histology variables.}
#' }
#'
#' Variable labels are attached using `labelled::set_variable_labels()`.
#'
#' @export
standardize_clin <- function(df, cohort_name, 
                             id_col, wes_id_col, wgs_id_col, rna_id_col,
                             time_col, event_col, 
                             age_col, figo_col, hpv_col, histo_col, 
                             necrosis_col = NULL){
  
  out <- dplyr::as_tibble(df)%>%
    dplyr::transmute(
      cohort         = cohort_name,
      patient_id     = as.character(.data[[id_col]]),
      time           = as.numeric(.data[[time_col]]),
      event          = as.numeric(.data[[event_col]]),
      age_n          = as.numeric(.data[[age_col]]),
      figo_stage     = as.character(.data[[figo_col]]),
      hpv            = as.character(.data[[hpv_col]]),
      histo          = as.character(.data[[histo_col]]),
      necrosis_perc  = if (!is.null(necrosis_col)) 
        as.character(.data[[necrosis_col]]) else NA_character_, 
      wes_id         = as.character(.data[[wes_id_col]]),
      wgs_id         = as.character(.data[[wgs_id_col]]),
      rna_id         = as.character(.data[[rna_id_col]]), 
      dplyr::across(dplyr::everything(), ~.x)
    )%>%
    dplyr::mutate(
      cohort2 = case_when(cohort == "Bio-RAIDs"~"raids",
                          cohort == "TCGA-CESC"~"tcga"),
      age_c_f = as.factor(ifelse(age_n>50, ">50", "\u226450")),
      age = ifelse(age_n>50, 1, 0),
      figo_c_f = as.factor(case_when(
        gsub("Stage |A|B|1|2", "", figo_stage) %in% c("I","II")~"I/II", 
        gsub("Stage |A|B|1|2", "", figo_stage) %in% c("III","IV")~"III/IV")),
      figo = case_when(
        gsub("Stage |A|B|1|2", "", figo_stage) %in% c("I","II")~0,
        gsub("Stage |A|B|1|2", "", figo_stage) %in% c("III","IV")~1),
      hpv_negative = ifelse(is.na(hpv) | hpv =="Final type (hierarchical)", NA,
                            ifelse(!is.na(hpv) & hpv %in% c("NEG","negative"), 1, 0)),
      hpv_negative_f = factor(hpv_negative, 0:1, c("Positive","Negative")),
      necrosis = ifelse(is.na(necrosis_perc), NA,
                        ifelse(necrosis_perc == 0, 0, 1)),
      necrosis_f = factor(necrosis, 0:1, c("No","Yes")),
      histo_f = factor(
        ifelse(histo %in% c("Squamous cell carcinoma (epidermoid)",
                            "Cervical Squamous Cell Carcinoma"), 
               "Squamous cell carcinoma", 
               ifelse(histo == "Adeno-carcinoma" | grepl("Adenocarcinoma", histo),
                      "Adenocarcinoma", "Other")),
        levels = c("Squamous cell carcinoma", "Adenocarcinoma", "Other")), 
      
      clin = dplyr::case_when(
        cohort == "Bio-RAIDs" &
          !is.na(age) &
          !is.na(figo) &
          !is.na(hpv_negative) &
          !is.na(necrosis) ~ 1,
        
        cohort != "Bio-RAIDs" &
          !is.na(age) &
          !is.na(figo) &
          !is.na(hpv_negative) ~ 1,
        
        TRUE ~ NA_real_
      ),
      dna = ifelse(wes_id!="" & wgs_id!="", 1, NA),
      rna = ifelse(rna_id!="", 1, NA),
    )%>%
    labelled::set_variable_labels(
      patient_id = "Patient ID",
      age_n = "Age",
      age_c_f = "Age",
      age = "Age > 50 years",
      figo_c_f = "FIGO stage",
      figo = "FIGO Stages III and IV",
      hpv = "HPV status",
      hpv_negative = "HPV negative",
      hpv_negative_f = "HPV status",
      necrosis = "Necrosis",
      necrosis_f = "Necrosis",
      histo = "Histology",
      histo_f = "Histology"
    )%>%
    dplyr::select(cohort, cohort2, patient_id, wes_id, wgs_id, rna_id, 
                  time, event, clin, dna, rna, 
                  age_n, age_c_f, age, figo, figo_c_f, 
                  hpv, hpv_negative, hpv_negative_f, 
                  necrosis, necrosis_f, histo, histo_f)
  if (is.null(necrosis_col)) {
    out %>% dplyr::select(-c(necrosis, necrosis_f))
  } else { out }
}
