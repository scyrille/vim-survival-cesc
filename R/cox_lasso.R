
# ----------------------------- Documentation ------------------------------#

#' This file contains functions to fit Cox regression models with lasso 
#' penalties and likelihood-based boosting. 

# ----------------------------- Dependencies -------------------------------#

library(tidyverse)
library(labelled)
library(gt)
library(survival)
library(glmnet)   
library(grpreg)
library(SGL)  
library(prioritylasso)
library(gtools)
library(CoxBoost)
library(mboost)

# ------------------------------ Functions ---------------------------------#

# Statistical helpers -----------------------------------------------------

# 1. Lasso (glmnet)

#' Fit a Cox regression model via penalized maximum likelihood 
#' @param time numeric vector of observed follow-up times. 
#' @param event numeric vector of status indicators of whether an event was 
#' observed or not.
#' @param X data.frame of observed covariate values.
#' @param nfolds number of folds.
#' @param seed Optional seed for reproducibility.
#' 
fit_cox_lasso <- function(time, event, X, nfolds = 5, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  Y <- cbind(time, status = event)
  if (class(X)[1]!="matrix"){X <- as.matrix(X)}
  cv.glmnet <- glmnet::cv.glmnet(X, 
                                 Y, 
                                 family = "cox",  
                                 alpha = 1, 
                                 nfolds = nfolds)
  glmnet <- glmnet::glmnet(X, 
                           Y, 
                           family = "cox", 
                           alpha = 1)
  
  fit <- list(cv.glmnet = cv.glmnet,
              glmnet = glmnet)
  
  #' Get coefficients 
  .get_coef_cox_lasso <- function(fit){
    
    coef <- fit[["glmnet"]]$beta[,which.min(fit[["cv.glmnet"]]$cvm)]
    variable <- names(coef)
    tbl_coef <- data.frame(variable, coef, row.names = NULL)
    return(tbl_coef)
  }
  
  coef <- .get_coef_cox_lasso(fit) %>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  # Selected features 
  coef_select <- coef %>%
    dplyr::filter(coef!=0)
  var_select <- coef_select$variable
  
  list(fit = fit,
       coef = coef,
       coef_select = coef_select,
       var_select = var_select)
}

# 2. Group Lasso (grpreg)

#'Fit a penalized Cox model with grouped covariates.

#' @param time numeric vector of observed follow-up times. 
#' @param event numeric vector of status indicators of whether an event was 
#' observed or not.
#' @param X data.frame of observed covariate values.
#' @param group vector describing the grouping of the covariates in the 
#' Group Lasso Cox model.  
#' @param nfolds number of folds.
#' @param seed Optional seed for reproducibility.
#' 
fit_cox_group_lasso <- function(time, event, X, group, nfolds = 5, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  Y <- cbind(time, status = event)
  if (class(X)[1]!="matrix"){X <- as.matrix(X)}

  cv.grpsurv <- grpreg::cv.grpsurv(X = X, 
                                   y = Y, 
                                   group = group, 
                                   penalty = "grLasso",
                                   nfolds = nfolds, 
                                   returnY = T)
  grpsurv <- grpreg::grpsurv(X = X, 
                             y = Y, 
                             group = group,
                             penalty = "grLasso")
  
  fit <- list(cv.grpsurv = cv.grpsurv,
              grpsurv = grpsurv)
  
  # Get coefficients 
  .get_coef_cox_group_lasso <- function(fit){
    
    coef <- fit[["cv.grpsurv"]]$fit$beta[,which.min(fit[["cv.grpsurv"]]$cve)]
    variable <- names(coef)
    tbl_coef <- data.frame(variable, coef, row.names = NULL)
    return(tbl_coef)
  }
  
  coef <- .get_coef_cox_group_lasso(fit) %>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  # Selected features 
  coef_select <- coef %>%
    dplyr::filter(coef!=0)
  var_select <- coef_select$variable
  
  list(fit = fit,
       coef = coef,
       coef_select = coef_select,
       var_select = var_select)
}

# 3. Sparse group lasso (SGL)

#' Fits a regularized Cox model via penalized maximum likelihood that 
#' encourages the sparsity both on a group and within the group level.

#' @param time numeric vector of observed follow-up times. 
#' @param event numeric vector of status indicators of whether an event was 
#' observed or not.
#' @param X data.frame of observed covariate values.
#' @param index vector indicating group membership of each covariate in the 
#' Sparse Group Lasso Cox model. 
#' @param nfolds number of folds.
#' @param seed Optional seed for reproducibility.
#' 
fit_cox_sparse_group_lasso <- function(time, 
                                       event,
                                       X, 
                                       index, 
                                       nfolds = 5,
                                       seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  if (class(X)[1]!="matrix"){X <- as.matrix(X)}
  data <- list(x = X, time = time, status = event)
  cvSGL <- SGL::cvSGL(data = data, 
                      index = index, 
                      type = "cox", 
                      # nlam = 100,
                      nfold = nfolds)
  SGL <- SGL::SGL(data = data, 
                  index = index, 
                  type = "cox"#, 
                  # nlam = 100
                  )
  fit <- list(cvSGL = cvSGL, 
              SGL = SGL)
  
  # Get coefficients 
  .get_coef_cox_sparse_group_lasso <- function(fit){
    
    coef <- fit[["cvSGL"]]$fit$beta[, which.min(fit[["cvSGL"]]$lldiff)]
    variable <- names(fit[["SGL"]][["X.transform"]][["X.means"]])
    tbl_coef <- data.frame(variable, coef, row.names = NULL)
    return(tbl_coef)
  }
  
  coef <- .get_coef_cox_sparse_group_lasso(fit) %>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  # Selected features 
  coef_select <- coef %>%
    dplyr::filter(coef!=0)
  var_select <- coef_select$variable
  
  list(fit = fit,
       coef = coef,
       coef_select = coef_select,
       var_select = var_select)
}

# 4. Priority lasso (prioritylasso)

#' Fit successive Lasso models for several ordered blocks of variables and 
#' takes the predicted values as an offset for the next block.

#' @param time numeric vector of observed follow-up times. 
#' @param event numeric vector of status indicators of whether an event was 
#' observed or not.
#' @param X data.frame of observed covariate values.
#' @param blocks_list list of blocks containing the indices of variables in 
#' the priority Lasso Cox model. 
#' @param nfolds number of folds.
#' @param seed Optional seed for reproducibility.
#' 
fit_cox_priority_lasso <- function(time, 
                                   event, 
                                   X,  
                                   blocks_list,  
                                   nfolds = 5,
                                   block1_penalization = TRUE,
                                   seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  if (class(X)[1]!="matrix"){X <- as.matrix(X)}
  cvm_prioritylasso <- 
    prioritylasso::cvm_prioritylasso(X = X,  
                                     Y = survival::Surv(time, event),
                                     family = "cox", 
                                     type.measure = "deviance", 
                                     blocks.list = blocks_list, 
                                     nfolds = nfolds,
                                     block1.penalization = block1_penalization)
  fit <- cvm_prioritylasso
  
  .get_coef_priority_lasso <- function(fit){
    coef <- fit[["coefficients"]]
    variable <- names(coef)
    tbl_coef <- data.frame(variable, coef, row.names = NULL)
  }
  
  coef <- .get_coef_priority_lasso(fit) %>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  # Selected features 
  coef_select <- coef %>%
    dplyr::filter(coef!=0)
  var_select <- coef_select$variable
  
  list(fit = fit,
       coef = coef,
       coef_select = coef_select,
       var_select = var_select)
}

#' Generate ordered block lists for priority lasso
#'
#' This internal helper function generates all possible ordered combinations
#' (permutations) of input variable index groups and formats them as named
#' block lists (e.g., `bp1`, `bp2`, `bp3`). It is primarily designed for use
#' in block-based modeling approaches such as priority lasso, where different 
#' block orderings may affect model fitting.
#'
#' @param ... One or more vectors of variable indices (e.g., clinical, DNA,
#'   RNA indices). Each input corresponds to a block of variables.
#'
#' @return A list of lists. Each element corresponds to one permutation of
#'   the input blocks, where:
#' \itemize{
#'   \item each sublist represents a specific ordering of blocks,
#'   \item block names are assigned as `"bp1"`, `"bp2"`, ..., according to
#'   their position in the ordering.
#' }
.make_blocks_priority_lasso <- function(...) {
  x <- list(...)
  # x <- enquos(...)
  perms <- gtools::permutations(length(x), length(x))
  
  map(seq_len(nrow(perms)), function(i) {
    out <- x[perms[i, ]]
    names(out) <- paste0("bp", seq_along(out))
    out
  })
}

# Plotting and table helpers ----------------------------------------------

#' Plot selected coefficients
#'
#' This function plots the non-zero coefficients selected from a fitted
#' CoxBoost model at the optimal boosting step. Coefficients are displayed as
#' horizontal bars and can optionally be colored by data type
#' (clinical, genomic, or transcriptomic).
#'
#' @param fit_obj A fitted model object returned by `fit_coxboost()`. It must
#'   contain a component named `coef_select`, which is expected to be a
#'   data.frame or tibble with at least the following columns:
#'   \itemize{
#'     \item `variable`: variable names,
#'     \item `label`: variable labels to display on the plot,
#'     \item `coef`: estimated coefficients.
#'   }
#'
#' @param data_type Logical; if `TRUE`, bars are colored according to variable
#'   type (`"Clinical"`, `"Genomic"`, or `"Transcriptomic"`). If `FALSE`,
#'   bars are colored according to the sign of the coefficient (positive or
#'   negative). Default is `FALSE`.
#'
#' @return A `ggplot2` object.
#' 
#' @export
plot_coef <- function(fit_obj, data_type = F){
  
  coef <- fit_obj$coef_select %>% 
    arrange(coef)%>%
    mutate(label = as_factor(label))%>%
    arrange(coef)
  
  if (data_type) {
    
    coef <- coef %>%
      dplyr::mutate(
        data = case_when(
          variable %in% c("age","figo","necrosis","hpv_negative")~
            "Clinical",
          grepl(paste(c("hrd","tmb","msi","cnv_","snv_","fusion_",
                        "altered_","sig_group","genomic_pathway"), 
                      collapse = "|"), variable) ~"Genomic",
          grepl("rna_seq_|hallmark_", variable) ~"Transcriptomic"
          )
        )
    
    color <- coef %>%
      dplyr::mutate(color = case_when(data == "Clinical"~"#0AB329",
                                      data == "Genomic"~"#FF5733",
                                      data == "Transcriptomic"~"#0771A2"))%>%
      dplyr::select(color)%>%
      pull()
    
    plot <- ggplot(coef, aes(x=label, y=coef, label=coef)) +
      geom_bar(stat='identity', aes(fill=factor(data)), width=.5) +
      scale_fill_manual("", 
                        breaks = c("Clinical","Genomic","Transcriptomic"),
                        values = c("#0AB329","#FF5733","#0771A2"))+
      labs(x = "Variable", y = "Coefficient")+
      coord_flip()+
      theme_classic()+
      theme(legend.position = "none",
            axis.text.y = element_text(color = color),  
            text = element_text(size = 13),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 13))
  } else {
    
    coef <- coef %>%
      dplyr::mutate(sign = if_else(coef>0, "p","n"))
                    
    plot <- ggplot(coef, aes(x=label, y=coef, label=coef)) +
      geom_bar(stat='identity', aes(fill=sign), width=.5) +
      scale_fill_manual(values = c("p"="#f8766d", "n"="#00ba38")) +
      labs(x = "Variable", y = "Coefficient")+
      coord_flip()+
      theme_classic()+
      theme(legend.position = "none",
            text = element_text(size = 13),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 13))
  }
  
  plot
}

#' Summary table with number of selected features

#'@param fit_obj

tbl_n_select_coef <- function(fit_obj){

  coef_select <- fit_obj$coef_select
  n_clin = nrow(coef_select %>%
                  dplyr::filter(variable %in% c("age",
                                                "figo",
                                                "hpv_negative",
                                                "necrosis")))
  n_DNA = nrow(coef_select %>%
                 dplyr::filter(str_starts(variable,
                                          paste(c("altered",
                                                  "cnv",
                                                  "snv",
                                                  "fusion",
                                                  "tmb",
                                                  "msi",
                                                  "hrd",
                                                  "sig_group",
                                                  "genomic_pathway"),
                                                collapse = "|"))))

  n_RNA = nrow(coef_select %>%
                 dplyr::filter(str_starts(variable,
                                          paste(c("rna_seq",
                                                  "hallmark"),
                                                collapse = "|"))))

  n_tot = n_clin + n_DNA + n_RNA
  tibble(n_clin, n_DNA, n_RNA, n_tot)
}

#' Build a summary table of selected coefficients across data types
#'
#' This function takes a named list of summary tables and combines them into a
#' single formatted `gt` table. Each element of `tbl_list` is expected to be a
#' data frame containing at least the columns:
#' `n_clin`, `n_DNA`, `n_RNA`, and `n_tot`.
#'
#' The names of `tbl_list` are used to define the data type represented by each
#' table, and are converted into a human-readable factor with a predefined order.
#' DNA and RNA counts are replaced by an em dash when the corresponding data type
#' does not include DNA or RNA features.
#'
#' Expected list names are:
#' - `"gene_clin_dna"`
#' - `"gene_clin_rna"`
#' - `"gene_clin_dna_rna"`
#' - `"pathway_clin_dna"`
#' - `"pathway_clin_rna"`
#' - `"pathway_clin_dna_rna"`
#'
#' @param tbl_list A named list of data frames or tibbles. Each element should
#'   contain summary counts for selected coefficients.
#'
#' @return A `gt_tbl` object.
.custom_list_tbl_n_select_coef <- function(tbl_list) {
  
  data_levels <- c(
    "gene_clin_dna",
    "gene_clin_rna",
    "gene_clin_dna_rna",
    "pathway_clin_dna",
    "pathway_clin_rna",
    "pathway_clin_dna_rna"
  )
  
  data_labels <- c(
    "Clinical + DNA (gene-level)",
    "Clinical + RNA (gene-level)",
    "Clinical + DNA (gene-level) + RNA (gene-level)",
    "Clinical + DNA (pathway-level)",
    "Clinical + RNA (pathway-level)",
    "Clinical + DNA (pathway-level) + RNA (pathway-level)"
  )
  
  tbl_list %>%
    purrr::imap(~ .x %>%
                  dplyr::mutate(
                    data = factor(.y, levels = data_levels, labels = data_labels),
                    .before = n_clin
                  ) %>%
                  dplyr::mutate(
                    n_DNA = dplyr::if_else(!grepl("DNA", data), "\u2014", as.character(n_DNA)),
                    n_RNA = dplyr::if_else(!grepl("RNA", data), "\u2014", as.character(n_RNA))
                  )
    ) %>%
    dplyr::bind_rows() %>%
    gt::gt() %>%
    gt::cols_label(
      data   = gt::md("**Data type**"),
      n_clin = gt::md("**n. Clinical**"),
      n_DNA  = gt::md("**n. DNA**"),
      n_RNA  = gt::md("**n. RNA**"),
      n_tot  = gt::md("**Total**")
    ) %>%
    gt::cols_align(align = "left")
}
