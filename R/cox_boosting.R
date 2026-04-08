
# ----------------------------- Documentation ------------------------------#

#' This file contains functions to fit Cox regression models with model-based 
#' and likelihood-based boosting. 

# ----------------------------- Dependencies -------------------------------#

library(tidyverse)
library(labelled)
library(gt)
library(survival)
library(CoxBoost)
library(mboost)

# ------------------------------ Functions ---------------------------------#

# Statistical helpers -----------------------------------------------------

# Likelihood-based boosting 

#' Fit a Cox model by likelihood-based boosting.

#' @param time Numeric vector of observed follow-up times.
#' @param event Integer/numeric vector (0/1) indicating event occurrence.
#' @param X A data.frame (or tibble) of covariates (one row per subject).
#' @param mandatory Character vector of covariate names to be fitted unpenalized.
#'   Must be subset of colnames(X). Use character(0) if none.
#' @param nfolds Number of CV folds (K).
#' @param maxstepno Maximum number of boosting steps evaluated in CV.
#' @param cv_type CV type passed to cvcb.control (default: "verweij").
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with fitted model and tidy coefficient tables.
fit_cox_lboost <- function(time, 
                           event, 
                           X,
                           mandatory = NULL, 
                           nfolds = 10,
                           stepno = 200,
                           cv_type = "verweij",
                           seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  data <- data.frame(time = time, event = event, X)
  
  # Fit CoxBoost model 
  fit <- CoxBoost::iCoxBoost(
    formula = survival::Surv(time, event) ~., 
    data = data,
    cv = CoxBoost::cvcb.control(K=nfolds, type = cv_type),
    stepno = stepno,
    mandatory = mandatory)
  
  # Optimal step 
  opt_step <- fit$cv.res$optimal.step
  
  #' Get coefficients at optimal step
  .get_coef_cox_lboost <- function(fit){
    
    coef <- coef(fit, at.step = fit[["cv.res"]][["optimal.step"]])
    variable <- names(coef)
    tbl_coef <- data.frame(variable, coef, row.names = NULL)
    return(tbl_coef)
  }
  
  coef <- .get_coef_cox_lboost(fit) %>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  # Selected features 
  coef_select <- coef %>%
    dplyr::filter(coef!=0)
  var_select <- coef_select$variable
  
  # Get coefficients across boosting steps
  .get_coef_cox_lboost_step <- function(fit, 
                                        steps, 
                                        vars = NULL) {
    
    map_dfr(steps, function(s) {
      b <- coef(fit, at.step = s)
      
      if (!is.null(vars)) {
        b <- b[names(b) %in% vars]
      }
      
      tibble(
        variable = names(b),
        coef     = as.numeric(b),
        selected = coef,
        step     = s
      )
    })
  }
  # Coef path
  coef_step <- .get_coef_cox_lboost_step(fit, steps = 0:opt_step)%>%
    dplyr::left_join(lab_map, by = "variable") %>%
    dplyr::mutate(label = dplyr::coalesce(label, variable))
  
  list(fit = fit, 
       opt_step = fit$cv.res$optimal.step, 
       coef = coef,
       coef_select = coef_select,
       var_select = var_select,
       coef_step = coef_step)
}

# Model-based boosting 

#' Fit a Cox model by model-based boosting.

#' @param time numeric vector of observed follow-up times. 
#' @param event numeric vector of status indicators of whether an event was 
#' observed or not.
#' @param X data.frame of observed covariate values.
#' @param nfolds number of folds.
#' @param seed Optimal seed for reproducibility. 

fit_cox_mboost <- function(time, event, X, nfolds = 10, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  lab_map <- .get_var_labels(X)
  
  data <- data.frame(time = time, event = event, X)
  
  # time <- time 
  # event <- event 
  # data <- data.frame(time, event, X)
  # 
  glmboost <- mboost::glmboost(survival::Surv(time, event) ~., 
                               data = data, 
                               family = CoxPH())
  cvrisk <- mboost::cvrisk(glmboost, 
                           folds = mboost::cv(model.weights(glmboost), 
                                              type = "kfold", B = nfolds))
  fit <- list(glmboost = glmboost,
              cvrisk = cvrisk)
  
  .get_coef_cox_mboost <- function(fit){
    coef <- coef(fit[["glmboost"]][mboost::mstop(fit[["cvrisk"]])])
    variable <- names(coef)
    tbl_coef <- data.frame(variable, coef, row.names = NULL)%>%
      dplyr::filter(variable!="(Intercept)")
  }
  
  if (mstop(fit[["cvrisk"]])!=0) {
    coef <- .get_coef_cox_mboost(fit) %>%
      dplyr::left_join(lab_map, by = "variable") %>%
      dplyr::mutate(label = dplyr::coalesce(label, variable))
    
    # Selected features 
    coef_select <- coef %>%
      dplyr::filter(coef!=0)
    var_select <- coef_select$variable
    
  } else {
    
    coef <- NULL
    coef_select <- NULL
    var_select <- NULL 
  }
  
  list(fit = fit,
       coef = coef,
       coef_select = coef_select,
       var_select = var_select)
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

#' Plot coefficient updates at a specific boosting step
#'
#' @param fit_obj Output of fit_coxboost()
#' @param step_show Step number to display (integer)
#' 
#' @return A ggplot object
plot_coef_updates <- function(fit_obj,
                              step_show) {
  
  dat <- fit_obj$coef_step %>%
    dplyr::filter(variable %in% fit_obj$var_select)%>%
    mutate(type = case_when(
      variable %in% c("age","figo","necrosis","hpv_negative")~
        "Clinical",
      grepl(paste(c("hrd","tmb","msi","cnv_","snv_","fusion_",
                    "altered_","sig_group","genomic_pathway"),
                  collapse = "|"), variable) ~"Genomic",
      grepl("rna_seq_|hallmark_", variable) ~"Transcriptomic"
    )) %>%
    arrange(label, step) %>%
    group_by(label) %>%
    mutate(coef_prev = lag(coef, default = 0),
           delta = coef - coef_prev) %>%
    ungroup() %>%
    filter(step == step_show) %>%
    arrange(coef)%>%
    mutate(label=as_factor(label))%>%
    arrange(coef)
  
  color <- dat %>%
    mutate(color = case_when(type == "Clinical"~"#0AB329",
                             type == "Genomic"~"#FF5733",
                             type == "Transcriptomic"~"#0771A2"))%>%
    dplyr::select(color)%>%
    pull()
  
  ggplot(dat, aes(x = label)) +
    # barre = beta courant
    geom_col(aes(y = coef, fill = factor(type)), width = 0.5) +
    # segment = update (entre prev et curr)
    geom_segment(aes(xend = label,
                     y = 0,
                     yend = coef_prev),
                 linewidth = 11, colour = "black", alpha = 0.4) +
    scale_fill_manual("",
                      breaks = c("Clinical","Genomic","Transcriptomic"),
                      values = c("#0AB329","#FF5733","#0771A2")) +
    scale_y_continuous(limits = c(-1.5,1.5),
                       breaks = seq(-1.5, 1.5, .5))+
    coord_flip() +
    labs(#title = paste("Boosting step n°", step_show),
      x = "Variable", y = "Coefficient") +
    theme_classic() +
    theme(axis.text.y = element_text(color = color),
          legend.position = "none",
          text = element_text(size = 13),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 13))
}

#' Plot coefficient paths 
#'
#' @param fit_obj Output of fit_cox_boost()
#' @param step_show Step number to display (integer)
#' 
#' @return A ggplot object.
#' 
#' @export 
plot_coef_paths <- function(fit) {
  
  dat <- fit$coef_step %>%
    dplyr::filter(variable %in% fit$var_select)%>%
    mutate(data = case_when(
      variable %in% c("age","figo","necrosis","hpv_negative")~
        "Clinical",
      grepl(paste(c("hrd","tmb","msi","cnv_","snv_","fusion_",
                    "altered_","sig_group","genomic_pathway"),
                  collapse = "|"), variable)~
        "Genomic",
      grepl("rna_seq_|hallmark_", variable) ~
        "Transcriptomic"
    ))
  
  ggplot(dat, aes(x = step, y = coef, group = label, color = data)) +
    # geom_hline(yintercept = 0) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    # geom_text(
    #   data = dat,
    #   aes(label = label),
    #   vjust = 1.5,
    #   size = 3,
    #   show.legend = FALSE
    # ) +
    # scale_x_continuous(minor_breaks = 1, limits = c(0, max(dat$step)+0.5))+
    scale_x_continuous(breaks = seq(0, max(dat$step), 1), 
                       limits = c(0, max(dat$step)))+
    scale_color_manual(values = c(
      Clinical = "#0AB329",
      Genomic = "#FF5733",
      Transcriptomic = "#0771A2",
      `Non selected`="gray"
    )) +
    labs(x = "Boosting step", y = expression(beta)) +
    theme_classic()+
    theme(legend.position = "none")
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

