
# ----------------------------- Documentation ------------------------------#

#' This file contains functions to:
#' 1. Make survival tables and survival curves
#' 2. Fit a Cox proportional hazards model
#' 3. Compute global tests

# ----------------------------- Dependencies -------------------------------#

library(tidyverse)   
library(janitor)
library(survival)      
library(survminer)    
library(gtsummary)
library(globaltest)

# ------------------------------ Functions --------------------------------#

#' Plot Kaplan–Meier survival curves (null or adjusted model)
#'
#' This function produce survival curves, either:
#'  - null model:   Surv(time, status) ~ 1       (no groups, no legend)
#'  - adjusted model: Surv(time, status) ~ group (groups as strata, legend shown)
#'
#' @param formula A survival formula, e.g. `Surv(time, status) ~ 1` or
#'   `Surv(time, status) ~ arm`.
#' @param data A data.frame containing the variables in `formula`.
#' @param title Plot title.
#' @param legend Position of the legend for adjusted models
#'   (e.g. `"top"`, `"bottom"`, `"right"`, `"left"`, `"none"`).
#' @param legend.title Legend title. 
#' @param xlab Label for x-axis (time).
#' @param ylab Label for y-axis (survival probability).
#' @param xlim Numeric length-2 vector for x-axis limits.
#' @param break.x.by Numeric: distance between x-axis breaks.
#' @param conf.int Logical: show confidence intervals.
#' @param risk.table Character or logical: type of risk table
#'   (e.g. `"nrisk_cumcensor"`, `TRUE`, or `"none"`).
#' @param risk.table.y.text Logical: show group labels in risk table y-axis.
#' @param risk.table.fontsize Numeric: font size for risk table.
#' @param tables.height Numeric: relative height of risk table vs plot.
#' @param surv.scale Character: `"percent"` or `"default"`.
#' @param axes.offset Logical: add offset to axes (ggsurvplot argument).
#' @param gg_theme A ggplot2 theme object for the main plot.
#' @param tables.theme A ggplot2 theme object for the risk table.
#' @param pval Logical: add log-rank p-value (for adjusted models).
#' @param ... Additional arguments passed to `ggsurvplot()`.
#'
#' @return A `ggsurvplot` object (list with ggplot & table components).

plot_surv <- function(formula,
                      data, 
                      title = "",
                      legend = "top", 
                      legend.title = "",
                      xlab = "Time in months",
                      ylab = "Progression-free survival",
                      xlim = c(0,37),  
                      break.x.by = 6,
                      conf.int = F, 
                      risk.table = "nrisk_cumcensor",
                      risk.table.y.text = F,
                      risk.table.fontsize = 4,
                      tables.height = 0.13,
                      surv.scale = "percent",
                      axes.offset = T,
                      gg_theme = theme_bw()+
                        theme(axis.title = element_text(size = 9),
                              axis.text  = element_text(size = 9)),
                      tables.theme = theme_cleantable()+
                        theme(plot.title = element_text(size = 10)),
                      pval = F,
                      ...){
  
  fit <- surv_fit(formula, data, match.fd = FALSE)
  plot <- vector(mode = "list", length(fit))
  
  for (i in seq_along(fit)){
    if (grepl("null_model", names(fit)[i])){ 
      
      plot[[i]] <- ggsurvplot(
        fit = fit[[i]], 
        legend = "none", 
        title = title, 
        xlab = xlab,
        ylab = ylab,
        xlim = xlim,  
        break.x.by = break.x.by,
        conf.int = conf.int, 
        risk.table = risk.table,
        risk.table.y.text = risk.table.y.text,
        risk.table.fontsize = risk.table.fontsize, 
        tables.height = tables.height,
        gg_theme = gg_theme, 
        tables.theme = tables.theme,
        surv.scale = surv.scale,
        axes.offset = axes.offset,
        ...
      ) } else {
        
        names(fit[[i]]$strata) <- sub(".*?=", "", names(fit[[i]]$strata)) 
        
        plot[[i]] <- 
          ggsurvplot(
            fit = fit[[i]], 
            xlab = xlab,
            ylab = ylab,
            title = title,
            legend = legend,
            legend.title = legend.title[i],
            xlim = xlim,  
            break.x.by = break.x.by,
            conf.int = conf.int, 
            risk.table = risk.table,
            risk.table.y.text = risk.table.y.text,
            risk.table.fontsize = risk.table.fontsize, 
            gg_theme = gg_theme, 
            tables.theme = tables.theme,
            surv.scale = surv.scale,
            axes.offset = axes.offset,
            pval = pval,
            ...
          )
      }
  }
  names(plot) <- gsub("::", "_", names(fit))
  plot
}

#' Create stratified Cox proportional hazards regression tables
#'
#' This function fits Cox proportional hazards regression models within levels
#' of a stratification variable and summarizes the results using `gtsummary`.
#' It supports both univariable and multivariable Cox models and returns a
#' stratified regression table.
#'
#' For each stratum, the function computes a Cox model and formats the output
#' with hazard ratios, 95\% confidence intervals, p-values, number of events,
#' and sample sizes.
#'
#' @param data A data.frame or tibble containing survival variables, the
#'   stratification variable, and covariates.
#'
#' @param time A character string specifying the name of the survival time
#'   variable.
#'
#' @param event A character string specifying the name of the event indicator
#'   variable.
#'
#' @param strata A character string specifying the name of the variable used
#'   for stratification. Separate Cox regression tables are computed within
#'   each level of this variable.
#'
#' @param covars A character vector containing the names of the covariates to
#'   include in the Cox model.
#'
#' @param model A character string indicating the model type. Must be either
#'   `"univ"` for univariable Cox regression or `"multi"` for multivariable
#'   Cox regression. Default is `c("univ", "multi")`.
#'
#' @param exponentiate Logical; if `TRUE`, regression coefficients are
#'   exponentiated and reported as hazard ratios. Default is `TRUE`.
#'
#' @param tidy_fun A function used to tidy model results before table
#'   formatting. Default is `broom.helpers::tidy_parameters`.
#'
#' @return A `gtsummary` stratified table object combining Cox regression
#'   results across strata.
#'
#' @export
tbl_cox_strata <- function(data,
                           time,
                           event,
                           strata,
                           covars,
                           model = c("univ", "multi"),
                           exponentiate = TRUE,
                           tidy_fun = broom.helpers::tidy_parameters) {
  model <- match.arg(model)
  
  df <- data %>%
    dplyr::select(all_of(c(time, event, strata, covars)))
  
  .tbl_fun <- function(.x) {
    if (model == "univ") {
      .x %>%
        gtsummary::tbl_uvregression(
          method = survival::coxph,
          y = survival::Surv(rlang::eval_tidy(time), rlang::eval_tidy(event)),
          include = covars,
          exponentiate = exponentiate,
          hide_n = TRUE,
          tidy_fun = tidy_fun
        )
    } else {
      
      fml <- stats::as.formula(
        paste0(
          "survival::Surv(",
          time, ", ",
          event,
          ") ~ ",
          paste(vapply(covars, as.character, character(1)), collapse = " + ")
        )
      )
      
      survival::coxph(fml, data = .x) %>%
        gtsummary::tbl_regression(
          exponentiate = exponentiate,
          hide_n = TRUE,
          tidy_fun = tidy_fun
        )
    }
  }
  
  df %>%
    gtsummary::tbl_strata(
      strata = !!strata,
      .tbl_fun = ~ .tbl_fun(.x) %>%
        gtsummary::bold_labels() %>%
        gtsummary::add_global_p() %>%
        gtsummary::add_nevent(location = "level") %>%
        gtsummary::add_n(location = "level") %>%
        gtsummary::modify_column_merge(
          pattern = "{estimate} [{conf.low} - {conf.high}]",
          rows = !is.na(estimate)
        ) %>%
        gtsummary::modify_header(label = "", estimate = "**HR [95% CI]**") %>%
        gtsummary::bold_p(t = 0.05)
    )
}


#' Compute global tests for clinical, DNA, and RNA feature sets in survival 
#' analysis
#'
#' This function performs global association tests between groups of predictors
#' and a survival outcome using the `globaltest` package. It computes separate
#' global tests for clinical variables, gene-level DNA features, pathway-level
#' DNA features, gene-level RNA features, and pathway-level RNA features.
#'
#' Predictor groups are defined by explicitly provided clinical variables and by
#' column name prefixes for DNA and RNA features.
#'
#' @param data_gene A data.frame or tibble containing patient-level clinical,
#'   gene-level DNA, gene-level RNA, and survival variables.
#'
#' @param data_pathway A data.frame or tibble containing patient-level
#'   pathway-level DNA, pathway-level RNA, and survival variables.
#'
#' @param time A character string specifying the name of the survival time
#'   variable.
#'
#' @param event A character string specifying the name of the event indicator
#'   variable.
#'
#' @param clin A character vector containing the names of the clinical
#'   variables to include in the clinical global test.
#'
#' @param dna_pattern A character string used to identify gene-level DNA
#'   feature columns in `data_gene`, typically via a common prefix.
#'
#' @param dna_pathway_pattern A character string used to identify pathway-level
#'   DNA feature columns in `data_pathway`, typically via a common prefix.
#'
#' @param rna_pattern A character string used to identify gene-level RNA
#'   feature columns in `data_gene`, typically via a common prefix.
#'
#' @param rna_pathway_pattern A character string used to identify pathway-level
#'   RNA feature columns in `data_pathway`, typically via a common prefix.
#'
#' @return A tibble with one row per feature set tested, containing:
#' \itemize{
#'   \item `data_type`: the feature group tested,
#'   \item `x_cov`: the number of covariates included in the test,
#'   \item `statistic`: the observed global test statistic,
#'   \item `expected`: the expected value of the statistic under the null,
#'   \item `std_dev`: the standard deviation under the null,
#'   \item `p_value`: the p-value of the global test.
#' }
#'
#' @export
compute_global_test <- function(data_gene,
                                data_pathway, 
                                time, 
                                event, 
                                clin, 
                                dna_pattern, 
                                dna_pathway_pattern, 
                                rna_pattern, 
                                rna_pathway_pattern){
  
  dna_vars <- data_gene %>% 
    dplyr::select(starts_with(dna_pattern))%>%
    colnames()
  dna_pathway_vars <- data_pathway %>% 
    dplyr::select(starts_with(dna_pathway_pattern))%>%
    colnames()
  rna_vars <- data_gene %>% 
    dplyr::select(starts_with(rna_pattern))%>%
    colnames()
  rna_pathway_vars <- data_pathway %>% 
    dplyr::select(starts_with(rna_pathway_pattern))%>%
    colnames()
  
  .gt_surv <- function(.data, X) {
    Y <- survival::Surv(.data[[time]], .data[[event]])
    X <- .data %>% 
      dplyr::select(all_of(X))%>%
      as.matrix()
    janitor::clean_names(data.frame(globaltest::gt(Y, X)@result))%>%
      dplyr::relocate(x_cov, statistic, expected, std_dev, p_value)
  }
  
  gt_clin <- tibble(
    data_type = "Clinical",
    .gt_surv(data_gene, clin)
    )
  gt_dna <- tibble(
    data_type = "DNA (gene-level)",
    .gt_surv(data_gene, dna_vars)
    )
  gt_dna_pathway <- tibble(
    data_type = "DNA (pathway-level)",
    .gt_surv(data_pathway, dna_pathway_vars)
    )
  gt_rna <- tibble(
    data_type = "RNA (gene-level)",
    .gt_surv(data_gene, rna_vars)
    )
  gt_rna_pathway <- tibble(
    data_type = "RNA (pathway-level)",
    .gt_surv(data_pathway, rna_pathway_vars)
    )
  
  list(gt_clin, 
       gt_dna, 
       gt_dna_pathway,
       gt_rna, 
       gt_rna_pathway)%>%
    bind_rows()
}
