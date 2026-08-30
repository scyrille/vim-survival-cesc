
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
                      xlim = c(0,25),  
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

plot_surv_table <- function(
    formula,
    data,
    title = "",
    legend = "top",
    legend.title = "",
    xlab = "Time in months",
    ylab = "Progression-free survival",
    xlim = c(0, 32),
    break.x.by = 6,
    annot_table = TRUE,
    surv_times = c(12, 24),
    surv_label = "PFS",
    annot_table_x = 0.02,
    annot_table_y = 0.05,
    annot_table_width = 0.88,
    annot_table_height = 0.22,
    risk.table = "nrisk_cumcensor",
    risk.table.y.text = TRUE,
    risk.table.fontsize = 2.8,
    tables.height = 0.16,
    risk.table.col = "strata",
    conf.int = FALSE,
    surv.scale = "percent",
    axes.offset = TRUE,
    palette, 
    pval = FALSE
) {
  
  surv_times <- sort(unique(surv_times))
  
  surv_time_labels <- format(
    surv_times,
    trim = TRUE,
    scientific = FALSE
  )
  
  fit <- survminer::surv_fit(
    formula = formula,
    data = data,
    match.fd = FALSE
  )
  
  has_strata <- !is.null(fit$strata)
  
  if (has_strata) {
    names(fit$strata) <- sub(
      "^.*?=",
      "",
      names(fit$strata)
    )
  }
  
  # Survival curve 
  
  gg <- survminer::ggsurvplot(
    fit = fit,
    title = title,
    legend = if (has_strata) legend else "none",
    legend.title = legend.title,
    xlab = xlab,
    ylab = ylab,
    xlim = xlim,
    break.x.by = break.x.by,
    conf.int = conf.int,
    risk.table = risk.table,
    risk.table.y.text = risk.table.y.text,
    risk.table.fontsize = risk.table.fontsize,
    risk.table.col = risk.table.col,
    tables.height = tables.height,
    censor.shape = 124,
    censor.size = 1.5,
    size = 0.8,
    palette = palette, 
    ggtheme = ggplot2::theme_classic(base_size = 8) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        
        axis.text = ggplot2::element_text(
          size = 7.5,
          colour = "black"
        ),
        axis.title.x = ggplot2::element_text(
          size = 8,
          face = "bold", 
          margin = ggplot2::margin(t = 4)
        ),
        axis.title.y = ggplot2::element_text(
          size = 8,
          face = "bold",
          margin = ggplot2::margin(r = 4)
        ),
        
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.width = grid::unit(0.8, "cm"),
        legend.spacing.x = grid::unit(0.15, "cm"),
        
        plot.margin = ggplot2::margin(
          t = 3,
          r = 5,
          b = 3,
          l = 3
        )
      ),
    
    tables.theme = survminer::theme_cleantable() +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        plot.background = ggplot2::element_rect(
          fill = "white",
          colour = NA
        ),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          size = 7.5,
          face = "bold"
        )
      ),
    surv.scale = surv.scale,
    axes.offset = axes.offset,
    pval = pval
  )
  
  if (!isTRUE(annot_table)) {
    return(gg)
  }
  
  
  # N and N events 
  fit_table <- summary(fit)$table
  
  if (is.null(dim(fit_table))) {
    
    n_column <- intersect(
      c("records", "n.max", "n.start"),
      names(fit_table)
    )[1]
    
    if (is.na(n_column)) {
      stop("Impossible d'identifier le nombre de patients.")
    }
    
    base_df <- tibble::tibble(
      Strata = "Overall",
      N = as.integer(
        fit_table[[n_column]]
      ),
      Events = as.integer(
        fit_table[["events"]]
      )
    )
    
  } else {
    
    base_df <- fit_table %>%
      as.data.frame() %>%
      tibble::rownames_to_column(
        var = "Strata"
      )
    
    n_column <- intersect(
      c("records", "n.max", "n.start"),
      names(base_df)
    )[1]
    
    if (is.na(n_column)) {
      stop("Impossible d'identifier le nombre de patients.")
    }
    
    base_df <- base_df %>%
      dplyr::transmute(
        Strata = sub(
          "^.*=",
          "",
          Strata
        ),
        `No. of\npatients` = as.integer(
          .data[[n_column]]
        ),
        `No. of\nevents` = as.integer(
          .data[["events"]]
        )
      )
  }
  
  # Survival rates 
  
  rate_column_names <- paste0(
    surv_time_labels,
    "-month ",
    surv_label,
    "\n[95% CI]"
  )
  
  rates_long <- dplyr::bind_rows(
    lapply(
      seq_along(surv_times),
      function(i) {
        
        requested_time <- surv_times[i]
        
        surv_at_time <- summary(
          fit,
          times = requested_time,
          extend = TRUE
        )
        
        if (is.null(surv_at_time$strata)) {
          
          strata_at_time <- rep(
            "Overall",
            length(surv_at_time$surv)
          )
          
        } else {
          
          strata_at_time <- sub(
            "^.*=",
            "",
            as.character(surv_at_time$strata)
          )
        }
        
        missing_estimate <- (
          is.na(surv_at_time$surv) |
            is.na(surv_at_time$lower) |
            is.na(surv_at_time$upper)
        )
        
        formatted_rate <- ifelse(
          missing_estimate,
          "NE",
          sprintf(
            "%.1f%% [%.1f–%.1f]",
            100 * surv_at_time$surv,
            100 * surv_at_time$lower,
            100 * surv_at_time$upper
          )
        )
        
        tibble::tibble(
          Strata = strata_at_time,
          rate_column = rate_column_names[i],
          rate = formatted_rate
        )
      }
    )
  )
  
  rates_df <- rates_long %>%
    dplyr::mutate(
      rate_column = factor(
        rate_column,
        levels = rate_column_names
      )
    ) %>%
    tidyr::pivot_wider(
      id_cols = Strata,
      names_from = rate_column,
      values_from = rate,
      names_sort = FALSE
    )
  
  # Final table
  
  summary_df <- base_df %>%
    dplyr::left_join(
      rates_df,
      by = "Strata"
    ) %>%
    dplyr::select(
      Strata,
      `No. of\npatients`,
      `No. of\nevents`,
      dplyr::all_of(rate_column_names)
    )
  
  # Colors 
  plot_build <- ggplot2::ggplot_build(
    gg$plot
  )
  
  palette <- unique(
    plot_build$data[[1]]$colour
  )
  
  palette <- palette[
    !is.na(palette)
  ]
  
  if (length(palette) == 0L) {
    palette <- "black"
  }
  
  # Table construction 
  summary_grob <- gridExtra::tableGrob(
    summary_df,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = 7,
      core = list(
        fg_params = list(
          fontsize = 7
        ),
        bg_params = list(
          fill = "white",
          col = NA
        ),
        padding = grid::unit(
          c(1, 1.5),
          "mm"
        )
      ),
      colhead = list(
        fg_params = list(
          col = "black",
          fontface = "bold",
          fontsize = 7
        ),
        bg_params = list(
          fill = "white",
          col = NA
        ),
        padding = grid::unit(
          c(1, 1.5),
          "mm"
        )
      )
    )
  )
  
  # Column widths 
  summary_grob$widths <- grid::unit(
    c(
      1.20,                         # Strata
      1,                            # N
      1,                            # Events
      rep(1.8, length(surv_times))  # Survival rates
    ),
    "null"
  )
  
  # Colour all cells 
  table_layout <- summary_grob$layout
  
  body_cells <- which(
    table_layout$name == "core-fg"
  )
  
  body_rows <- sort(
    unique(table_layout$t[body_cells])
  )
  
  for (i in seq_along(body_rows)) {
    
    row_cells <- body_cells[
      table_layout$t[body_cells] == body_rows[i]
    ]
    
    row_color <- palette[
      (i - 1) %% length(palette) + 1
    ]
    
    for (cell in row_cells) {
      
      is_strata_column <- (
        table_layout$l[cell] == 1
      )
      
      summary_grob$grobs[[cell]]$gp <-
        grid::gpar(
          col = row_color,
          fontsize = 7#,
          # fontface = if (
          #   is_strata_column
          # ) {
          #   "bold"
          # } else {
          #   "plain"
          # }
        )
    }
  }
  
  # Table height
  if (is.null(annot_table_height)) {
    
    annot_table_height <- max(
      0.5,
      min(
        0.40,
        0.2 * (nrow(summary_df) + 1)
      )
    )
  }
  
  # Insert table
  
  plot_x_limits <- if (is.null(xlim)) {
    range(fit$time, na.rm = TRUE)
  } else {
    xlim
  }
  
  plot_x_range <- diff(plot_x_limits)
  
  table_xmin <- plot_x_limits[1] +
    annot_table_x * plot_x_range
  
  table_xmax <- plot_x_limits[1] +
    min(
      1,
      annot_table_x + annot_table_width
    ) * plot_x_range
  
  table_ymin <- annot_table_y
  
  table_ymax <- min(
    0.98,
    annot_table_y + annot_table_height
  )
  
  # Add the table to the ggsurvplot 
  
  gg$plot <- gg$plot +
    ggplot2::annotation_custom(
      grob = summary_grob,
      xmin = table_xmin,
      xmax = table_xmax,
      ymin = table_ymin,
      ymax = table_ymax
    )
  
  if (!is.null(gg$table)) {
    
    gg$table <- gg$table +
      ggplot2::labs(
        x = "",
        title = "No. at risk (censored)"
      )
    
    gg$table$theme$axis.title.y <- ggplot2::element_blank()
    gg$table$theme$axis.text.x <- ggplot2::element_blank()
    gg$table$theme$axis.ticks.x <- ggplot2::element_blank()
    gg$table$theme$axis.title.x <- ggplot2::element_blank()
    gg$table$theme$legend.position <- "none"
    
    risk_table_height <- max(
      tables.height,
      min(
        0.35,
        0.10 + 0.04 * nrow(summary_df)
      )
    )
  
    aligned_plots <- cowplot::align_plots(
      gg$plot,
      gg$table,
      align = "v",
      axis = "lr"
    )
    
    final_plot <- cowplot::plot_grid(
      plotlist = aligned_plots,
      ncol = 1,
      rel_heights = c(
        1 - risk_table_height,
        risk_table_height
      )
    )
    
  } else {
    
    final_plot <- gg$plot
  }
  
  final_grob <- ggplot2::ggplotGrob(
    final_plot
  )

  list(
    plot = final_plot,
    grob = final_grob,
    ggs = gg,
    table = summary_df,
    table_grob = summary_grob,
    fit = fit
  )
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
tbl_cox <- function(data,
                    time,
                    event,
                    strata = NULL,
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
          include = all_of(covars),
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
  
  if (!is.null(strata)) {
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
  } else {
    df %>% 
      .tbl_fun() %>%
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
  }
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
