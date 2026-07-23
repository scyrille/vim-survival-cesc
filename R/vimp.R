
#'@description Model-agnostic, algorithm-free variable importance for 
#'survival analysis (`survML`)

# ----------------------------- Documentation ------------------------------#

#'This file contains functions to compute variable importance estimates via
#'an algorithm-agnostic estimator of variable importance (`survML`). 

# ----------------------------- Dependencies -------------------------------#

library(dplyr)
library(stringr)
library(gt)
library(ggplot2)
library(geomtextpath)
library(egg)
library(survival)  
library(survML)  
library(SuperLearner)
library(glmnet)
library(ranger)
library(xgboost)

# ------------------------------ Functions ---------------------------------#

# Statistical helpers -----------------------------------------------------

# Algorithm-agnostic variable importance via survML 

## Variable importance relative to all features 
#'@param time Numeric vector of observed follow-up times. 
#'@param event Numeric vector of status indicators of whether an event was 
#'observed or not.
#'@param X Data.frame of observed covariate values.
#'@param feature_groups named vector describing the grouping of the 
#'covariates.  
#'@param landmark_times Numeric vector giving landmark times at which to 
#'estimate variable importance. 
#'@param type Type of variable importance to compute. Options include 
#'"accuracy", "AUC", "Brier", "R-squared", "C-index", and 
#'"survival_time_MSE". 
#'@param bin_size Size of time grid used for estimating nuisance parameters
#'@param SL.library Character vector of prediction algorithms included 
#'in the SuperLearner package.
#'@param nfolds Number of cross-validation folds. 
#'@param cf_fold_num Number of cross-fitting folds. 

compute_vimp_survML_full <- function(time, 
                                     event, 
                                     X, 
                                     feature_groups, 
                                     landmark_times = c(12, 24, 36), 
                                     type           = "AUC", 
                                     bin_size       = 0.5,
                                     SL.library     = c("SL.mean",
                                                        "SL.ranger",
                                                        "SL.glmnet",
                                                        "SL.xgboost"), 
                                     nfolds         = 5,
                                     cf_fold_num    = 5,
                                     seed           = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  output <- vector("list", length = length(feature_groups))
  
  # Estimate VIMP for the first feature group 
  parse_indx <- function(x) as.numeric(stringr::str_split(x, "\\s*,\\s*")[[1]])
  
  message(paste("Estimating", type, "importance of",names(feature_groups)[1],"..."))
  output[[1]] <- survML::vim(
    type = type,
    time = time,
    event = event,
    X = X,
    landmark_times = landmark_times,
    restriction_time = max(time[event == 1]), 
    large_feature_vector = 1:ncol(X),
    small_feature_vector = (1:ncol(X))[-parse_indx(feature_groups[[1]])],
    conditional_surv_generator_control = list(SL.library = SL.library,
                                              bin_size = bin_size,
                                              V = nfolds),
    large_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    small_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    cf_fold_num = cf_fold_num,
    sample_split = TRUE,
    scale_est = TRUE)
  
  saved_conditional_surv_preds <- output[[1]]$conditional_surv_preds
  saved_large_oracle_preds <- output[[1]]$large_oracle_preds
  saved_folds <- output[[1]]$folds
  saved_approx_times <- output[[1]]$approx_times
  pooled_output <- output$result   
  
  # Iterate over other feature groups
  output[2:length(feature_groups)] <- purrr::map2(
    feature_groups[2:length(feature_groups)],
    names(feature_groups)[2:length(feature_groups)],
    ~ {
      message(paste("Estimating", type, "importance of", .y, "..."))
      
      survML::vim(
        type = type,
        time = time,
        event = event,
        X = X,
        landmark_times = landmark_times,
        restriction_time = max(time[event == 1]),
        approx_times = saved_approx_times,
        large_feature_vector = 1:ncol(X),
        small_feature_vector = (1:ncol(X))[-parse_indx(.x)],
        conditional_surv_preds = saved_conditional_surv_preds,
        large_oracle_preds = saved_large_oracle_preds,
        cf_folds = saved_folds$cf_folds,
        ss_folds = saved_folds$ss_folds,
        small_oracle_generator_control = list(SL.library = SL.library,
                                              V = nfolds),
        sample_split = TRUE,
        scale_est = TRUE
      )
    }
  )
  
  output %>%
    purrr::imap(
      ~ {
        .x$result <- .x$result %>%
          dplyr::mutate(
            indx = feature_groups[[.y]],
            variable = names(feature_groups)[[.y]]
          )
        .x
      }
    ) %>% 
    purrr::set_names(names(feature_groups))
}


## Variable importance relative to base model

#'@param time Numeric vector of observed follow-up times. 
#'@param event Numeric vector of status indicators of whether an event was 
#'observed or not.
#'@param X Data.frame of observed covariate values.
#'@param base_features A named vector describing the baseline set of 
#'covariates.  
#'@param feature_groups named vector describing the grouping of the 
#'covariates.  
#'@param landmark_times Numeric vector giving landmark times at which to 
#'estimate variable importance. 
#'@param type Type of variable importance to compute. Options include 
#'"accuracy", "AUC", "Brier", "R-squared", "C-index", and 
#'"survival_time_MSE". 
#'@param bin_size Size of time grid used for estimating nuisance parameters
#'@param SL.library Character vector of prediction algorithms included 
#'in the SuperLearner package.
#'@param nfolds Number of cross-validation folds. 
#'@param cf_fold_num Number of cross-fitting folds. 

compute_vimp_survML_base <- function(time, 
                                     event, 
                                     X, 
                                     base_features,
                                     feature_groups, 
                                     landmark_times = 24, 
                                     type           = "AUC", 
                                     bin_size       = 0.1,
                                     SL.library, 
                                     nfolds         = 5,
                                     cf_fold_num    = 5,
                                     seed           = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  #'Estimating variable importance relative to all features 
  #'This step is necessary to estimate the conditional survival function
  #'estimates given ALL features 
  
  message("Estimating the conditional survival function estimates given all features...")
  output_base <- vim(
    type = type,
    time = time,
    event = event,
    X = X,
    landmark_times = landmark_times,
    large_feature_vector = 1:ncol(X),
    small_feature_vector = (1:ncol(X))[-as.numeric(base_features)],
    conditional_surv_generator_control = list(SL.library = SL.library,
                                              V = nfolds,
                                              bin_size = bin_size),
    large_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    small_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    cf_fold_num = cf_fold_num,
    sample_split = TRUE,
    scale_est = TRUE)
  
  saved_conditional_surv_preds <- output_base$conditional_surv_preds
  saved_folds <- output_base$folds
  saved_approx_times <- output_base$approx_times
  
  #'Estimating variable importance relative to the baseline set of variables
  feature_groups <- 
    lapply(feature_groups, 
           function(x) paste0(c(x, base_features), 
                              collapse = ","))%>% unlist()
  output <- vector("list", length = length(feature_groups))
  
  # Estimate VIMP for the first feature group 
  parse_indx <- function(x) as.numeric(stringr::str_split(x, "\\s*,\\s*")[[1]])
  
  message(paste("Estimating", type,"importance of", names(feature_groups[1]),"..."))
  
  output[[1]] <- vim(
    type = type,
    time = time,
    event = event,
    X = X,
    landmark_times = landmark_times,
    approx_times = saved_approx_times,
    large_feature_vector = parse_indx(feature_groups[[1]]),
    small_feature_vector =  as.numeric(base_features),
    conditional_surv_preds = saved_conditional_surv_preds,
    large_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    small_oracle_generator_control = list(SL.library = SL.library,
                                          V = nfolds),
    cf_folds = saved_folds$cf_folds,
    ss_folds = saved_folds$ss_folds,
    sample_split = TRUE,
    scale_est = TRUE)
  saved_small_oracle_preds <- output[[1]]$small_oracle_preds
  
  # Iterate over other feature groups
  output[2:length(feature_groups)] <- purrr::map2(
    feature_groups[2:length(feature_groups)],
    names(feature_groups)[2:length(feature_groups)],
    ~ {
      message(paste("Estimating", type, "importance of", .y, "..."))
      
      survML::vim(
        type = type,
        time = time,
        event = event,
        X = X,
        landmark_times = landmark_times,
        approx_times = saved_approx_times,
        large_feature_vector = parse_indx(.x),
        small_feature_vector = as.numeric(base_features),
        conditional_surv_preds = saved_conditional_surv_preds,
        small_oracle_preds = saved_small_oracle_preds,
        large_oracle_generator_control = list(SL.library = SL.library, 
                                              V = nfolds),
        cf_folds = saved_folds$cf_folds,
        ss_folds = saved_folds$ss_folds,
        sample_split = TRUE,
        scale_est = TRUE
      )
    }
  )
  
  output %>%
    purrr::imap(
      ~ {
        .x$result <- .x$result %>%
          dplyr::mutate(
            indx = feature_groups[[.y]],
            variable = names(feature_groups)[[.y]]
          )
        .x
      }
    ) %>% 
    purrr::set_names(names(feature_groups))
}

#' Build the Super Learner library for VIM relative to base model 
make_SL_library <- function(
    ranger_trees = c(250, 500, 1000),
    xgb_trees = c(250, 500, 1000),
    xgb_depth = c(1, 2, 4)
) {

  # Ridge (glmnet)
  ridge_learner <- SuperLearner::create.Learner(
    "SL.glmnet", 
    detailed_names = TRUE, 
    params = list(
      alpha = 0,
      standardize = TRUE
    ),
    name_prefix = "ridge",
    env = .GlobalEnv
  )
  
  # Ranger 
  ranger_learner <- SuperLearner::create.Learner(
    "SL.ranger",
    tune = list(
      num.trees = ranger_trees
    ),
    name_prefix = "ranger",
    env = .GlobalEnv
  )
  
  # XGBoost 
  xgboost_learner <- SuperLearner::create.Learner(
    "SL.xgboost",
    tune = list(
      ntrees = xgb_trees,
      max_depth = xgb_depth
    ),
    name_prefix = "xgboost",
    env = .GlobalEnv
  )
  
  # Final library 
  c("SL.mean",
    ridge_learner$names,
    ranger_learner$names,
    xgboost_learner$names)
  
}

#' Input of `compute_vimp_survML_full` function
make_input_vimp_survML_full <- function(data, 
                                        var_clin,
                                        dna_prefix, 
                                        rna_prefix){
  
  X <- data %>%
    dplyr::select(all_of(var_clin), 
                  starts_with(unlist(strsplit(dna_prefix, split = "[|]"))),
                  starts_with(unlist(strsplit(rna_prefix, split = "[|]"))))
  
  feature_groups <- paste0(which(names(X) %in% names(X)))%>%
    set_names(var_label(X, unlist = T))

  list(time = data$time,
       event = data$event,
       X = X, 
       feature_groups = feature_groups)
}

#' Input of `compute_vimp_survML_base` function
make_input_vimp_survML_base <- function(data, 
                                        var_clin,
                                        dna_prefix, 
                                        rna_prefix){
  
  X <- data %>%
    dplyr::select(all_of(var_clin), 
                  starts_with(unlist(strsplit(dna_prefix, split = "[|]"))),
                  starts_with(unlist(strsplit(rna_prefix, split = "[|]"))))
  
  base_features <- paste0(which(names(X) %in% var_clin))%>%
    set_names(var_label(X %>% dplyr::select(all_of(var_clin)),
                        unlist = T))
  
  pathway_features <- paste0(which(!names(X) %in% var_clin))%>%
    set_names(var_label(X %>% dplyr::select(-all_of(var_clin)),
                        unlist = T))
  
  group_DNA_pathway_features <- get_pathways_process() %>%
    dplyr::filter(data_type == "DNA_pathways") %>%
    dplyr::select(process, variable) %>%
    dplyr::transmute(process = paste0(process, " (DNA)"), variable) %>%
    dplyr::group_by(process) %>%
    dplyr::summarise(variable = list(variable), .groups = "drop") %>%
    tibble::deframe()%>%
    purrr::imap_chr(~paste(which(as.character(var_label(X)) %in% .x), collapse = ",")
    )
  
  group_RNA_pathway_features <- get_pathways_process() %>%
    dplyr::filter(data_type == "RNA_pathways") %>%
    dplyr::transmute(process = paste0(process, " (RNA)"), variable) %>%
    dplyr::group_by(process) %>%
    dplyr::summarise(variable = list(variable), .groups = "drop") %>%
    tibble::deframe()%>%
    purrr::imap_chr(~paste(which(as.character(var_label(X)) %in% .x), collapse = ",")
  )
  
  feature_groups <- c(pathway_features, 
                      group_DNA_pathway_features,
                      group_RNA_pathway_features)
  
  list(time = data$time,
       event = data$event,
       X = X, 
       base_features = base_features, 
       feature_groups = feature_groups)
}

# Plotting and table helpers ----------------------------------------------

#' Get the variable importance estimates.

get_vimp_est <- function(fit, 
                         landmark_time = 24, 
                         method = "survML"){
  
  if (method == "survML"){
    
    vim_est <- do.call(rbind, lapply(fit, function(x) x$result))%>%
      {rownames(.) <- NULL; .}%>%
      left_join(get_pathways_process(), by = "variable")%>%
      mutate(data_type = case_when(grepl("\\(DNA\\)", variable)~"DNA_processes",
                                   grepl("\\(RNA\\)", variable)~"RNA_processes",
                                   TRUE~data_type),
             variable = gsub(" \\(DNA\\)| \\(RNA\\)", "", variable),
             label = ifelse(is.na(label), variable, label)) %>%
      dplyr::filter(!is.na(data_type) & landmark_time %in% !!landmark_time)%>%
      dplyr::group_by(data_type) %>%
      dplyr::mutate(p_adj = p.adjust(p, method = "BH") ) %>%
      dplyr::ungroup()
    
  } else if (method == "survex"){
    
    dna_pathway <- fit[["model_parts_pathways"]]$result %>%
      dplyr::select(`_times_`, all_of(as.character(dna_labels)), `_permutation_`)%>%
      dplyr::mutate(data_type = "DNA_pathways")%>%
      dplyr::filter(`_times_`%in% !!landmark_time & `_permutation_` == 0)%>%
      pivot_longer(-c(`_times_`,`_permutation_`, data_type), 
                   names_to = "label", 
                   values_to = "est")
    
    dna_process <- fit[["model_parts_pathways_groups"]]$result %>%
      dplyr::select(`_times_`, matches("\\(DNA\\)"), `_permutation_`)%>%
      dplyr::mutate(data_type = "DNA_processes")%>%
      dplyr::filter(`_times_`%in% !!landmark_time & `_permutation_` == 0)%>%
      pivot_longer(-c(`_times_`,`_permutation_`, data_type), 
                   names_to = "label", 
                   values_to = "est")
    
    rna_pathway <- fit[["model_parts_pathways"]]$result %>%
      dplyr::select(`_times_`, all_of(as.character(rna_labels)), `_permutation_`)%>%
      dplyr::mutate(data_type = "RNA_pathways")%>%
      dplyr::filter(`_times_`%in% !!landmark_time & `_permutation_` == 0)%>%
      pivot_longer(-c(`_times_`,`_permutation_`, data_type), 
                   names_to = "label", 
                   values_to = "est")
    
    rna_process <- fit[["model_parts_pathways_groups"]]$result %>%
      dplyr::select(`_times_`, matches("\\(RNA\\)"), `_permutation_`)%>%
      dplyr::mutate(data_type = "RNA_processes")%>%
      dplyr::filter(`_times_` %in% !!landmark_time & `_permutation_` == 0)%>%
      pivot_longer(-c(`_times_`,`_permutation_`, data_type), 
                   names_to = "label", 
                   values_to = "est")
    
    vim_est <- list(dna_pathway, dna_process, rna_pathway, rna_process)%>%
      purrr::reduce(full_join, by = names(dna_pathway))%>%
      dplyr::rename(landmark_time = `_times_`)%>%
      dplyr::mutate(label = gsub(" \\(DNA\\)| \\(RNA\\)","", label),
                    est = ifelse(est<0, 0, est))%>%
      left_join(get_pathways_process(), by = c("label","data_type"))
    
  }
  return(vim_est)
}

## Plot VIM estimates
plot_vimp_est <- function(
    vims,
    method = "survML",
    xlab = NULL,
    ylab = "Variable",
    type = c("barplot", "dotplot"),
    process_panel = FALSE,
    fill_by = NULL,
    fill_values = NULL,
    fill_label = NULL,
    p_text_color = "#F8766D",
    vline0 = TRUE,
    vline_color = "red",
    vline_linetype = "dotted",
    vline_linewidth = 0.7,
    limits = c(0, 0.5),
    legend.position = "top",
    legend.direction = "horizontal"
) {
  
  type <- match.arg(type)
  
  if (is.null(xlab)) {
    xlab <- if (method == "survML") {
      "VIM ± 95% CI"
    } else {
      "VIM"
    }
  }
  
  label_col <- "label"
  
  required_cols <- if (method == "survML") {
    c("est", "cil_1sided", "ciu", "p", label_col)
  } else {
    c("est", label_col)
  }
  
  missing_cols <- setdiff(required_cols, names(vims))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  df <- vims %>%
    dplyr::mutate(
      ord_group = factor(.data[[label_col]])
    )
  
  if (method == "survML") {
    df <- df %>%
      dplyr::mutate(
        est_lab = paste0(
          sprintf("%.2f", est),
          " [",
          sprintf("%.2f", cil_1sided),
          "-",
          sprintf("%.2f", ciu),
          "] ",
          gtools::stars.pval(p)
        ),
        p_val = gtools::stars.pval(p)
      )
  }
  
  # Fill variable
  if (is.null(fill_by)) {
    aes_fill <- NULL
    
  } else {
    fill_by <- rlang::as_name(rlang::ensym(fill_by))
    aes_fill <- fill_by
    
    if (is.null(fill_label)) {
      fill_label <- dplyr::case_when(
        fill_by == "cohort" ~ "Cohort",
        TRUE ~ fill_by
      )
    }
    
    if (is.null(fill_values) && fill_by == "cohort") {
      fill_values <- c(
        "TCGA-CESC" = "#377eb8",
        "Bio-RAIDs" = "#1b9e77"
      )
    }
  }
  
  base_theme <- function(
    plot,
    legend.position,
    legend.direction
  ) {
    plot +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(size = 18),
        axis.text.y = ggplot2::element_text(size = 18),
        legend.position = legend.position,
        legend.direction = legend.direction,
        legend.title = ggplot2::element_text(size = 18),
        legend.text = ggplot2::element_text(size = 18),
        strip.text = ggplot2::element_text(
          face = "bold",
          size = 18
        ),
        strip.placement = "outside",
        strip.background = ggplot2::element_rect(
          fill = "grey90",
          colour = NA
        ),
        plot.margin = ggplot2::margin(8, 8, 8, 8),
        panel.spacing = grid::unit(0.8, "lines"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          colour = "grey80",
          fill = NA,
          linewidth = 0.4
        )
      )
  }
  
  add_common_elements <- function(
    plot,
    facet,
    ncol,
    nrow,
    legend.position,
    legend.direction,
    scale_type
  ) {
    
    if (isTRUE(vline0)) {
      plot <- plot +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = vline_linetype,
          colour = vline_color,
          linewidth = vline_linewidth
        )
    }
    
    if (!is.null(aes_fill) && !is.null(fill_values)) {
      if (scale_type == "colour") {
        plot <- plot +
          ggplot2::scale_colour_manual(values = fill_values)
      } else {
        plot <- plot +
          ggplot2::scale_fill_manual(values = fill_values)
      }
    }
    
    if (isTRUE(facet)) {
      plot <- plot +
        ggplot2::facet_wrap(
          ~ process,
          nrow = nrow,
          ncol = ncol,
          scales = "free_y"
        )
    }
    
    base_theme(
      plot,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
  
  .dotplot_vimp <- function(
    df,
    facet = FALSE,
    ncol = NULL,
    nrow = NULL,
    legend.position,
    legend.direction
  ) {
    
    dodge <- ggplot2::position_dodge(
      width = 0.8,
      preserve = "single"
    )
    
    plot <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["est"]],
        y = .data[["ord_group"]],
        colour = if (!is.null(aes_fill)) {
          .data[[aes_fill]]
        } else {
          NULL
        }
      )
    ) +
      ggplot2::geom_point(
        position = dodge,
        size = 3
      )
    
    if (method == "survML") {
      plot <- plot +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = .data[["cil_1sided"]],
            xmax = .data[["ciu"]]
          ),
          position = dodge,
          linewidth = 0.9
        ) #+
        # ggplot2::geom_text(
        #   ggplot2::aes(
        #     x = .data[["ciu"]],
        #     label = .data[["p_val"]]
        #   ),
        #   hjust = -0.2,
        #   size = 6,
        #   colour = p_text_color
        # )
    }
    
    plot <- plot +
      ggplot2::scale_x_continuous(limits = limits) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::labs(
        x = xlab,
        y = ylab,
        colour = fill_label
      )
    
    add_common_elements(
      plot = plot,
      facet = facet,
      ncol = ncol,
      nrow = nrow,
      legend.position = legend.position,
      legend.direction = legend.direction,
      scale_type = "colour"
    )
  }
  
  .barplot_vimp <- function(
    df,
    facet = FALSE,
    ncol = NULL,
    nrow = NULL,
    legend.position,
    legend.direction
  ) {
    
    dodge <- ggplot2::position_dodge(
      width = 0.8,
      preserve = "single"
    )
    
    plot <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = .data[["est"]],
        y = .data[["ord_group"]],
        fill = if (!is.null(aes_fill)) {
          .data[[aes_fill]]
        } else {
          NULL
        }
      )
    ) +
      ggplot2::geom_col(
        width = 0.55,
        position = dodge,
        colour = NA,
        alpha = 0.95
      )
    
    if (method == "survML") {
      plot <- plot +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = .data[["cil_1sided"]],
            xmax = .data[["ciu"]]
          ),
          position = dodge,
          linewidth = 0.4,
          colour = "black"
        ) #+
        # ggplot2::geom_text(
        #   ggplot2::aes(
        #     x = .data[["ciu"]],
        #     label = .data[["p_val"]]
        #   ),
        #   hjust = -0.2,
        #   size = 6,
        #   colour = p_text_color
        # )
    }
    
    plot <- plot +
      ggplot2::scale_x_continuous(limits = limits) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::labs(
        x = xlab,
        y = ylab,
        fill = fill_label
      )
    
    add_common_elements(
      plot = plot,
      facet = facet,
      ncol = ncol,
      nrow = nrow,
      legend.position = legend.position,
      legend.direction = legend.direction,
      scale_type = "fill"
    )
  }
  
  fun <- switch(
    type,
    dotplot = .dotplot_vimp,
    barplot = .barplot_vimp
  )
  
  if (isTRUE(process_panel)) {
    
    if (!"process" %in% names(df)) {
      stop(
        "`process_panel = TRUE` requires a `process` column.",
        call. = FALSE
      )
    }
    
    df_p1 <- df %>%
      dplyr::filter(
        process %in% c(
          "Cellular component",
          "DNA damage",
          "Pathway"
        )
      )
    
    df_p2 <- df %>%
      dplyr::filter(
        process %in% c(
          "Proliferation",
          "Development",
          "Epigenetics"
        )
      )
    
    df_p3 <- df %>%
      dplyr::filter(
        process %in% c(
          "Metabolic",
          "Immune"
        )
      )
    
    df_p4 <- df %>%
      dplyr::filter(process == "Signaling")
    
    p1 <- fun(
      df_p1,
      facet = TRUE,
      nrow = 3,
      ncol = 1,
      legend.position = legend.position,
      legend.direction = legend.direction
    ) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      )
    
    p2 <- fun(
      df_p2,
      facet = TRUE,
      nrow = 2,
      ncol = 1,
      legend.position = "none",
      legend.direction = legend.direction
    ) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      )
    
    p3 <- fun(
      df_p3,
      facet = TRUE,
      nrow = 2,
      ncol = 1,
      legend.position = "none",
      legend.direction = legend.direction
    )
    
    p4 <- fun(
      df_p4,
      facet = TRUE,
      nrow = 1,
      ncol = 1,
      legend.position = "none",
      legend.direction = legend.direction
    )
    
    egg::ggarrange(
      p1,
      p2,
      p3,
      p4,
      ncol = 2,
      nrow = 2,
      widths = c(1, 1),
      heights = c(1, 1)
    )
    
  } else {
    
    fun(
      df = df,
      facet = FALSE,
      nrow = NULL,
      ncol = NULL,
      legend.position = legend.position,
      legend.direction = legend.direction
    )
  }
}

# Summary tables of all VIM estimates
tbl_vimp <- function(
    vims,
    compare = TRUE
) {
  
  # ------------------------------------------------------------------#
  # Detect method
  # ------------------------------------------------------------------#
  method <- dplyr::case_when(
    all(
      c(
        "label", "est", "cil_1sided",
        "ciu", "p", "p_adj"
      ) %in% names(vims)
    ) ~ "survML",
    
    all(c("variable", "est") %in% names(vims)) ~ "survex",
    
    TRUE ~ NA_character_
  )
  
  if (is.na(method)) {
    stop("Unrecognized VIM results format.")
  }
  
  # ------------------------------------------------------------------#
  # Checks
  # ------------------------------------------------------------------#
  required_vims <- c(
    "data_type",
    "landmark_time",
    "est"
  )
  
  if (compare) {
    required_vims <- c(required_vims, "cohort")
  }
  
  missing_vims <- setdiff(required_vims, names(vims))
  
  if (length(missing_vims) > 0) {
    stop(
      "Missing columns in `vims`: ",
      paste(missing_vims, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------------#
  # Grouping variables
  # ------------------------------------------------------------------#
  group_vars <- c(
    if (compare) "cohort",
    "data_type",
    "landmark_time"
    
  )
  
  # ------------------------------------------------------------------#
  # Formatting helpers
  # ------------------------------------------------------------------#
  .format_survml_table <- function(dat) {
    
    dat %>%
      gt::gt() %>%
      gt::cols_label(
        label = gt::md("**Pathway**"),
        est_ci = gt::md("**VIM [95% CI]**"),
        p = gt::md("***p*-value**"),
        p_adj = gt::md("***q*-value**")
      ) %>%
      gt::fmt(
        columns = c(p, p_adj),
        fns = gtsummary::style_pvalue
      ) %>%
      gt::tab_style_body(
        style = gt::cell_text(weight = "bold"),
        columns = c(p, p_adj),
        fn = \(x) !is.na(x) & x < 0.05
      ) %>%
      gt::tab_source_note(
        paste(
          c(
            "VIM: variable importance point estimate",
            "CI: confidence interval",
            "q-value: Benjamini–Hochberg-adjusted p-value"
          ),
          collapse = "; "
        )
      )
  }
  
  .format_survex_table <- function(dat) {
    
    dat %>%
      gt::gt() %>%
      gt::cols_label(
        label = gt::md("**Pathway**"),
        est = gt::md("**VIM**")
      ) %>%
      gt::fmt_number(
        columns = est,
        decimals = 3
      ) %>%
      gt::tab_source_note(
        "VIM: variable importance point estimate"
      )
  }
  
  # ------------------------------------------------------------------#
  # Prepare data
  # ------------------------------------------------------------------#
  if (method == "survML") {
    
    prepared_vims <- vims %>%
      dplyr::mutate(
        est_ci = sprintf(
          "%.3f [%.3f–%.3f]",
          est,
          cil_1sided,
          ciu
        )
      )
    
  } else {
    
    prepared_vims <- vims 
  }
  
  # ------------------------------------------------------------------#
  # Separate tables: compare = FALSE
  # ------------------------------------------------------------------#
  if (!compare) {
    
    tables <- group_split_custom(
      prepared_vims,
      !!!rlang::syms(group_vars)
    )
    
    if (method == "survML") {
      
      return(
        purrr::map(
          tables,
          ~ .x %>%
            dplyr::select(
              label,
              est_ci,
              p,
              p_adj
            ) %>%
            .format_survml_table()
        )
      )
    }
    
    return(
      purrr::map(
        tables,
        ~ .x %>%
          dplyr::select(
            label,
            est
          ) %>%
          .format_survex_table()
      )
    )
  }
  
  # ------------------------------------------------------------------#
  # Combinations to compare
  # ------------------------------------------------------------------#
  combinations <- prepared_vims %>%
    dplyr::distinct(
      data_type,
      landmark_time
    ) %>%
    dplyr::arrange(
      data_type,
      landmark_time
    )%>%
    dplyr::filter(!is.na(data_type))
  
  # ------------------------------------------------------------------#
  # Comparison tables: survML
  # ------------------------------------------------------------------#
  if (method == "survML") {
    
    tables <- purrr::pmap(
      combinations,
      function(data_type, landmark_time) {
        
        dat <- prepared_vims %>%
          dplyr::filter(
            data_type == !!data_type,
            landmark_time == !!landmark_time
          ) %>%
          dplyr::mutate(
            cohort = dplyr::recode(
              cohort,
              "Bio-RAIDs" = "raids",
              "TCGA-CESC" = "tcga"
            )
          ) %>%
          dplyr::select(
            cohort,
            label,
            est_ci,
            p,
            p_adj
          ) %>%
          tidyr::pivot_wider(
            id_cols = label,
            names_from = cohort,
            values_from = c(est_ci, p, p_adj),
            names_sep = "_"
          ) %>%
          dplyr::select(
            label,
            dplyr::any_of(
              c(
                "est_ci_raids",
                "p_raids",
                "p_adj_raids",
                "est_ci_tcga",
                "p_tcga",
                "p_adj_tcga"
              )
            )
          )
        
        p_columns <- intersect(
          c(
            "p_raids",
            "p_adj_raids",
            "p_tcga",
            "p_adj_tcga"
          ),
          names(dat)
        )
        
        tbl <- dat %>%
          gt::gt()
        
        if (any(endsWith(names(dat), "_raids"))) {
          tbl <- tbl %>%
            gt::tab_spanner(
              columns = gt::ends_with("_raids"),
              label = gt::md("**Bio-RAIDs**")
            )
        }
        
        if (any(endsWith(names(dat), "_tcga"))) {
          tbl <- tbl %>%
            gt::tab_spanner(
              columns = gt::ends_with("_tcga"),
              label = gt::md("**TCGA-CESC**")
            )
        }
        
        tbl %>%
          gt::cols_label(
            label = gt::md("**Pathway**"),
            est_ci_raids = gt::md("**VIM [95% CI]**"),
            p_raids = gt::md("***p*-value**"),
            p_adj_raids = gt::md("***q*-value**"),
            est_ci_tcga = gt::md("**VIM [95% CI]**"),
            p_tcga = gt::md("***p*-value**"),
            p_adj_tcga = gt::md("***q*-value**"),
            .fn = identity
          ) %>%
          gt::fmt(
            columns = dplyr::all_of(p_columns),
            fns = gtsummary::style_pvalue
          ) %>%
          gt::tab_style_body(
            style = gt::cell_text(weight = "bold"),
            columns = dplyr::all_of(p_columns),
            fn = \(x) !is.na(x) & x < 0.05
          ) %>%
          gt::tab_source_note(
            paste(
              c(
                "VIM: variable importance measure",
                "CI: confidence interval",
                "q-value: Benjamini–Hochberg-adjusted p-value"
              ),
              collapse = "; "
            )
          )
      }
    )
    
  } else {
    
    # ----------------------------------------------------------------#
    # Comparison tables: survex
    # ----------------------------------------------------------------#
    tables <- purrr::pmap(
      combinations,
      function(data_type, landmark_time) {
        
        dat <- prepared_vims %>%
          dplyr::filter(
            data_type == !!data_type,
            landmark_time == !!landmark_time
          ) %>%
          dplyr::mutate(
            cohort = dplyr::recode(
              cohort,
              "Bio-RAIDs" = "raids",
              "TCGA-CESC" = "tcga"
            )
          ) %>%
          dplyr::select(
            label,
            cohort,
            est
          ) %>%
          tidyr::pivot_wider(
            id_cols = label,
            names_from = cohort,
            values_from = est,
            names_prefix = "vimp_"
          )
        
        tbl <- dat %>%
          gt::gt()
        
        # if ("vimp_raids" %in% names(dat)) {
        #   tbl <- tbl %>%
        #     gt::tab_spanner(
        #       columns = vimp_raids,
        #       label = gt::md("**Bio-RAIDs**")
        #     )
        # }
        # 
        # if ("vimp_tcga" %in% names(dat)) {
        #   tbl <- tbl %>%
        #     gt::tab_spanner(
        #       columns = vimp_tcga,
        #       label = gt::md("**TCGA-CESC**")
        #     )
        # }
        
        tbl %>%
          gt::cols_label(
            label = gt::md("**Pathway**"),
            vimp_raids = gt::md("**Bio-RAIDs**"),
            vimp_tcga = gt::md("**TCGA-CESC**"),
            .fn = identity
          ) %>%
          gt::fmt_number(
            columns = gt::starts_with("vimp_"),
            decimals = 3
          ) %>%
          gt::tab_source_note(
            "VIM: variable importance measure"
          )
      }
    )
  }
  
  names(tables) <- paste(
    combinations$data_type,
    combinations$landmark_time,
    sep = "_"
  )
  
  purrr::compact(tables)
}

# Summary tables of Top-N ranked VIM estimates
tbl_top10_vimp <- function(
    vims,
    compare = TRUE,
    n_top = 10,
    digits = 3
) {
  
  # ------------------------------------------------------------------#
  # Detect method
  # ------------------------------------------------------------------#
  method <- dplyr::case_when(
    all(
      c(
        "label", "est", "cil_1sided",
        "ciu", "p", "p_adj"
      ) %in% names(vims)
    ) ~ "survML",
    
    all(c("variable", "est") %in% names(vims)) ~ "survex",
    
    TRUE ~ NA_character_
  )
  
  if (is.na(method)) {
    stop("Unrecognized VIM results format.")
  }
  
  # ------------------------------------------------------------------#
  # Checks
  # ------------------------------------------------------------------#
  required_vims <- c(
    "data_type",
    "landmark_time",
    "est"
  )
  
  if (compare) {
    required_vims <- c(required_vims, "cohort")
  }
  
  missing_vims <- setdiff(required_vims, names(vims))
  
  if (length(missing_vims) > 0) {
    stop(
      "Missing columns in `vims`: ",
      paste(missing_vims, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------------#
  # Select Top-N variables
  # ------------------------------------------------------------------#
  group_vars <- c(
    if (compare) "cohort",
    "data_type",
    "landmark_time"
  )
  
  top_vims <- vims %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_vars))
    ) %>%
    dplyr::arrange(
      dplyr::desc(est),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      rank = dplyr::row_number()
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------------#
  # Standardise label column
  # ------------------------------------------------------------------#
  if (method == "survML") {
    
    top_vims <- top_vims %>%
      dplyr::mutate(
        est_ci = sprintf(
          paste0("%.",digits,"f [%.",digits,"f–%.",digits,"f]"),
          # "%.3f [%.3f–%.3f]",
          est,
          cil_1sided,
          ciu
        )
      )
    
  } 
  
  # ------------------------------------------------------------------#
  # Shared pathways by landmark and data type
  # ------------------------------------------------------------------#
  if (compare) {
    
    common_features <- top_vims %>%
      dplyr::group_by(
        landmark_time,
        data_type,
        label
      ) %>%
      dplyr::summarise(
        n_cohorts = dplyr::n_distinct(cohort),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_cohorts > 1)
    
  } else {
    
    common_features <- NULL
  }
  
  # ------------------------------------------------------------------#
  # Formatting helpers
  # ------------------------------------------------------------------#
  .survml_notes <- function(add_common_note = FALSE) {
    
    notes <- c(
      "VIM: variable importance point estimate",
      "CI: confidence interval",
      "q-value: Benjamini–Hochberg-adjusted p-value"
    )
    
    if (add_common_note) {
      notes <- c(
        notes,
        paste0(
          "\u2020 Pathway ranked among the top ",
          n_top,
          " in both Bio-RAIDs and TCGA-CESC"
        )
      )
    }
    
    paste(notes, collapse = "; ")
  }
  
  .survex_notes <- function(add_common_note = FALSE) {
    
    notes <- "VIM: variable importance measure"
    
    if (add_common_note) {
      notes <- c(
        notes,
        paste0(
          "\u2020 Pathway ranked among the top ",
          n_top,
          " in both Bio-RAIDs and TCGA-CESC"
        )
      )
    }
    
    paste(notes, collapse = "; ")
  }
  
  .format_survml_table <- function(dat, add_common_note = FALSE) {
    
    dat %>%
      gt::gt() %>%
      gt::cols_label(
        rank = gt::md("**Rank**"),
        label = gt::md("**Pathway**"),
        est_ci = gt::md("**VIM [95% CI]**"),
        p = gt::md("***p*-value**"),
        p_adj = gt::md("***q*-value**")
      ) %>%
      gt::fmt(
        columns = c(p, p_adj),
        fns = gtsummary::style_pvalue
      ) %>%
      gt::tab_style_body(
        style = gt::cell_text(weight = "bold"),
        columns = c(p, p_adj),
        fn = \(x) !is.na(x) & x < 0.05
      ) %>%
      gt::tab_source_note(
        .survml_notes(add_common_note)
      )
  }
  
  .format_survex_table <- function(dat, add_common_note = FALSE) {
    
    dat %>%
      gt::gt() %>%
      gt::cols_label(
        rank = gt::md("**Rank**"),
        label = gt::md("**Pathway**"),
        est = gt::md("**VIM**")
      ) %>%
      gt::fmt_number(
        columns = est,
        decimals = digits
      ) %>%
      gt::tab_source_note(
        .survex_notes(add_common_note)
      )
  }
  
  # ------------------------------------------------------------------#
  # Separate tables: compare = FALSE
  # ------------------------------------------------------------------#
  if (!compare) {
    
    tables <- group_split_custom(
      top_vims,
      !!!rlang::syms(group_vars)
    )
    
    if (method == "survML") {
      
      return(
        purrr::map(
          tables,
          ~ .x %>%
            dplyr::select(
              rank,
              label,
              est_ci,
              p,
              p_adj
            ) %>%
            .format_survml_table()
        )
      )
    }
    
    return(
      purrr::map(
        tables,
        ~ .x %>%
          dplyr::select(
            rank,
            label,
            est
          ) %>%
          .format_survex_table()
      )
    )
  }
  
  # ------------------------------------------------------------------#
  # Landmark × data-type combinations
  # ------------------------------------------------------------------#
  combinations <- top_vims %>%
    dplyr::distinct(
      data_type,
      landmark_time
    ) %>%
    dplyr::arrange(
      data_type,
      landmark_time
    )%>%
    dplyr::filter(!is.na(data_type))
  
  # ------------------------------------------------------------------#
  # Build one comparison table
  # ------------------------------------------------------------------#
  .make_comparison_table <- function(
    data_type,
    landmark_time
  ) {
    
    current_common <- common_features %>%
      dplyr::filter(
        data_type == !!data_type,
        landmark_time == !!landmark_time
      ) %>%
      dplyr::pull(label)
    
    add_common_note <- grepl(
      "pathways",
      data_type,
      ignore.case = TRUE
    ) && length(current_common) > 0
    
    dat <- top_vims %>%
      dplyr::filter(
        data_type == !!data_type,
        landmark_time == !!landmark_time
      ) %>%
      dplyr::mutate(
        label = dplyr::if_else(
          add_common_note & label %in% current_common,
          paste0(label, "\u2020"),
          label
        ),
        cohort = dplyr::recode(
          cohort,
          "Bio-RAIDs" = "raids",
          "TCGA-CESC" = "tcga"
        )
      )
    
    if (method == "survML") {
      
      dat <- dat %>%
        dplyr::select(
          cohort,
          rank,
          label,
          est_ci,
          p,
          p_adj
        ) %>%
        tidyr::pivot_wider(
          id_cols = rank,
          names_from = cohort,
          values_from = c(label, est_ci, p, p_adj),
          names_glue = "{.value}_{cohort}"
        ) %>%
        dplyr::select(
          rank,
          dplyr::any_of(
            c(
              "label_raids",
              "est_ci_raids",
              "p_raids",
              "p_adj_raids",
              "label_tcga",
              "est_ci_tcga",
              "p_tcga",
              "p_adj_tcga"
            )
          )
        )
      
      p_columns <- intersect(
        c(
          "p_raids",
          "p_adj_raids",
          "p_tcga",
          "p_adj_tcga"
        ),
        names(dat)
      )
      
      tbl <- dat %>%
        gt::gt()
      
      if (any(endsWith(names(dat), "_raids"))) {
        tbl <- tbl %>%
          gt::tab_spanner(
            columns = gt::ends_with("_raids"),
            label = gt::md("**Bio-RAIDs**")
          )
      }
      
      if (any(endsWith(names(dat), "_tcga"))) {
        tbl <- tbl %>%
          gt::tab_spanner(
            columns = gt::ends_with("_tcga"),
            label = gt::md("**TCGA-CESC**")
          )
      }
      
      return(
        tbl %>%
          gt::cols_label(
            rank = gt::md("**Rank**"),
            label_raids = gt::md("**Pathway**"),
            est_ci_raids = gt::md("**VIM [95% CI]**"),
            p_raids = gt::md("***p*-value**"),
            p_adj_raids = gt::md("***q*-value**"),
            label_tcga = gt::md("**Pathway**"),
            est_ci_tcga = gt::md("**VIM [95% CI]**"),
            p_tcga = gt::md("***p*-value**"),
            p_adj_tcga = gt::md("***q*-value**"),
            .fn = identity
          ) %>%
          gt::fmt(
            columns = dplyr::all_of(p_columns),
            fns = gtsummary::style_pvalue
          ) %>%
          gt::tab_style_body(
            style = gt::cell_text(weight = "bold"),
            columns = dplyr::all_of(p_columns),
            fn = \(x) !is.na(x) & x < 0.05
          ) %>%
          gt::tab_source_note(
            .survml_notes(add_common_note)
          )
      )
    }
    
    dat <- dat %>%
      dplyr::select(
        cohort,
        rank,
        label,
        est
      ) %>%
      tidyr::pivot_wider(
        id_cols = rank,
        names_from = cohort,
        values_from = c(label, est),
        names_glue = "{.value}_{cohort}"
      ) %>%
      dplyr::select(
        rank,
        dplyr::any_of(
          c(
            "label_raids",
            "est_raids",
            "label_tcga",
            "est_tcga"
          )
        )
      )
    
    tbl <- dat %>%
      gt::gt()
    
    if (any(endsWith(names(dat), "_raids"))) {
      tbl <- tbl %>%
        gt::tab_spanner(
          columns = gt::ends_with("_raids"),
          label = gt::md("**Bio-RAIDs**")
        )
    }

    if (any(endsWith(names(dat), "_tcga"))) {
      tbl <- tbl %>%
        gt::tab_spanner(
          columns = gt::ends_with("_tcga"),
          label = gt::md("**TCGA-CESC**")
        )
    }
    
    tbl %>%
      gt::cols_label(
        rank = gt::md("**Rank**"),
        label_raids = gt::md("**Pathway**"),
        est_raids = gt::md("**VIM**"),
        label_tcga = gt::md("**Pathway**"),
        est_tcga = gt::md("**VIM**"),
        .fn = identity
      ) %>%
      gt::fmt_number(
        columns = gt::ends_with(c("_raids", "_tcga")) &
          gt::starts_with("est"),
        decimals = digits
      ) %>%
      gt::tab_source_note(
        .survex_notes(add_common_note)
      )
  }
  
  # ------------------------------------------------------------------#
  # Generate comparison tables
  # ------------------------------------------------------------------#
  tables <- purrr::pmap(
    combinations,
    .make_comparison_table
  )
  
  names(tables) <- paste(
    combinations$data_type,
    combinations$landmark_time,
    sep = "_"
  )
  
  purrr::compact(tables)
}


# Univariate Cox regression of top-10 ranked pathways 
tbl_cox_univ_top10_vimp <- function(
    vims,
    data,
    compare = TRUE,
    n_top = 10
) {
  
  # ------------------------------------------------------------------#
  # Checks
  # ------------------------------------------------------------------#
  required_vims <- c(
    "data_type",
    "landmark_time",
    "variable_name",
    "est"
  )
  
  if (compare) {
    required_vims <- c(required_vims, "cohort")
  }
  
  missing_vims <- setdiff(required_vims, names(vims))
  
  if (length(missing_vims) > 0) {
    stop(
      "Missing columns in `vims`: ",
      paste(missing_vims, collapse = ", ")
    )
  }
  
  required_data <- c("time", "event")
  
  if (compare) {
    required_data <- c(required_data, "cohort")
  }
  
  missing_data <- setdiff(required_data, names(data))
  
  if (length(missing_data) > 0) {
    stop(
      "Missing columns in `data`: ",
      paste(missing_data, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------------#
  # Select Top-N pathways
  # ------------------------------------------------------------------#
  group_vars <- c(
    if (compare) "cohort",
    "data_type",
    "landmark_time"
  )
  
  top_vims <- vims %>%
    dplyr::filter(
      grepl("pathways", data_type, ignore.case = TRUE)
    ) %>%
    dplyr::group_by(
      dplyr::across(dplyr::all_of(group_vars))
    ) %>%
    dplyr::arrange(
      dplyr::desc(est),
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      rank = dplyr::row_number()
    ) %>%
    dplyr::ungroup()
  
  # ------------------------------------------------------------------#
  # Pathways shared between cohorts, by landmark and data type
  # ------------------------------------------------------------------#
  if (compare) {
    
    common_pathways <- top_vims %>%
      dplyr::group_by(
        landmark_time,
        data_type,
        variable_name
      ) %>%
      dplyr::summarise(
        n_cohorts = dplyr::n_distinct(cohort),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_cohorts > 1)
    
    notes <- c(
      "HR: hazard ratio",
      "CI: confidence interval",
      "q-value: Benjamini–Hochberg-adjusted p-value",
      paste0(
        "\u2020 Pathway ranked among the top ",
        n_top,
        " in both Bio-RAIDs and TCGA-CESC"
      )
    )
    
  } else {
    
    common_pathways <- NULL
    
    notes <- c(
      "HR: hazard ratio",
      "CI: confidence interval",
      "q-value: Benjamini–Hochberg-adjusted p-value"
    )
  }
  
  # ------------------------------------------------------------------#
  # Formatting helper
  # ------------------------------------------------------------------#
  .format_table <- function(tbl) {
    
    tbl %>%
      gtsummary::remove_source_note() %>%
      gtsummary::remove_abbreviation() %>%
      gtsummary::remove_footnote_header() %>%
      gtsummary::modify_source_note(
        paste(notes, collapse = "; ")
      ) %>%
      gtsummary::modify_column_alignment(
        columns = dplyr::everything(),
        align = "left"
      )
  }
  
  # ------------------------------------------------------------------#
  # Fit one set of univariable Cox models
  # ------------------------------------------------------------------#
  .fit_one_table <- function(x) {
    
    vars <- x %>%
      dplyr::pull(variable_name) %>%
      unique()
    
    landmark_value <- unique(x$landmark_time)
    data_type_value <- unique(x$data_type)
    
    analysis_data <- data
    
    if (compare) {
      
      cohort_value <- unique(x$cohort)
      
      analysis_data <- analysis_data %>%
        dplyr::filter(cohort == cohort_value)
    }
    
    # Labels specific to the current landmark and data type
    current_labels <- labelled::var_label(analysis_data)
    
    if (compare) {
      
      current_common <- common_pathways %>%
        dplyr::filter(
          landmark_time == landmark_value,
          data_type == data_type_value
        ) %>%
        dplyr::pull(variable_name)
      
      current_labels <- purrr::imap(
        current_labels,
        ~ if (.y %in% current_common) {
          paste0(.x, "\u2020")
        } else {
          .x
        }
      )
    }
    
    analysis_data <- analysis_data %>%
      labelled::set_variable_labels(
        .labels = current_labels
      )
    
    analysis_data %>%
      gtsummary::tbl_uvregression(
        method = survival::coxph,
        y = survival::Surv(time, event),
        include = dplyr::all_of(vars),
        exponentiate = TRUE,
        hide_n = TRUE
      ) %>%
      gtsummary::add_global_p() %>%
      gtsummary::bold_p() %>%
      gtsummary::add_q(method = "BH") %>%
      gtsummary::bold_p(q = TRUE) %>%
      gtsummary::modify_column_merge(
        pattern = "{estimate} [{conf.low}–{conf.high}]",
        rows = !is.na(estimate)
      ) %>%
      gtsummary::modify_table_body(
        ~ .x %>%
          dplyr::mutate(
            rank = dplyr::row_number(),
            .before = label
          )
      ) %>%
      gtsummary::modify_header(
        rank = "**Rank**",
        label = "**Pathway**",
        estimate = "**HR [95% CI]**"
      )
  }
  
  # ------------------------------------------------------------------#
  # Split and fit
  # ------------------------------------------------------------------#
  tables <- group_split_custom(
    top_vims,
    !!!rlang::syms(group_vars)
  ) %>%
    purrr::map(.fit_one_table)
  
  # ------------------------------------------------------------------#
  # Return separate tables
  # ------------------------------------------------------------------#
  if (!compare) {
    
    return(
      purrr::map(
        tables,
        .format_table
      )
    )
  }
  
  # ------------------------------------------------------------------#
  # Merge Bio-RAIDs and TCGA-CESC tables
  # ------------------------------------------------------------------#
  combinations <- top_vims %>%
    dplyr::distinct(
      data_type,
      landmark_time
    ) %>%
    dplyr::arrange(
      data_type,
      landmark_time
    )
  
  merged_tables <- purrr::pmap(
    combinations,
    function(landmark_time, data_type) {
      
      table_names <- paste(
        c("Bio-RAIDs", "TCGA-CESC"),
        data_type,
        landmark_time,
        sep = "_"
      )
      
      missing_tables <- setdiff(
        table_names,
        names(tables)
      )
      
      if (length(missing_tables) > 0) {
        warning(
          "Missing table(s): ",
          paste(missing_tables, collapse = ", ")
        )
      }
      
      selected_tables <- tables[
        intersect(table_names, names(tables))
      ]
      
      if (length(selected_tables) == 0) {
        return(NULL)
      }
      
      spanners <- names(selected_tables) %>%
        stringr::str_extract(
          "Bio-RAIDs|TCGA-CESC"
        ) %>%
        paste0("**", ., "**")
      
      gtsummary::tbl_merge(
        tbls = selected_tables,
        tab_spanner = spanners,
        merge_vars = "rank"
      ) %>%
        # gtsummary::modify_caption(
        #   paste0(
        #     "**",
        #     landmark_time,
        #     "-month horizon — ",
        #     data_type,
        #     "**"
        #   )
        # ) %>%
        .format_table()
    }
  )
  
  names(merged_tables) <- paste(
    combinations$data_type,
    combinations$landmark_time,
    sep = "_"
  )
  
  purrr::compact(merged_tables)
}

# Overlap of top-10 ranked pathways 
tbl_top_overlap_vimp <- function(vims, k = 10, B = 10000, seed = 123){
  
  .top_overlap <- function(vims_raids, vims_tcga, k, B, seed, modality){
    
    set.seed(123)
    
    # Spearman correlation
    spearman_rho <- cor(vims_raids, vims_tcga, method = "spearman")
    spearman_p <- cor.test(vims_raids, vims_tcga,
                           method = "spearman", exact = FALSE)$p.value
    
    universe <- intersect(names(vims_raids), names(vims_tcga))
    
    # Top-k
    ranking_raids <- head(sort(vims_raids, decreasing = T), k)
    ranking_tcga <- head(sort(vims_tcga, decreasing = T), k)
    
    # Number of overlapping features 
    overlap <- intersect(names(ranking_raids), names(ranking_tcga))
    n_overlap <- length(overlap)
    
    # Percentage overlap 
    pct_overlap <- sprintf("%.0f", 100 * n_overlap / k)
    
    # Overlap features
    overlap_features <- paste(overlap, collapse = ",\n")
    
    # Jaccard index
    jaccard <- sprintf("%.2f", n_overlap /
                         length(union(names(ranking_raids), 
                                      names(ranking_tcga))))
    
    # Permutation test to assess whether the observed Top-10 overlap 
    # was greater than expected by chance             
    null_overlap <- replicate(
      B,
      {
        random_top <- sample(universe, k)
        
        length(intersect(names(ranking_raids), random_top))
      }
    )
    
    p_perm <- (sum(null_overlap >= n_overlap) + 1) / (B + 1)
    
    metric <- c(
      "Spearman's ρ",
      "Spearman's p-value", 
      paste0("Top-", k , " overlap (%)"), 
      paste0("Top-", k , " overlapping features"), 
      "Jaccard index", 
      paste0("Permutation p-value (top-",k, " overlap)")
    )

    tibble(spearman_rho = sprintf("%.3f", spearman_rho), 
           spearman_p = gtsummary::style_pvalue(spearman_p), 
           # n_overlap = sprintf("%.0f", n_overlap), 
           pct_overlap = pct_overlap,
           overlap_features = overlap_features,
           jaccard = jaccard, 
           p_perm = gtsummary::style_pvalue(p_perm))%>%
      pivot_longer(everything(), names_to = "Metric", values_to = modality)%>%
      dplyr::mutate(Metric = metric)
  }
  
  vims_split <- vims %>%
    dplyr::filter(!grepl("processes", data_type))%>%
    group_split_custom(cohort, data_type, landmark_time)%>%
    purrr::map(
      ~.x %>% 
        dplyr::select(label, est)%>%
        tibble::deframe()
    )%>%
    purrr::list_flatten()
  
  dna_top_k <- .top_overlap(vims_raids = vims_split$`Bio-RAIDs_dna_pathways_24`, 
                            vims_tcga  = vims_split$`TCGA-CESC_dna_pathways_24`,
                            k = k, B = B, seed = seed, modality = "DNA")
  rna_top_k <- .top_overlap(vims_raids = vims_split$`Bio-RAIDs_rna_pathways_24`, 
                            vims_tcga  = vims_split$`TCGA-CESC_rna_pathways_24`,
                            k = k, B = B, seed = seed, modality = "RNA")
  overall_top_k <- .top_overlap(vims_raids = c(vims_split$`Bio-RAIDs_dna_pathways_24`,
                                               vims_split$`Bio-RAIDs_rna_pathways_24`), 
                                vims_tcga  = c(vims_split$`TCGA-CESC_dna_pathways_24`,
                                               vims_split$`TCGA-CESC_rna_pathways_24`),
                                k = k, B = B, seed = seed, modality = "Overall")
  
  list(dna_top_k, rna_top_k, overall_top_k)%>%
    purrr::reduce(full_join, by = "Metric")%>%
    gt::gt()%>%
    gt::tab_style(style = list(cell_text(weight = "bold")),
                  locations = cells_body(columns = 1))%>%
    gt::tab_style(style = list(cell_text(weight = "bold")),
                  locations = cells_title())%>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels(everything())
    )
}

# Compare VIM point estimates with Spearman correlation 
scatter_plot_vimp <- function(vims,
                              annotate_x = 0.11,
                              annotate_overall_y = 0.25,
                              annotate_data_type_y = c(0.24,0.23)) {
  
  dat <- vims %>%
    dplyr::select(variable, data_type, cohort, est) %>%
    dplyr::filter(data_type %in% c("dna_pathways","rna_pathways"))%>%
    dplyr::mutate(data_type = toupper(gsub("_pathways","", data_type)))%>%
    tidyr::pivot_wider(
      names_from  = cohort,
      values_from = est,
      names_prefix = "est_"
    ) %>%
    dplyr::mutate(data_type = gsub("_", " ", data_type)) %>%
    dplyr::filter(!is.na(`est_Bio-RAIDs`) & !is.na(`est_TCGA-CESC`))
  
  # --- Spearman overall
  rho_all <- suppressWarnings(cor(dat$`est_Bio-RAIDs`, dat$`est_TCGA-CESC`, 
                                  method = "spearman"))
  p_all   <- suppressWarnings(cor.test(dat$`est_Bio-RAIDs`, dat$`est_TCGA-CESC`, 
                                       method = "spearman", exact = FALSE)$p.value)
  n_all   <- nrow(dat)
  
  lab_all <- paste0(
    "Overall: \u03c1=", sprintf("%.3f", rho_all),
    ", ", gtsummary::style_pvalue(p_all, prepend_p = T),
    ", n=", n_all
  )
  
  # --- Spearman par data_type + positions d'annotation
  ann <- dat %>%
    dplyr::group_by(data_type) %>%
    dplyr::summarise(
      n = dplyr::n(),
      rho = cor(`est_Bio-RAIDs`, `est_TCGA-CESC`, 
                method = "spearman"),
      p   = cor.test(`est_Bio-RAIDs`, `est_TCGA-CESC`, 
                     method = "spearman", exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    dplyr::arrange(data_type) %>%                 # ordre stable
    dplyr::mutate(
      label = paste0(
        data_type, ": ",
        "\u03c1=", sprintf("%.3f", rho),
        ", ", gtsummary::style_pvalue(p, prepend_p = T),
        ", n=", n
      ),
      y = annotate_data_type_y, 
    )
  
  
  ggplot2::ggplot(dat, ggplot2::aes(x = `est_Bio-RAIDs`, y = `est_TCGA-CESC`, 
                                    color = data_type)) +
    ggplot2::geom_point(alpha = 0.8) +
    geomtextpath::geom_labelsmooth(
      ggplot2::aes(label = data_type),
      fill = "white",
      method = "lm", formula = y ~ x,
      size = 3, linewidth = 1, boxlinewidth = 0.4
    ) +
    ggplot2::geom_smooth(
      data = dat,
      aes(x = `est_Bio-RAIDs`, y = `est_TCGA-CESC`),
      inherit.aes = FALSE,
      method = "lm", formula = y ~ x,
      color = "black", linewidth = 1.1, linetype = 2
    ) +
    
    # Annotation overall 
    ggplot2::annotate("text", 
                      x = annotate_x, 
                      y = annotate_overall_y, 
                      label = lab_all,
                      hjust = -0.05, #vjust = 1.2, 
                      size = 3.2) +
    # Annotations by data type 
    ggplot2::geom_text(
      data = ann,
      aes(x = annotate_x, y = y, label = label, color = data_type),
      inherit.aes = FALSE,
      hjust = -0.05, #vjust = ann$y,
      size = 3,
      show.legend = FALSE
    )  +
    ggplot2::labs(
      x = "VIM (Bio-RAIDs)",
      y = "VIM (TCGA-CESC)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = "none")
}

