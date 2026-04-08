
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
library(patchwork)
library(survival)  
library(survML)  
library(SuperLearner)

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
  
  group_pathway_features <- get_pathways_process() %>%
    dplyr::filter(data_type == "RNA_pathways") %>%
    dplyr::transmute(process, variable = tolower(variable)) %>%
    dplyr::group_by(process) %>%
    dplyr::summarise(variable = list(variable), .groups = "drop") %>%
    tibble::deframe()%>%
    purrr::imap_chr(~paste(which(names(X) %in% .x), collapse = ",")
  )
  feature_groups <- c(pathway_features, group_pathway_features)
  
  list(time = data$time,
       event = data$event,
       X = X, 
       base_features = base_features, 
       feature_groups = feature_groups)
}

# Plotting and table helpers ----------------------------------------------

#' Get the variable importance estimates.
#'@param compute_vimp The object returned by compute_vimp_*()

get_vimp_survML_est <- function(compute_vimp){
  do.call(rbind, lapply(compute_vimp, function(x) x$result))%>%
    {rownames(.) <- NULL; .}%>%
    left_join(get_pathways_process(), by = "variable")%>%
    mutate(data_type = ifelse(is.na(data_type), "RNA_processes", data_type))
}

## Plot VIM estimates
plot_vimp_survML_est <- function(vims,
                                 xlab = "VIM ± 95% CI",
                                 ylab = "Variable",
                                 type = c("barplot", "dotplot"),
                                 fill_by = NULL,         
                                 fill_values = NULL,      
                                 fill_label = NULL,
                                 p_text_color = "#F8766D",
                                 vline0 = TRUE,
                                 vline_color = "red",
                                 vline_linetype = "dotted",
                                 vline_linewidth = 0.7,
                                 limits = c(0, 0.4),
                                 legend.position = "top",
                                 legend.direction = "horizontal") {
  
  type <- match.arg(type)
  
  df <- vims %>%
    # dplyr::arrange(dplyr::desc(est)) %>%
    dplyr::mutate(
      # ord_group = forcats::fct_reorder(variable, est),
      ord_group = as.factor(variable),
      est_lab = paste0(
        format(round(est, 2), nsmall = 2), " [",
        format(round(cil, 2), nsmall = 2), "-",
        format(round(ciu, 2), nsmall = 2), "] ",
        gtools::stars.pval(p)
      ),
      p_val = gtools::stars.pval(p)
    )
  
  df_p1 <- df %>% filter(process %in% c("Cellular component",
                                        "DNA damage", "Pathway"))
  df_p2 <- df %>% filter(process %in% c("Proliferation","Development"))
  df_p3 <- df %>% filter(process %in% c("Metabolic","Immune"))
  df_p4 <- df %>% filter(process %in% c("Signaling"))
  
  # Default fill handling
  if (is.null(fill_by)) {
    aes_fill <- NULL
  } else {
    # allow unquoted or character; normalize to string
    fill_by <- rlang::as_name(rlang::ensym(fill_by))
    aes_fill <- fill_by
    
    if (is.null(fill_label)) {
      fill_label <- dplyr::case_when(
        fill_by == "cohort" ~ "Cohort",
        # fill_by == "figo"   ~ "FIGO stage",
        TRUE                ~ fill_by
      )
    }
    
    if (is.null(fill_values)) {
      fill_values <- dplyr::case_when(
        fill_by == "cohort" ~ TRUE,
        # fill_by == "figo"   ~ TRUE,
        TRUE                ~ FALSE
      )
      fill_values <- if (isTRUE(fill_values)) {
        if (fill_by == "cohort") {
          c("TCGA-CESC" = "#377eb8", 
            "Bio-RAIDs" = "#1b9e77")
        } 
      } else {
        NULL
      }
    }
  }
  
  base_theme <- function(plot, legend.position = legend.position, 
                         legend.direction = legend.direction) {
    plot + 
      theme_minimal(base_size = 18) +
      theme(
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = legend.position,
        legend.direction = legend.direction, 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        strip.text = element_text(face = "bold", size = 18),
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey90", colour = NA),
        plot.margin = margin(8, 8, 8, 8),
        panel.spacing = unit(0.8, "lines"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA, 
                                    linewidth = 0.4)
      )
  } 
  
  .dotplot_vimp <- function(df, xlab, ylab, aes_fill, fill_values, 
                            p_text_color, vline0, 
                            facet = F, ncol = NULL, nrow = NULL, 
                            legend.position, legend.direction){
    
    plot <- ggplot(df, aes(x = .data[["est"]], y = .data[["ord_group"]], 
                   colour = if (!is.null(aes_fill)) .data[[aes_fill]]
                   else NULL)) +
      geom_point(position = position_dodge(width = 0.8, preserve = "single"),
                 size = 3) +
      geom_errorbarh(aes(xmin = .data[["cil"]], xmax = .data[["ciu"]]),
                     position = position_dodge(width = 0.8, preserve = "single"),
                     linewidth = 0.9) +
      scale_x_continuous(limits = limits)+
      scale_y_discrete(limits = rev)+
      labs(x = xlab, y = ylab, colour = fill_label) +
      geom_text(aes(x = .data[["ciu"]], label = .data[["p_val"]], 
                    hjust = -0.2),
                size = 6, col = p_text_color)
    
    if (isTRUE(vline0)) {
      plot <- plot + geom_vline(
        xintercept = 0,
        linetype = vline_linetype,
        color = vline_color,
        linewidth = vline_linewidth
      )
    }
    
    if (!is.null(aes_fill) && !is.null(fill_values)) {
      plot <- plot + scale_colour_manual(values = fill_values)
    }
    
    if (isTRUE(facet)){
      plot <- plot + 
        facet_wrap(~.data[["process"]],
                   nrow = nrow, ncol = ncol, scales = "free_y")
    }
    
    base_theme(plot, legend.position, legend.direction)
  }
  
  .barplot_vimp <- function(df, xlab, ylab, aes_fill, fill_values, 
                            p_text_color, vline0, 
                            facet = F, ncol = NULL, nrow = NULL, 
                            legend.position, legend.direction
                            ){
    
    plot <- ggplot(df, aes(x = .data[["est"]], 
                           y = .data[["ord_group"]], 
                           fill = if (!is.null(aes_fill)) 
                             .data[[aes_fill]] else NULL)) +
      geom_col(
        width = 0.55,
        position = position_dodge(width = 0.8, preserve = "single"),
        colour = NA,
        alpha = 0.95
      ) +
      geom_errorbarh(
        aes(xmin = .data[["cil"]], xmax = .data[["ciu"]]),
        position = position_dodge(width = 0.8, preserve = "single"),
        linewidth = 0.4,
        colour = "black"
      )+
      scale_x_continuous(limits = limits)+
      scale_y_discrete(limits = rev)+
      labs(x = xlab, y = ylab, fill = fill_label) +
      geom_text(aes(x = .data[["ciu"]], label = .data[["p_val"]], 
                    hjust = -0.2),
                size = 6, col = p_text_color)
    
    if (isTRUE(vline0)) {
      plot <- plot + geom_vline(
        xintercept = 0,
        linetype = vline_linetype,
        color = vline_color,
        linewidth = vline_linewidth
      )
    }
    if (!is.null(aes_fill) && !is.null(fill_values)) {
      plot <- plot + scale_fill_manual(values = fill_values)
    }
    
    if (isTRUE(facet)){
      plot <- plot + 
        facet_wrap(~.data[["process"]],
                   nrow = nrow, ncol = ncol, scales = "free_y")
    }
  
    base_theme(plot, legend.position, legend.direction)
  }

  if (type == "dotplot" & unique(vims$data_type)!="RNA_pathways") {
    
    plot <- .dotplot_vimp(df, xlab, ylab, aes_fill, fill_values,
                          p_text_color, vline0, 
                          legend.position = legend.position, 
                          legend.direction = legend.direction)
    
  } else if (type == "dotplot" & unique(vims$data_type)=="RNA_pathways")  {
    
    p1 <- .dotplot_vimp(df_p1, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 3, ncol = 1,
                        legend.position = legend.position, 
                        legend.direction = legend.direction)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p2 <- .dotplot_vimp(df_p2, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 2, ncol = 1,
                        legend.position = "none", 
                        legend.direction = legend.direction)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p3 <- .dotplot_vimp(df_p3, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 2, ncol = 1,
                        legend.position = "none", 
                        legend.direction = legend.direction)
    p4 <- .dotplot_vimp(df_p4, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 1, ncol = 1,
                        legend.position = "none", 
                        legend.direction = legend.direction)
    
    plot <- egg::ggarrange(p1, p2, p3, p4,
                           ncol = 2, nrow = 2,
                           widths = c(1, 1), heights = c(1, 1))
    
  } else if (type == "barplot"  & unique(vims$data_type)!="RNA_pathways") { 
    
    plot <- .barplot_vimp(df, xlab, ylab, aes_fill, fill_values, 
                          p_text_color, vline0, 
                          facet = F, nrow = NULL, ncol = NULL, 
                          legend.position, legend.direction) 

  } else if (type == "barplot"  & unique(vims$data_type)=="RNA_pathways"){
  
    p1 <- .barplot_vimp(df_p1, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 3, ncol = 1,
                        legend.position = legend.position, legend.direction)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p2 <- .barplot_vimp(df_p2, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 2, ncol = 1,
                        legend.position = "none", legend.direction)+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
    p3 <- .barplot_vimp(df_p3, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 2, ncol = 1,
                        legend.position = "none", legend.direction)
    p4 <- .barplot_vimp(df_p4, xlab, ylab, aes_fill, fill_values,
                        p_text_color, vline0,
                        facet = T, nrow = 1, ncol = 1,
                        legend.position = "none", legend.direction)
    
    plot <- egg::ggarrange(p1, p2, p3, p4,
                           ncol = 2, nrow = 2,
                           widths = c(1, 1), heights = c(1, 1))
  }
  
  plot
  
}

## Summary table of VIM estimates.
#' @param vims Tibble returned by compute_vimp_*()

tbl_vimp_survML <- function(vims){
  
  vims %>%
    dplyr::mutate(
      vimp = paste0(format(round(est, 3), nsmall = 3)," [",
                    format(round(cil_1sided, 3), nsmall = 3), "-",
                    format(round(ciu, 3), nsmall = 3), "]"),
      across(c(large_predictiveness, small_predictiveness), 
             ~format(round(.x , 2), nsmall = 2)))%>%
    dplyr::select(variable, 
                  large_predictiveness, small_predictiveness, 
                  vimp, p)%>%
    gt::gt()%>%
    gt::cols_label(
      variable = "", 
      large_predictiveness = md(paste("**V_full**")), 
      small_predictiveness = md(paste("**V_reduced**")), 
      vimp = md("**VIM [95% CI]**"),
      p = md("**p-value**"))%>%
    gt::tab_style_body(style = cell_text(weight = "bold"),
                   columns = p,
                   fn = function(x) x < 0.05)%>%
    gt::fmt(columns = p, 
            fns = function(x) ifelse(
              x < 0.05, formatC(x, format = "f", digits = 3),
              ifelse(x < 0.001, "<0.001",
                     formatC(x, format = "f", digits = 2))))%>%
    gt::tab_style(style = cell_text(weight = "bold"),
                  locations = cells_body(columns = variable))
}


# Compare VIM point estimates

scatter_plot_vimp <- function(df){
  
  dat <- df %>%
    dplyr::select(variable, data_type, cohort, est) %>%
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
                                       method = "spearman")$p.value)
  n_all   <- nrow(dat)
  
  lab_all <- paste0(
    "Overall: \u03c1=", sprintf("%.2f", rho_all),
    ", p=", ifelse(p_all < 0.001, "<0.001", formatC(p_all, format = "f", digits = 3)),
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
                     method = "spearman")$p.value,
      .groups = "drop"
    ) %>%
    dplyr::arrange(data_type) %>%                 # ordre stable
    dplyr::mutate(
      label = paste0(
        data_type, ": ",
        "\u03c1=", sprintf("%.2f", rho),
        ", p", ifelse(p < 0.001, "<0.001",
                      paste0("=", formatC(p, format = "f", digits = 3))),
        ", n=", n
      ),
      y = c(0.24,0.23,0.22), 
    )
  
  
  ggplot2::ggplot(dat, ggplot2::aes(x = `est_Bio-RAIDs`, y = `est_TCGA-CESC`, 
                                    color = data_type)) +
    ggplot2::geom_point(alpha = 0.8) +
    geom_labelsmooth(
      ggplot2::aes(label = data_type),
      fill = "white",
      method = "lm", formula = y ~ x,
      size = 3, linewidth = 1, boxlinewidth = 0.4
    ) +
    geom_smooth(
      data = dat,
      aes(x = `est_Bio-RAIDs`, y = `est_TCGA-CESC`),
      inherit.aes = FALSE,
      method = "lm", formula = y ~ x,
      color = "black", linewidth = 1.1, linetype = 2
    ) +
    
    # Annotation overall 
    ggplot2::annotate("text", x = 0.11, y = 0.25, label = lab_all,
                      hjust = -0.05, #vjust = 1.2, 
                      size = 3.2) +
    # Annotations by data type 
    geom_text(
      data = ann,
      aes(x = 0.11, y = y, label = label, color = data_type),
      inherit.aes = FALSE,
      hjust = -0.05, #vjust = ann$y,
      size = 3,
      show.legend = FALSE
    )  +
    ggplot2::labs(
      x = "VIM point estimates (Bio-RAIDs)",
      y = "VIM point estimates (TCGA-CESC)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::guides(color = "none")
}


# Get top-10 ranked features
get_top10_vimp_survML <- function(vims){
  vims %>%
    dplyr::arrange(desc(est))%>%
    dplyr::filter(row_number() %in% 1:10)%>%
    mutate(n = row_number(),
           est_ci = paste0(format(round(est, 3), nsmall = 3)," [",
                           format(round(cil_1sided, 3), nsmall = 3), "-",
                           format(round(ciu, 3), nsmall = 3), "]"))%>%
    dplyr::select(n, variable, est_ci, p)%>%
    gt::gt()%>%
    gt::cols_label(
      n =  md("**N°**"), 
      variable =  md("**Pathway**"), 
      est_ci =  md("**VIM [95% CI]**"),
      p = md("**p-value**"))%>%
    tab_style_body(style = cell_text(weight = "bold"),
                   columns = p,
                   fn = function(x) x < 0.05)%>%
    fmt(columns = p, 
        fns = function(x) ifelse(
          x < 0.05, formatC(x, format = "f", digits = 3),
          ifelse(x < 0.001, "<0.001",
                 formatC(x, format = "f", digits = 2))))%>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(columns = variable))
}


#' Overlap of top-10 ranked pathways 
get_overlap_top10_vimp_survML <- function(vims){
  top10 <- vims %>%
    map(~.x %>% 
          dplyr::filter(!grepl("RNA_processes", data_type))%>%
          group_split_custom(data_type, landmark_time)%>%
          map(~.x %>% dplyr::arrange(desc(est))%>%
                dplyr::filter(row_number() %in% 1:10)%>%
                dplyr::select(variable)%>%
                pull())%>%
          set_names(c("DNA-based pathways","RNA-based pathways")))
  
  map(names(top10[[1]]), function(modality) {
    
    raids <- top10$`Bio-RAIDs`[[modality]][1:10]
    tcga <- top10$`TCGA-CESC`[[modality]][1:10]
    
    tibble(
      modality = modality,
      n_overlap = length(intersect(raids, tcga)),
      pct_overlap = 100 * n_overlap / 10,
      overlap_features = paste(intersect(raids, tcga), collapse = ",\n")
    )
  }) %>%
    bind_rows()%>%
    gt::gt()%>%
    gt::cols_label(
      modality = md("**Data type**"),
      n_overlap = md("**Number of overlapping features (top 10)**"),
      pct_overlap = md("**Overlap (%)**"),
      overlap_features = md("**Overlapping features**"))
}

