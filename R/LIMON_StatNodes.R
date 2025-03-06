#' Individualized Microbial Node Variable selection
#'
#' This function takes a LIMON output object and pulls out all the significant interactions by your dependent variable of interest. First the model will fit a linear regression with a user specified dependent variable, for every unique interaction between taxa. Then it will filter to the results only significant based on the users p-value threshold. These results will be saved to "Significant_Interactions".
#'
#' @param LIMON_obj A LIMON object after the LIMON_IndEdges() step
#' @param dependent User specified dependent variable. For Linear Regression, must be a continuous variable. For Logistic Regression, it must be binary (0,1). For Multinomial, it must have at least three categories.
#' @param node_type Option to specify what node characteristic should be used for variable selection. Default/only currently supported option is "Edge_weight"
#' @param pval Option to change the p value from the linear models to filter the models by. Default is no threshold applied (pvalue of 1).
#' @param adjpval Option to change the BH adjusted pvalue filtering threshold. Default is no threshold applied (adjusted pvalue of 1).
#' @param time Option to specify what the time column is labeled as. default is "Time"
#' @param method Default is "lm" for linear regression. Other options are "logit" for logistic regression or "multinom" for multinomial
#' @param timepoints Option to provide a vector of which timepoints to run the analysis on instead of all timepoints
#' @param plot_nodes Option to plot the significant interactions for linear or logistic models. Default FALSE
#' @param estimate Option of an absolute value to filter the interaction edge weight coefficient estimate by
#' @param custom_colors Option to specify either an R color palette or supply a vector of colors to use in the plot_nodes. Default is NULL and uses Set3 from RColorBrewer.
#' @param point_size Dot size for the plot option. Default 7
#' @return A LIMON object with "Significant Interactions" which describes each model that passed thresholding, and "Nodes_Data" which will return a data table containing the metdata and only those taxa that were in the Significant_Interactions Table.
#'
#' @export
LIMON_StatNodes <- function (LIMON_obj, dependent, node_type = "Edge_weight", pval = 1,
                              adjpval = 1, time = "Time", method = "lm", timepoints = NULL,
                              plot_nodes = FALSE, estimate = NULL, custom_colors = NULL,
                              point_size = 7)
{
  library(broom)
  library(dplyr)
  library(nnet)
  library(tidyr)
  library(ggplot2)

  # Check messages to ensure the data is present
  if (is.null(LIMON_obj[["Merged_Edge_Table"]])) {
    print("Merged_Edge_Table is not present in the LIMON object. Must run LIMON_IndEdges() first.")
    return(LIMON_obj)
  }


  # Error Message if they select a linear model but dont have a continuous value with more than 3 values
  if (method == "lm") {
    # Warning if dependent is not binary
    unique_vals <- unique(LIMON_obj[["Merged_Edge_Table"]][[dependent]])
    if (length(unique_vals) <= 3 || !is.numeric(unique_vals)) {
      stop(paste("Error: The dependent variable", dependent,
                 "must be numeric and continuous (not categorical data) for a linear model. Consider using a different model for this variable"))}}



  # Error Message if they select logit and don't have a binary
  if (method == "logit") {
    # Warning if dependent is not binary
    unique_vals <- unique(LIMON_obj[["Merged_Edge_Table"]][[dependent]])
    if (length(unique_vals) != 2 || !all(unique_vals %in% c(0, 1))) {
      stop(paste("Error: The dependent variable", dependent,
                 "must be binary (0 and 1) for a logistic regression. Consider converting it to a binary variable or use the multinomial model if more than 2 categories."))}}


  # Error Message if they select multinom and don't have more than three categories
  if (method == "multinom") {
    # Warning if dependent is not multinomial
    unique_vals <- unique(LIMON_obj[["Merged_Edge_Table"]][[dependent]])
    if (length(unique_vals) <= 2) {
      stop(paste("Error: The dependent variable", dependent,
                 "must have 3 or more levels overall for a multinomial model. Consider using a different model"))}}

  # Option to filter down to timepoints of interest
  if (!is.null(timepoints)) {
    LIMON_obj[["Merged_Edge_Table"]] <- LIMON_obj[["Merged_Edge_Table"]] %>%
      filter(.data[[time]] == timepoints)
  }

  # Prep the Data First
  if (node_type == "Edge_weight") {
    # get the data
    Edge_Table <- LIMON_obj[["Merged_Edge_Table"]]

    # list to store all of the model data in
    node_results <- list()

    # filter to each time level
    for (time_level in unique(Edge_Table[[time]])) {
      edge_data_time_full <- Edge_Table %>% filter(Edge_Table[[time]] == time_level)

      interactions_per_time <- unique(edge_data_time_full$Interaction)
      for (interaction in interactions_per_time) {
        edge_data_time <- edge_data_time_full %>% filter(Interaction == interaction)

        # check if there is enough data for that interaction at that time
        if (nrow(edge_data_time) < 10 & method == "lm") {
          print(paste("Skipping linear regression for interaction", interaction, "at time", time_level,
                      "because data has less than 10 observations"))
          next
        }

        if (nrow(edge_data_time) < 20 & method == "logit") {
          print(paste("Skipping logistic regression for interaction", interaction, "at time", time_level,
                      "because data has less than 20 observations"))
          next
        }

        if (nrow(edge_data_time) < 30 & method == "multinom") {
          print(paste("Skipping multinomial for interaction", interaction, "at time", time_level,
                      "because data has less than 30 observations"))
          next
        }

        # Run the models here
        # 1. Linear Model
        if (method == "lm") {
          formula <- as.formula(paste(dependent, "~ Edge_weight"))
          set.seed(12345)
          model <- lm(formula, data = edge_data_time)
          model_summary <- broom::tidy(model)
          model_summary$Interaction <- interaction
          model_summary$Time_Level <- time_level
          model_summary$Model_Type <- "Linear Regression"
          model_summary$Model_SampleSize <- stats::nobs(model)
          node_results[[paste(interaction, time_level,
                              sep = "_")]] <- model_summary
        }


        # 2. Logistic Model
        if (method == "logit") {
          # Count unique levels
          unique_levels <- length(unique(edge_data_time[[dependent]]))
          # Skip if more than 2 levels
          if (unique_levels != 2) {
            print(paste("Skipping interaction", interaction,
                        "because dependent variable has", unique_levels,
                        "levels. Logistic regression requires exactly 2. Consider using multinomial for more than two or filtering to timepoints that have 2 options."))
            next } # Skip to the next interaction

          # Proceed with logistic regression if valid
          formula <- as.formula(paste(dependent, "~ Edge_weight"))
          set.seed(12345)
          model <- glm(formula, data = edge_data_time, family = binomial(link = "logit"))
          model_summary <- broom::tidy(model)
          model_summary$Interaction <- interaction
          model_summary$Time_Level <- time_level
          model_summary$Model_Type <- "Logistic Regression"
          model_summary$Model_SampleSize <- stats::nobs(model)
          node_results[[paste(interaction, time_level, sep = "_")]] <- model_summary
        }



        # 3. Multinomial Model
        if (method == "multinom") {
          # Skip if less than 2 levels
          # Count unique levels of the dependent variable
          edge_data_time[[dependent]] <- as.factor(edge_data_time[[dependent]])
          unique_levels <- length(unique(edge_data_time[[dependent]]))
          if (unique_levels <= 2) {
            print(paste("Skipping interaction", interaction,
                        "because dependent variable has less than two levels and multinomial requires 3. Consider using logistic regression with binary outcome"))
            next }  # Skip to the next interaction
          formula <- as.formula(paste(dependent, "~Edge_weight"))
          set.seed(12345)
          model <- nnet::multinom(formula, data = edge_data_time, trace = FALSE)
          model_summary <- broom::tidy(model, conf.int = TRUE)
          model_summary$Interaction <- interaction
          model_summary$Time_Level <- time_level
          model_summary$Model_Type <- "Multinomial"
          model_summary$Model_SampleSize <- stats::nobs(model)
          node_results[[paste(interaction, time_level, sep = "_")]] <- model_summary
        }
      }
    }

    # Summarize the Data
    Nodes <- bind_rows(node_results)
    Nodes <- Nodes %>% dplyr::filter(term == "Edge_weight")

    # Adjust p-value for all models
    Nodes <- Nodes %>% dplyr::mutate(p.adjusted = stats::p.adjust(p.value, method = "BH"))
    # Filter by p.value and or adjusted p value
    Nodes <- Nodes %>% dplyr::filter(p.value <=  pval) %>%
      dplyr::filter(p.adjusted <=  adjpval)


    # Finish summarizing the data
    sig_interactions <- Nodes$Interaction
    sig_edges <- Edge_Table %>% filter(Interaction %in% sig_interactions)
    unique_nodes <- unique(c(sig_edges$Source, sig_edges$Sink))
    Corrected_Counts <- LIMON_obj[["Corrected_Counts"]]
    SigNodes_Data <- Corrected_Counts[, unique_nodes, drop = FALSE]
    sample_data <- LIMON_obj[["SampleData"]]
    sample_data <- sample_data %>% dplyr::select(all_of(time),
                                                 all_of(dependent))
    SigNodes_Data <- SigNodes_Data %>% rownames_to_column(var = "SampleID")
    sample_data <- sample_data %>% rownames_to_column(var = "SampleID")
    SigNodes_Data <- merge(sample_data, SigNodes_Data, by = "SampleID",
                           all = TRUE)
    SigNodes_Data <- SigNodes_Data %>% column_to_rownames("SampleID")
  }

  # Return error message if they made a different request
  else {
    print("Only takes Edge Weights right now")
    Nodes <- NULL
    SigNodes_Data <- NULL
  }

  # Return Data
  LIMON_obj$Significant_Interactions <- Nodes
  LIMON_obj$Nodes_data <- SigNodes_Data


  # Plot the nodes ONLY for linear or logistic models
  if (plot_nodes == TRUE) {
    if (method %in% c("lm", "logit")) {
      sig_nodes <- as.data.frame(LIMON_obj[["Significant_Interactions"]]) %>%
        filter(term == "Edge_weight") %>%
        filter(p.value <= pval) %>%
        filter(abs(estimate) >= estimate)

      # Define color scale based on custom_colors
      if (is.null(custom_colors)) {
        color_scale <- scale_color_distiller(palette = "Set3", name = "Time")
      } else if (is.character(custom_colors) && length(custom_colors) == 1) {
        color_scale <- scale_color_distiller(palette = custom_colors, name = "Time")
      } else if (is.vector(custom_colors)) {
        color_scale <- scale_color_gradientn(colors = custom_colors, name = "Time")
      } else {
        stop("Invalid custom_colors argument. Provide a palette name or a named vector of colors.")
      }

      # Generate plot
      plot <- sig_nodes %>%
        arrange(estimate) %>%
        mutate(Interaction = factor(Interaction, levels = unique(Interaction[order(estimate)]))) %>%
        ggplot(aes(x = Interaction, y = estimate, color = as.numeric(Time_Level))) +
        geom_segment(aes(x = Interaction, xend = Interaction, y = 0, yend = estimate), color = "gray") +
        geom_point(size = point_size) +
        coord_flip() + theme_bw() + xlab("Interaction") +
        ylab("Estimate") + ggtitle("Effect Size of Microbial Interactions per Timepoint") +
        color_scale + theme(axis.text.x = element_text(color = "black", family = "Arial", size = 11),
                            axis.text.y = element_text(color = "black", family = "Arial", size = 11))

      print(plot)

    } else if (method == "multinom") {
      stop("Plotting function only available for Linear and Logistic Regression results, not Multinomial. Refer to LIMON tutorial for examples of how to plot these findings.")
    }
  }

  return(LIMON_obj)
}
