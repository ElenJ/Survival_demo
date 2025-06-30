#' Plot log(-log) survival curves for all covariates
#'
#' For each variable in \code{covariates}, this function plots log(-log) Kaplan-Meier survival curves.
#' Continuous (numeric, >2 unique values) variables are dichotomized at the median before plotting.
#' Categorical or binary variables are used as-is. Each plot is displayed in turn.
#'
#' @param covariates Character vector. Column names of covariates to plot (excluding time/event).
#' @param df Data frame. Must contain the covariates and the time/status columns referenced in \code{SurvObj}.
#' @param SurvObj Character string. Expression for the survival object, e.g. \code{"Surv(time, status)"}.
#'
#' @return Invisibly returns NULL (plots are produced as side-effect).
#' @importFrom survminer surv_fit ggsurvplot
#' @importFrom stats as.formula median
#' @examples
#' plot_all_vars_loglog(c("age", "sex"), lung, "Surv(time, status)")
plot_all_vars_loglog <- function(covariates, df, SurvObj){
  # cox model with all covariates
  cox_formula <- as.formula(paste0(SurvObj, "~", paste(covariates, collapse="+")))
  cox_model <- coxph(cox_formula, data=df)
  #perform a statistical test on the proportional hazards assumption
  ph_assumption_test <- cox.zph(cox_model, transform=rank)
  
  # plotting part
  for (i in 1:length(covariates)){
    #extract proportional hazards test pVal
    if(covariates[i] %in% rownames(ph_assumption_test$table)){
      ph_pval <- ph_assumption_test$table[covariates[i], "p"]
      ph_pval_text <- paste0("PH test p = ", signif(ph_pval, 3))
    } else {
      ph_pval_text = "PH test p = NA"
    }

    # If the variable is continuous:numeric, and not just 2 levels, dichotomize at median
    if (is.numeric(df[,covariates[i]]) && length(unique(df[,covariates[i]])) > 2) {
      print(covariates[i])
      medval <- median(df[[covariates[i]]], na.rm = TRUE)
      # Create a high/low grouping variable
      groupname <- paste0(covariates[i], "_group")
      df[[groupname]] <- ifelse(df[[covariates[i]]] <= medval, "Low", "High")
      df[[groupname]] <- factor(df[[groupname]], levels = c("Low", "High"))
      # Build the Surv formula
      myformula <- as.formula(paste0(SurvObj, " ~ ", groupname))
      kmfit <- surv_fit(formula = myformula, data = df)
      # Plot log(-log) survival for median split
      p <- ggsurvplot(
        kmfit,
        data = df,
        fun = "cloglog",
        legend.labs = c("Low", "High")
      ) + 
        ggtitle(paste("KM by", covariates[i], "(dichotomized at median)"))
      
      print(p$plot + 
        annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, 
                 label = ph_pval_text, size = 5, fontface = "italic"))

      # drop the added column to prevent it being replotted in the next elif
      df <- df[, !(colnames(df) %in% groupname)]
    } else if (!is.numeric(df[,covariates[i]]) | length(unique(df[,covariates[i]])) == 2) {
      print(covariates[i])
      # For binary or categorical variables, plot as-is
      myformula <- as.formula(paste0(SurvObj, " ~ ", covariates[i]))
      kmfit <- surv_fit(formula = myformula, data = df)

      # Plot log(-log) survival for median split
      p <- ggsurvplot(
        kmfit,
        data = df,
        fun = "cloglog",
        legend.labs = c("Low", "High")
      ) + 
        ggtitle(paste("KM by", covariates[i]))
      
      print(p$plot + 
              annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, 
                       label = ph_pval_text, size = 5, fontface = "italic"))
    }
  }
  invisible(NULL)
}


#' Select Covariates Based on PH Assumption via cox.zph
#'
#' Identifies and groups covariates from a Cox proportional hazards model as either satisfying or violating the proportional hazards (PH) assumption using the Schoenfeld residuals test.  
#' 
#' Internally calls `cox.zph()` to test the PH assumption for each covariate, returning those with p-values greater than 0.01 as non-violating and those with p-values less than or equal to 0.01 as violating. The "GLOBAL" test result, if significant, is excluded from the violating list.
#'
#' @param cox_model An object of class \code{coxph}, the result of fitting a Cox proportional hazards model.
#'
#' @return A named list containing:
#' \describe{
#'   \item{cofactors}{A list with two elements:}
#'   \itemize{
#'     \item{\code{non_violating}}{A character vector of covariate names satisfying the PH assumption (p > 0.01).}
#'     \item{\code{violating}}{A character vector of covariate names violating the PH assumption (p ≤ 0.01), excluding "GLOBAL".}
#'   }
#' }
#'
#' @details
#' The function performs the proportional hazards test using Schoenfeld residuals (with "rank" transformation) for each covariate of the fitted Cox model. The significance level threshold is set at 0.01.
#'
#' @seealso \code{\link[survival]{coxph}}, \code{\link[survival]{cox.zph}}
#' @export
PH_cofactors <- function(cox_model){
  # perform PH test of Schoenfeld residuals
  ph_test <- cox.zph(cox_model, transform = "rank")
  # convert to table and filter out non-violating factors
  ph_test <- as.data.frame(ph_test$table)
  ph_test$cofactor <- rownames(ph_test)
  non_violating_cofactors <- ph_test$cofactor[ph_test$p > 0.01]
  violating_cofactors <- ph_test$cofactor[ph_test$p <= 0.01]
  violating_cofactors <- setdiff(violating_cofactors, "GLOBAL") # drop "GLOBAL", if GLOBAL is significant
  # return list
  results <- list(
    cofactors = list(
      non_violating = non_violating_cofactors,
      violating = violating_cofactors
    )
  )
  return(results)
}

#' Calculate Univariate Cox Proportional Hazards Models for Multiple Features
#'
#' This function fits univariate Cox proportional hazards models for each variable specified in \code{mycovariates} and returns a summary of estimated hazard ratios, p-values, and concordance indices for each feature.
#'
#' @param df A data frame containing the survival data, including survival time (\code{survt}), event status (\code{status}), and all candidate covariates.
#' @param mycovariates A character vector of covariate (feature) names to be evaluated one by one in separate Cox models.
#'
#' @details
#' For each feature in \code{mycovariates}, the function fits a univariate Cox proportional hazards model of the form \code{Surv(survt, status == 1) ~ feature}. 
#' It extracts the hazard ratio (HR), p-value, and Harrell's concordance index for each model and returns the results in a data frame.
#'
#' @return A data frame with one row per feature, including:
#' \describe{
#'   \item{variable}{Feature name}
#'   \item{HR}{Estimated hazard ratio for the feature}
#'   \item{p_val}{P-value for the feature's effect}
#'   \item{concordance}{Concordance index (C-index) for the fitted model}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' library(survival)
#' result <- calculate_cox_per_feature(your_df, c("age", "sex", "biomarker"))
#' print(result)
#' }
#'
#' @import survival
#' @export
calculate_cox_per_feature <- function(df, mycovariates){
  # Initialize vectors to store hazard ratios (HR), p-values, and concordance indices for variables
  HR <- c()
  p_val <- c()
  concordance <- c()
  variable <- c()
  
  # Initialize vectors to store hazard ratios (HR) and p-values for competing risks regression
  HR_variable_comp <- c()
  p_variable_comp <- c()
  
  # Initialize vectors to store adjusted hazard ratios (HR) and p-values
  HR_variable_adj <- c()
  p_variable_adj <- c()
  
  # Loop through each variable in the quality-controlled list
  for(i in 1:length(mycovariates)) {
    # Prepare the data for survival analysis
    tmp <- df[c(mycovariates[i], "status", "survt")]
    variable <- c(variable, mycovariates[i])
    
    # Create a survival object
    tmp$Y <- with(tmp, Surv(survt, status == 1))
    
    # Fit a Cox proportional hazards model with the current variable
    cox.model <- coxph(Y ~ tmp[, 1], data = tmp)
    #print(tidy(cox.model))
    
    # Store the hazard ratio (HR), p-value, and concordance index for the current variable
    HR <- c(HR, summary(cox.model)$coefficient[2])
    p_val <- c(p_val, summary(cox.model)$coefficient[5])
    concordance <- c(concordance, summary(cox.model)$concordance[1])
    
  }
  # Combine the results into a data frame
  result_cox <- as.data.frame(cbind(variable, HR, p_val, concordance))
  result_cox <- result_cox %>%
    mutate(across(-variable, as.numeric))
  return(result_cox)
}

#' Run Backward Selection for Cox Proportional Hazards Model
#'
#' This function performs backward variable selection for a Cox proportional hazards model, starting from a set of pre-specified candidate variables. At each step, the variable with the highest p-value (above a Bonferroni-like cutoff) is removed, and the model is re-fit until all remaining variables have statistically significant p-values.
#'
#' @param df A data frame containing the data, with columns for survival time (\code{survt}), event status (\code{status}), and model variables.
#' @param modelVars A character vector of variable names (columns in \code{df}) to be considered for model selection.
#'
#' @details
#' The function (optionally) imputes missing values in the candidate variables using mean imputation (for numeric variables), constructs a survival object, and fits the initial Cox model with all candidate variables.  
#' In each iteration, the variable with the largest p-value above the cutoff (\code{0.01 / N}, where \code{N} is the number of candidate variables) is removed and the model is re-fit. This continues until all remaining variables have p-values below the cutoff, resulting in a parsimonious model.
#'
#' @return A character vector of variable names retained in the final Cox model after backward selection.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(dplyr)
#' selected_vars <- run_cox_backward_selection(mydata, c("age", "sex", "biomarker1", "biomarker2"))
#' }
#'
#' @import survival
#' @importFrom dplyr mutate across all_of
#' @importFrom broom tidy
#' @export
run_cox_backward_selection <- function(df, modelVars){
  master <- df # create copy to impute NA, if necessary
  
  # (Optional) Impute missing values for numeric modelVars, e.g. mean imputation
  master <- master |> mutate(across(all_of(modelVars), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  
  # Create survival object
  master <- master %>%
    mutate(Y = Surv(survt, status == 1))
  
  # Build Cox formula
  cox_formula_str <- paste("Y ~", paste(modelVars, collapse = " + "))
  cox_formula <- as.formula(cox_formula_str)
  
  # Fit Cox model
  cox_model <- coxph(cox_formula, data = master)
  summary(cox_model)
  
  # Tidy model output
  tidy_cox <- tidy(cox_model) %>% drop_na()
  
  # Set p-value cutoff (Bonferroni-like)
  pCut <- 0.01 / length(modelVars)
  
  # Backward selection loop
  currently_selected <- tidy_cox$term
  repeat {
    # Refit model with remaining variables
    cox_formula_str <- paste("Y ~", paste(currently_selected, collapse = " + "))
    cox_formula <- as.formula(cox_formula_str)
    cox_model <- coxph(cox_formula, data = master)
    tidy_cox <- tidy(cox_model) %>% drop_na()
    
    # Find max p-value
    max_pval <- max(tidy_cox$p.value, na.rm = TRUE)
    if (max_pval <= pCut) {
      final_cox_formula <- cox_formula
      break
    }
    
    # Remove variable with highest p-value
    var_to_remove <- tidy_cox$term[which.max(tidy_cox$p.value)]
    currently_selected <- setdiff(currently_selected, var_to_remove)
  }
  # Print final model formula
  cat("Your final model is:", deparse(final_cox_formula), "\n")
  return(currently_selected)
}
