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
