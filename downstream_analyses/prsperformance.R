# prsperformance.R
# This script analyzes the output from dprism_t2dprs()

library(dplyr)
library(pROC)
library(fmsb)
library(ggplot2)
library(gridExtra)

#' Test PRS performance with Plots
#'
#' @param data Input dataframe.
#' @param outcome_col Name of the phenotype column (will be renamed to t2d).
#' @param covariates Vector of covariate names.
#' @param predictor_col Name of the main predictor (default "zmetascore").
#' @param output_prefix String. Prefix for saved plot filenames.
#' @param log_file String. Filename for the tracking log (default "model_tracking_log.txt").
#'
#' @return A dataframe of model statistics. Side effects: saves .png plots and appends to log file.
prsperformance <- function(data, outcome_col, covariates, 
                           predictor_col = "zmetascore", id_col = "IID", 
                           output_prefix = "Results", 
                           log_file = "model_tracking_log.txt") {
  
  # ---------------- 1. DATA PREP ----------------
  mod_data <- data %>%
    rename(t2d = all_of(outcome_col)) %>%
    select(all_of(id_col), t2d, all_of(predictor_col), all_of(covariates)) %>% 
    na.omit()
  
  if (!all(mod_data$t2d %in% c(0, 1))) stop("Outcome must be 0/1.")
  
  # Sample Sizes
  n_total <- nrow(mod_data)
  n_cases <- sum(mod_data$t2d == 1)
  n_controls <- sum(mod_data$t2d == 0)
  sample_size_label <- paste0("N = ", n_total, " (Cases: ", n_cases, ", Controls: ", n_controls, ")")
  
  # ---------------- 2. BASELINE MODEL ----------------
  base_formula <- as.formula(paste("t2d ~", paste(covariates, collapse = " + ")))
  baseline.model <- glm(base_formula, data = mod_data, family = binomial)
  
  mod_data$phat0 <- predict(baseline.model, type = 'response')
  baseline.ROC <- pROC::roc(mod_data$t2d, mod_data$phat0, quiet = TRUE)
  baseline.AUC <- as.numeric(baseline.ROC$auc)
  
  # AUC 95% CI
  base_ci <- ci.auc(baseline.ROC)
  base_auc_str <- paste0(round(baseline.AUC, 3), " (", round(base_ci[1], 3), "-", round(base_ci[3], 3), ")")
  
  baseline.R2 <- NagelkerkeR2(baseline.model)$R2
  
  # ---------------- 3. FULL MODEL (PRS) ----------------
  prs_formula <- as.formula(paste("t2d ~", predictor_col, "+", paste(covariates, collapse = " + ")))
  PRS.model <- glm(prs_formula, data = mod_data, family = binomial)
  
  mod_data$phat1 <- predict(PRS.model, type = 'response')
  PRS.ROC <- pROC::roc(mod_data$t2d, mod_data$phat1, quiet = TRUE)
  PRS.AUC <- as.numeric(PRS.ROC$auc)
  
  # AUC 95% CI
  prs_ci <- ci.auc(PRS.ROC)
  prs_auc_str <- paste0(round(PRS.AUC, 3), " (", round(prs_ci[1], 3), "-", round(prs_ci[3], 3), ")")
  
  PRS.R2 <- NagelkerkeR2(PRS.model)$R2
  
  # Incremental Stats
  PRS.AUCincr <- PRS.AUC - baseline.AUC
  PRS.R2incr <- PRS.R2 - baseline.R2
  
  # Coefficients
  coef_summary <- summary(PRS.model)$coefficients
  beta <- coef_summary[predictor_col, "Estimate"]
  se   <- coef_summary[predictor_col, "Std. Error"]
  p    <- coef_summary[predictor_col, "Pr(>|z|)"]
  z_score <- coef_summary[predictor_col, "z value"]
  log10p <- (pnorm(-abs(z_score), log.p = TRUE) + log(2)) / log(10) * -1
  or_str <- paste0(round(exp(beta), 2), " (", 
                   round(exp(beta - 1.96*se), 2), "-", 
                   round(exp(beta + 1.96*se), 2), ")")
  
  # ---------------- 4. PLOTTING ----------------
  
  # A. Density Plots
  plot_df <- mod_data
  plot_df$Status <- factor(plot_df$t2d, levels = c(0, 1), labels = c("Controls", "Cases"))
  
  p1 <- ggplot(plot_df, aes_string(x = predictor_col)) +
    geom_density(fill = "grey", alpha = 0.5) +
    labs(title = paste(predictor_col, "overall"), subtitle = paste0("Total N = ", n_total)) +
    theme_minimal()
  
  p2 <- ggplot(plot_df, aes_string(x = predictor_col, fill = "Status")) +
    geom_density(alpha = 0.5) +
    labs(title = paste(predictor_col, "by Status"), subtitle = paste0("Cases: ", n_cases, " | Controls: ", n_controls)) +
    scale_fill_manual(values = c("Controls" = "#00BFC4", "Cases" = "#F8766D")) +
    theme_minimal() + theme(legend.position = "bottom")
  
  ggsave(filename = paste0(output_prefix, "_Density.png"), plot = grid.arrange(p1, p2, nrow = 2), width = 7, height = 7)
  
  # B. ROC Plot with CIs in Legend
  png(filename = paste0(output_prefix, "_ROC.png"), width = 2000, height = 2000, res = 300)
  pROC::plot.roc(baseline.ROC, col = "grey20", lty = 2, main = paste("ROC:", output_prefix))
  pROC::plot.roc(PRS.ROC, col = "red", lty = 1, add = TRUE)
  mtext(sample_size_label, side = 3, line = 0.5, cex = 0.8)
  
  legend("bottomright", 
         legend = c(paste0("Base AUC: ", base_auc_str), 
                    paste0("Full AUC: ", prs_auc_str)),
         col = c("grey20", "red"), lty = c(2, 1), lwd = 2, bty = "n", cex = 0.9)
  dev.off()
  
  # ---------------- 5. LOGGING ----------------
  
  # Create a timestamp
  rundate <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Create a single line log entry
  log_entry <- paste(
    rundate,
    "| Run:", output_prefix,
    "| Outcome:", outcome_col,
    "| Predictor:", predictor_col,
    "| Covariates:", paste(covariates, collapse=","),
    "| N_Total:", n_total,
    "| Base_AUC:", base_auc_str,
    "| Full_AUC:", prs_auc_str,
    "| OR:", or_str,
    "| P-val:", formatC(p, format = "e", digits = 2),
    "\n"
  )
  
  # Check if file exists, if not create header
  if (!file.exists(log_file)) {
    cat("Timestamp | Run_ID | Outcome | Predictor | Covariates | N | Base_AUC_95CI | Full_AUC_95CI | OR_95CI | P_Val\n", file = log_file)
  }
  
  # Append entry to file
  cat(log_entry, file = log_file, append = TRUE)
  
  # ---------------- 6. RETURN RESULTS ----------------
  results <- data.frame(
    n_total = n_total, n_cases = n_cases, n_controls = n_controls,
    beta = beta, se = se, or = or_str, p = p, log10p = log10p,
    baseline.AUC = baseline.AUC, baseline.AUC.low = base_ci[1], baseline.AUC.up = base_ci[3],
    PRS.AUC = PRS.AUC, PRS.AUC.low = prs_ci[1], PRS.AUC.up = prs_ci[3],
    PRS.AUCincr = PRS.AUCincr,
    baseline.R2 = baseline.R2, PRS.R2 = PRS.R2, PRS.R2incr = PRS.R2incr,
    stringsAsFactors = FALSE
  )
  
  return(results)
}
