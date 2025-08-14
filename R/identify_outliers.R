identify_outliers <- function(targets_data) {
  cat("Identifying outliers in target variables...\n")
  target_vars <- c("fat_percent", "protein_percent", "scc_thous_per_ml", "lactose_percent")
  
  outlier_summary <- map_dfr(target_vars, ~ {
    var_data <- targets_data[[.x]]
    q1 <- quantile(var_data, 0.25, na.rm = TRUE)
    q3 <- quantile(var_data, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr
    
    outliers <- which(var_data < lower_bound | var_data > upper_bound)
    
    tibble(
      variable = .x,
      n_outliers = length(outliers),
      outlier_percentage = round(length(outliers) / length(var_data) * 100, 2),
      lower_bound = round(lower_bound, 3),
      upper_bound = round(upper_bound, 3),
      min_outlier = ifelse(length(outliers) > 0, round(min(var_data[outliers]), 3), NA),
      max_outlier = ifelse(length(outliers) > 0, round(max(var_data[outliers]), 3), NA)
    )
  })
  
  cat("Outlier Summary (using IQR method):\n")
  print(kable(outlier_summary, format = "markdown"))
  cat("\n")
  return(outlier_summary)
}