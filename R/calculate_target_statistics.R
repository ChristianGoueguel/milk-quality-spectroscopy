calculate_target_statistics <- function(nir_data) {
  cat("=== Descriptive Statistics for Milk Quality Targets ===\n\n")
  targets_data <- nir_data$lab_results %>%
    select(
      fat_percent, protein_percent, scc_thous_per_ml, lactose_percent,
      barnname, lact, grp, dim_days, daily_milking, milk_wt_lbs) %>%
    filter(
      !is.na(fat_percent), 
      !is.na(protein_percent), 
      !is.na(scc_thous_per_ml),
      !is.na(lactose_percent)
    )
  
  cat("Total samples with complete target data:", nrow(targets_data), "\n\n")
  target_vars <- c("fat_percent", "protein_percent", "scc_thous_per_ml", "lactose_percent")
  
  summary_stats <- targets_data %>%
    select(all_of(target_vars)) %>%
    summarise(
      across(everything(), list(
        n = ~ sum(!is.na(.x)),
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        min = ~ min(.x, na.rm = TRUE),
        max = ~ max(.x, na.rm = TRUE),
        q25 = ~ quantile(.x, 0.25, na.rm = TRUE),
        q75 = ~ quantile(.x, 0.75, na.rm = TRUE),
        iqr = ~ IQR(.x, na.rm = TRUE),
        cv = ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE) * 100,
        skewness = ~ moments::skewness(.x, na.rm = TRUE),
        kurtosis = ~ moments::kurtosis(.x, na.rm = TRUE)
      ), .names = "{.col}_{.fn}")
    ) %>%
    pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
    separate(stat, into = c("variable", "statistic"), sep = "_(?=[^_]*$)") %>%
    pivot_wider(names_from = statistic, values_from = value) %>%
    select(variable, n, mean, median, sd, cv, min, q25, q75, max, iqr, skewness, kurtosis)
  
  cat("Summary Statistics Table:\n")
  print(kable(summary_stats, digits = 3, format = "markdown"))
  cat("\n")
  
  return(list(
    summary_stats = summary_stats,
    targets_data = targets_data
  ))
}
