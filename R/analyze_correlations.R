analyze_correlations <- function(targets_data) {
  cat("Analyzing correlations between target variables...\n")
  target_vars <- c("fat_percent", "protein_percent", "scc_thous_per_ml", "lactose_percent")
  cor_matrix <- targets_data %>%
    select(all_of(target_vars)) %>%
    cor(use = "complete.obs")
  
  cat("Correlation Matrix:\n")
  print(round(cor_matrix, 3))
  cat("\n")
  
  corrplot(
    cor_matrix, 
    method = "color", 
    type = "upper", 
    order = "hclust",
    tl.cex = 0.8, 
    tl.col = "black",
    addCoef.col = "black",
    number.cex = 0.8,
    title = "Correlation Matrix of Target Variables",
    mar = c(0,0,2,0)
  )
  
  pairs_plot <- targets_data %>%
    select(all_of(target_vars)) %>%
    ggpairs(
      columnLabels = c("Fat %", "Protein %", "SCC (1000/ml)", "Lactose %"),
      title = "Pairwise Relationships Between Target Variables"
    ) +
    theme_minimal()
  
  print(pairs_plot)
  return(cor_matrix)
}