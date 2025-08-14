create_distribution_plots <- function(targets_data) {
  cat("Creating distribution plots...\n")
  target_vars <- c("fat_percent", "protein_percent", "scc_thous_per_ml", "lactose_percent")
  hist_plots <- map(target_vars, ~ {
    targets_data %>%
      ggplot(aes(x = !!sym(.x))) +
      geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7, fill = "steelblue", color = "white") +
      geom_density(color = "red", size = 1, alpha = 0.7) +
      labs(
        title = paste("Distribution of", str_replace(.x, "_", " ") %>% str_to_title()),
        x = str_replace(.x, "_", " ") %>% str_to_title(),
        y = "Density"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
  })
  names(hist_plots) <- target_vars
  combined_hist <- wrap_plots(hist_plots, ncol = 2)
  print(combined_hist)
  box_plots <- targets_data %>%
    select(all_of(target_vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = str_replace(variable, "_", " ") %>% str_to_title()) %>%
    ggplot(aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    facet_wrap(~variable, scales = "free", ncol = 2) +
    labs(
      title = "Box Plots of Target Variables",
      x = "Variable",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold")
    )
  print(box_plots)
  return(list(histograms = hist_plots, boxplots = box_plots))
}