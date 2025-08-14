analyze_by_groups <- function(targets_data) {
  cat("Analyzing target variables by categorical groups...\n")
  target_vars <- c("fat_percent", "protein_percent", "scc_thous_per_ml", "lactose_percent")
  by_lactation <- targets_data %>%
    group_by(lact) %>%
    summarise(
      n = n(),
      across(all_of(target_vars), list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      ), .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  cat("Summary by Lactation Number:\n")
  print(kable(by_lactation, digits = 3, format = "markdown"))
  cat("\n")
  
  by_group <- targets_data %>%
    group_by(grp) %>%
    summarise(
      n = n(),
      across(all_of(target_vars), list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      ), .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  cat("Summary by Group:\n")
  print(kable(by_group, digits = 3, format = "markdown"))
  cat("\n")
  
  violin_lact <- targets_data %>%
    select(lact, all_of(target_vars)) %>%
    pivot_longer(cols = all_of(target_vars), names_to = "variable", values_to = "value") %>%
    mutate(
      variable = str_replace(variable, "_", " ") %>% str_to_title(),
      lact = factor(lact)
    ) %>%
    ggplot(aes(x = lact, y = value, fill = lact)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Distribution of Target Variables by Lactation Number",
      x = "Lactation Number",
      y = "Value",
      fill = "Lactation"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
  
  print(violin_lact)
  
  violin_group <- targets_data %>%
    select(grp, all_of(target_vars)) %>%
    pivot_longer(cols = all_of(target_vars), names_to = "variable", values_to = "value") %>%
    mutate(variable = str_replace(variable, "_", " ") %>% str_to_title()) %>%
    ggplot(aes(x = grp, y = value, fill = grp)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.8) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    labs(
      title = "Distribution of Target Variables by Group",
      x = "Group",
      y = "Value",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))
  
  print(violin_group)
  
  return(list(
    by_lactation = by_lactation,
    by_group = by_group,
    violin_lact = violin_lact,
    violin_group = violin_group
  ))
}