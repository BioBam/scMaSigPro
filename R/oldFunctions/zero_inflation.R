# Check for total zero inflation
total_sparsity <- function(mat, precision = 2) {
  sparsity_value <- coop::sparsity(as.matrix(mat)) * 100
  sparsity_value <- round(sparsity_value, precision)
  return(sparsity_value)
}

# Check induced simulation
induced_sparsity <- function(sim_mat, true_mat, precision = 2) {
  sim_sparsity_value <- coop::sparsity(as.matrix(sim_mat)) * 100
  true_sparsity_value <- coop::sparsity(as.matrix(true_mat)) * 100
  induced_sparsity_value <- round((sim_sparsity_value - true_sparsity_value), precision)
  return(induced_sparsity_value)
}

# Per column sparsity
per_col_zero <- function(sim_mat, precision = 2) {
  zero_vec <- round(colSums(sim_mat == 0), 2)
  non_zero_vec <- round(colSums(sim_mat != 0), 2)
  total <- rep(nrow(sim_mat), ncol(sim_mat))
  sparsity_per_col <- zero_vec / total * 100

  return(sparsity_per_col)
}

# Plot_zero distribution
plot_zero_density_per_cell <- function(sim_mat, precision = 2) {
  plt.table <- data.frame(
    Num = per_col_zero(sim_mat),
    cell_id = names(per_col_zero(sim_mat))
  )
  # Plot Data
  p <- ggplot(plt.table, aes(x = Num)) +
    geom_density(
      color = "#f68a53",
      fill = "#f68a53", alpha = 0.5
    ) +
    geom_vline(aes(xintercept = mean(Num)), linetype = "dashed", color = "#139289") +
    theme_classic() +
    theme(
      legend.position = "none", strip.text = element_text(size = rel(2)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(size = 0.7, linetype = "dotted"),
      panel.grid.minor = element_line(size = 0.2)
    ) +
    ggtitle("Density of Zeros per cell") +
    ylab("Number of cells") +
    xlab("Percentage of Zero")
  print(p)
}

# Plot_zero distribution
plot_zero_hist_per_cell <- function(sim_mat, precision = 2) {
  plt.table <- data.frame(
    Num = per_col_zero(sim_mat),
    cell_id = names(per_col_zero(sim_mat))
  )
  # Plot Data
  p <- ggplot(plt.table, aes(x = Num)) +
    geom_histogram(
      binwidth = 0.5, color = "#f68a53",
      fill = "#f68a53", alpha = 0.5
    ) +
    geom_vline(aes(xintercept = mean(Num)), linetype = "dashed", color = "#139289") +
    theme_classic() +
    theme(
      legend.position = "none", strip.text = element_text(size = rel(2)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(size = 0.7, linetype = "dotted"),
      panel.grid.minor = element_line(size = 0.2)
    ) +
    ggtitle("Distribution of Zeros per cell") +
    # scale_x_continuous(breaks = seq(0,30, by = 2))+
    ylab("Number of cells") +
    xlab("Percentage of Zero")
  print(p)
}
