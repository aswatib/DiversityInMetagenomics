#         ESTIMATION OF ENTROPIES FOR SYNTHETIC DATA
#       FOR BOTH MAPPED COUNTS AND ALIGNMENT COUNTS
# ------------------------------------------------------------------------------
source("SYNTHETIC _MAPPED.R")
# ------------------------------------------------------------------------------
# ==============================================================================
#       INDIVIDUAL SAMPLE ENTROPIES BEFORE REMOVING VARIANCE TOWERS 
# ==============================================================================
# Define dataset names and corresponding data frames
datasets <- c("BacOriginal_Dirichlet_alpha10", "BacOriginal_Dirichlet_alpha1",
              "BacOriginal_sparseDirichlet", "BacOriginal_sparseDirichlet2",
              "BacATG_Dirichlet_alpha10", "BacATG_Dirichlet_alpha1", 
              "BacATG_sparseDirichlet", "BacATG_sparseDirichlet2")

# INDIVIDUAL ENTROPIES FOR EACH SAMPLE
# Loop over each dataset and plot manually -- for BacOriginal
for (i in 1:length(datasets)) {
  df <- get(datasets[i])
  df <- df[-1]
  print(dim(df))
  entropies_df <- estimate_individual_entropies(df)
  
  # Save the entropy dataframe for each set
  assign(paste0(datasets[i], "_individual_entropies"), entropies_df, envir = .GlobalEnv)
}

# ------------------------------------------------------------------------------
# Calculate average for all entropy methods across samples
combined_averages <- data.frame(
  BacATG_Set_A = colMeans(BacATG_Dirichlet_alpha10_individual_entropies),
  BacATG_Set_B = colMeans(BacATG_Dirichlet_alpha1_individual_entropies),
  BacATG_Set_C = colMeans(BacATG_sparseDirichlet_individual_entropies),
  BacATG_Set_D = colMeans(BacATG_sparseDirichlet2_individual_entropies),
  BacOriginal_Set_A = colMeans(BacOriginal_Dirichlet_alpha10_individual_entropies),
  BacOriginal_Set_B = colMeans(BacOriginal_Dirichlet_alpha1_individual_entropies),
  BacOriginal_Set_C = colMeans(BacOriginal_sparseDirichlet_individual_entropies),
  BacOriginal_Set_D = colMeans(BacOriginal_sparseDirichlet2_individual_entropies)
)

# Take averages for each entropy method across samples
combined_averages$Method <- rownames(combined_averages)

# Reshape the data to long format
long_averages <- combined_averages %>%
  pivot_longer(cols = -Method, names_to = "DataFrame", values_to = "Average_Entropy") %>%
  mutate(
    Set = gsub(".*_Set_", "", DataFrame),
    Database = gsub("_Set_.*", "", DataFrame)
  ) %>%
  mutate(Database = factor(Database), Set = factor(Set, levels = c("A", "B", "C", "D")))


# Define colors for the sets
colors <- c("A" = "darkblue", "B" = "royalblue",  "C" = "dodgerblue", "D" = "deepskyblue")

# Define shapes for the sets
shapes <- c("A" = 19, "B" = 11, "C" = 17, "D" = 15)

# Add newlines to method names for better readability
long_averages$Method <- gsub("\\.", "\n", long_averages$Method)

# Create the plot
average_plot <- ggplot(long_averages, aes(x = Method, y = Average_Entropy, color = Set, shape = Set, group = Set)) +
  geom_line(linetype = "solid", size = 0.5) +
  geom_point(size=3) +
  geom_hline(yintercept = 3.688879, linetype = "dashed", color = "darkorange", size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  facet_grid(~ Database) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16, vjust = -1),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0, hjust = 0.5),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "none",
    strip.text = element_text(size = 16) 
  ) +
  labs(
    title = "Average Entropy Values by Method",
    subtitle = "True entropy = 3.688879 (as natural log of the number of taxa, that is 40, present in each sample)\n\nBefore removing the variance towers",
    x = "Method",
    y = "Average Entropy values",
    color = "Set",
    shape = "Set"
  )

# Print the plot
print(average_plot)
# ------------------------------------------------------------------------------
# ==============================================================================
# INDIVIDUAL SAMPLE ENTROPIES AFTER REMOVING VARIANCE TOWERS 
# ==============================================================================
# Create vector of dataset names
datasets_after_var <- c("BacOriginal_Dirichlet_alpha10_after_var", 
                        "BacOriginal_Dirichlet_alpha1_after_var",
                        "BacOriginal_sparseDirichlet_after_var", 
                        "BacOriginal_sparseDirichlet2_after_var",
                        "BacATG_Dirichlet_alpha10_after_var", 
                        "BacATG_Dirichlet_alpha1_after_var",
                        "BacATG_sparseDirichlet_after_var", 
                        "BacATG_sparseDirichlet2_after_var")

# Loop over each dataset to estimate entropies
for (i in 1:length(datasets_after_var)) {
  df <- get(datasets_after_var[i])
  df <- df[-1]
  print(dim(df))
  entropies_df <- estimate_individual_entropies(df)
  
  # Save the entropy dataframe for each set
  assign(paste0(datasets_after_var[i], "_individual_entropies"), entropies_df, envir = .GlobalEnv)
}

# Calculate average for all entropy methods
combined_averages_after_var <- data.frame(
  BacATG_Set_A = colMeans(BacATG_Dirichlet_alpha10_after_var_individual_entropies),
  BacATG_Set_B = colMeans(BacATG_Dirichlet_alpha1_after_var_individual_entropies),
  BacATG_Set_C = colMeans(BacATG_sparseDirichlet_after_var_individual_entropies),
  BacATG_Set_D = colMeans(BacATG_sparseDirichlet2_after_var_individual_entropies),
  BacOriginal_Set_A = colMeans(BacOriginal_Dirichlet_alpha10_after_var_individual_entropies),
  BacOriginal_Set_B = colMeans(BacOriginal_Dirichlet_alpha1_after_var_individual_entropies),
  BacOriginal_Set_C = colMeans(BacOriginal_sparseDirichlet_after_var_individual_entropies),
  BacOriginal_Set_D = colMeans(BacOriginal_sparseDirichlet2_after_var_individual_entropies)
)

# Add a method column
combined_averages_after_var$Method <- rownames(combined_averages_after_var)

# Reshape the data to long format
long_averages_after_var <- combined_averages_after_var %>%
  pivot_longer(cols = -Method, names_to = "DataFrame", values_to = "Average_Entropy") %>%
  mutate(
    Set = gsub(".*_Set_", "", DataFrame),
    Database = gsub("_Set_.*", "", DataFrame)
  ) %>%
  mutate(Database = factor(Database), Set = factor(Set, levels = c("A", "B", "C", "D")))

# Add newlines to method names for better readability
long_averages_after_var$Method <- gsub("\\.", "\n", long_averages_after_var$Method)


# Create the plot
average_plot_after_var <- ggplot(long_averages_after_var, aes(x = Method, y = Average_Entropy, color = Set, shape = Set, group = Set)) +
  geom_line(linetype = "solid", size = 0.5) +
  geom_point(size=3) +
  geom_hline(yintercept = 3.688879, linetype = "dashed", color = "darkorange", size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  facet_grid(~ Database) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16, vjust = -2),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0, hjust = 0.5),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    strip.text = element_text(size = 16) 
  ) +
  labs(
    subtitle = "After removing the variance towers",
    x = "Method",
    y = "Average Entropy values",
    color = "Set",
    shape = "Set"
  )

# Print the plots
print(average_plot_after_var)
average_plot / average_plot_after_var

# Saved to Images as "AvgEntropyComparison1.pdf"
# Saved to Images as "AvgEntropyComparison2.pdf"
# Saved to Images as "AvgEntropyComparison.pdf"
# ==============================================================================
#                 ENTROPIES BASED ON ALIGNMENT COUNTS
# ==============================================================================
# Create a vector of dataset names for alignment based abundances
datasets_aln <- c("BacOriginal_Dirichlet_alpha10_before_var_aln",
                  "BacOriginal_Dirichlet_alpha1_before_var_aln",
                  "BacOriginal_sparseDirichlet_before_var_aln",
                  "BacOriginal_sparseDirichlet2_before_var_aln",
                  "BacATG_Dirichlet_alpha10_before_var_aln",
                  "BacATG_Dirichlet_alpha1_before_var_aln",
                  "BacATG_sparseDirichlet_before_var_aln",
                  "BacATG_sparseDirichlet2_before_var_aln",
                  
                  "BacOriginal_Dirichlet_alpha10_after_var_aln",
                  "BacOriginal_Dirichlet_alpha1_after_var_aln",
                  "BacOriginal_sparseDirichlet_after_var_aln",
                  "BacOriginal_sparseDirichlet2_after_var_aln",
                  
                  "BacATG_Dirichlet_alpha10_after_var_aln",
                  "BacATG_Dirichlet_alpha1_after_var_aln",
                  "BacATG_sparseDirichlet_after_var_aln",
                  "BacATG_sparseDirichlet2_after_var_aln")

# Estimate sample-wise entropy for each set
for (i in 1:length(datasets_aln)) {
  df <- get(datasets_aln[i])
  df <- df[-1]
  #alpha <- alphas[i]
  print(dim(df))
  entropies_df <- estimate_individual_entropies(df)
  
  # Save the entropy dataframe for each set
  assign(paste0(datasets_aln[i], "_individual_entropies"), entropies_df, envir = .GlobalEnv)
}

# BacOriginal alignment counts
combined_averages_aln <- data.frame(
  BacOriginal_Set_A_aln = colMeans(BacOriginal_Dirichlet_alpha10_before_var_aln_individual_entropies),
  BacOriginal_Set_B_aln = colMeans(BacOriginal_Dirichlet_alpha1_before_var_aln_individual_entropies),
  BacOriginal_Set_C_aln = colMeans(BacOriginal_sparseDirichlet_before_var_aln_individual_entropies),
  BacOriginal_Set_D_aln = colMeans(BacOriginal_sparseDirichlet2_before_var_aln_individual_entropies),
  
  BacATG_Set_A_aln = colMeans(BacATG_Dirichlet_alpha10_before_var_aln_individual_entropies),
  BacATG_Set_B_aln = colMeans(BacATG_Dirichlet_alpha1_before_var_aln_individual_entropies),
  BacATG_Set_C_aln = colMeans(BacATG_sparseDirichlet_before_var_aln_individual_entropies),
  BacATG_Set_D_aln = colMeans(BacATG_sparseDirichlet2_before_var_aln_individual_entropies)
)


# BacATG alignment counts
combined_averages_after_var_aln <- data.frame(
  BacOriginal_Set_A_after_var_aln = colMeans(BacOriginal_Dirichlet_alpha10_after_var_aln_individual_entropies),
  BacOriginal_Set_B_after_var_aln = colMeans(BacOriginal_Dirichlet_alpha1_after_var_aln_individual_entropies),
  BacOriginal_Set_C_after_var_aln = colMeans(BacOriginal_sparseDirichlet_after_var_aln_individual_entropies),
  BacOriginal_Set_D_after_var_aln = colMeans(BacOriginal_sparseDirichlet2_after_var_aln_individual_entropies),
  
  BacATG_Set_A_after_var_aln = colMeans(BacATG_Dirichlet_alpha10_after_var_aln_individual_entropies),
  BacATG_Set_B_after_var_aln = colMeans(BacATG_Dirichlet_alpha1_after_var_aln_individual_entropies),
  BacATG_Set_C_after_var_aln = colMeans(BacATG_sparseDirichlet_after_var_aln_individual_entropies),
  BacATG_Set_D_after_var_aln = colMeans(BacATG_sparseDirichlet2_after_var_aln_individual_entropies)
)

# Add method column to entropies for alignment counts with variance towers
combined_averages_aln$Method <- rownames(combined_averages_aln)

# Reshape the data to long format
long_averages_aln <- combined_averages_aln %>%
  pivot_longer(cols = -Method, names_to = "DataFrame", values_to = "Average_Entropy") %>%
  mutate(
    Set = gsub(".*_Set_(.*)_aln", "\\1", DataFrame),
    Database = gsub("_Set_.*", "", DataFrame)
  ) %>%
  mutate(Database = factor(Database), Set = factor(Set, levels = c("A", "B", "C", "D")))


# Add newlines to method names for better readability
long_averages_aln$Method <- gsub("\\.", "\n", long_averages_aln$Method)


# Create the plot
colors <- c("A" = "darkblue", "B" = "royalblue",  "C" = "dodgerblue", "D" = "deepskyblue")
average_plot_aln <- ggplot(long_averages_aln, aes(x = Method, 
                                                  y = Average_Entropy, 
                                                  color = Set, 
                                                  shape = Set, 
                                                  group = Set)) +
  geom_line(linetype = "solid", size = 0.5) +
  geom_point(size=3) +
  geom_hline(yintercept = 3.688879, linetype = "dashed", color = "darkorange", size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  facet_grid(~ Database) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16, vjust = -2),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0, hjust = 0.5),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "none",
    strip.text = element_text(size = 16) 
  ) +
  labs(title="Average Entropy Values By Method Using Alignment Counts",
       subtitle = "Before removing the variance towers",
       x = "Method",
       y = "Average Entropy values",
       color = "Set",
       shape = "Set"
  )

# Print the plot
print(average_plot_aln)
# ------------------------------------------------------------------------------
# Add method column to entropies for alignment counts without variance towers
combined_averages_after_var_aln$Method <- rownames(combined_averages_after_var_aln)

# Reshape the data to long format
long_averages_after_var_aln <- combined_averages_after_var_aln %>%
  pivot_longer(cols = -Method, names_to = "DataFrame", values_to = "Average_Entropy") %>%
  mutate(
    Set = gsub(".*_Set_(.*)_after_var_aln", "\\1", DataFrame),
    Database = gsub("_Set_.*", "", DataFrame)
  ) %>%
  mutate(Database = factor(Database), Set = factor(Set, levels = c("A", "B", "C", "D")))


# Add newlines to method names for better readability
long_averages_after_var_aln$Method <- gsub("\\.", "\n", long_averages_after_var_aln$Method)


# Create the plot
colors <- c("A" = "darkblue", "B" = "royalblue",  "C" = "dodgerblue", "D" = "deepskyblue")
average_plot_after_var_aln <- ggplot(long_averages_after_var_aln, aes(x = Method, 
                                                                      y = Average_Entropy, 
                                                                      color = Set, 
                                                                      shape = Set, 
                                                                      group = Set)) +
  geom_line(linetype = "solid", size = 0.5) +
  geom_point(size=3) +
  geom_hline(yintercept = 3.688879, linetype = "dashed", color = "darkorange", size = 1) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  facet_grid(~ Database) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title.x = element_text(size = 16, vjust = -2),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 15, angle = 0, vjust = 0, hjust = 0.5),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    strip.text = element_text(size = 16) 
  ) +
  labs(
    subtitle = "After removing the variance towers",
    x = "Method",
    y = "Average Entropy values",
    color = "Set",
    shape = "Set"
  )

# Print the plot
print(average_plot_after_var_aln)

average_plot_aln / average_plot_after_var_aln

# Saved to Images as "AvgEntropyComparison_Alignment.pdf"