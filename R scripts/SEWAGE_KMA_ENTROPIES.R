# ==============================================================================
#                 ENTROPY ESTIMATIONS FOR SEWAGE COUNT FROM KMA HITS
#                              WITH VARIANCE TOWERS
# ==============================================================================

# Load the KMA abundance data and estimate entropies
load("../RData Objects/kma_com.RData")
kma_entropies <- estimate_individual_entropies(kma_com)

# Save estimated entropies as RData object
save(kma_entropies, file = "../RData Objects/kma_entropies.RData")

# Load the object entropies object for KMA data when required
load("../RData Objects/kma_entropies.RData")
# ------------------------------------------------------------------------------
# ==============================================================================
#  DEVIATION OF SAMPLE-WISE ENTROPY VALUES FROM AVG ENTROPY FOR EACH METHOD
#                              WITH VARIANCE TOWERS
# ==============================================================================
# Calculate the absolute difference between entropy and the average entropy
entropies <- kma_entropies
avg_entropies <- colMeans(entropies)

# Subtract each value in the columns from the corresponding column mean and take the absolute difference
diff_entropies <- abs(sweep(entropies, 2, avg_entropies, "-"))
diff_entropies$Sample <- as.numeric(rownames(diff_entropies))

# Reshape the data to long format
long_diff_entropies <- diff_entropies %>%
  pivot_longer(cols = -Sample, names_to = "Method", values_to = "Difference")

long_diff_entropies$Sample <- as.factor(long_diff_entropies$Sample)
long_diff_entropies$Method <- as.factor(long_diff_entropies$Method)

# Define distinct colors for the methods
colors <- colorspace::sequential_hcl(6, "Blues", rev = TRUE)

# Reshape the data to long format
long_diff_entropies <- diff_entropies %>%
  pivot_longer(cols = -Sample, names_to = "Method", values_to = "Difference")

# Create the stacked bar plot using ggplot2
ggplot(long_diff_entropies, aes(x = factor(Sample), y = Difference, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Samples", y = "Absolute Differences between Entropy and Avg Entropy") +
  ggtitle(label = "Deviation of Sample-wise Entropy Values from Average Entropy for Each Method",
          subtitle = "Average entropy for each method is calculated as the mean entropy value across all samples.\n\nBased on KMA hits") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face="italic"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(vjust=0),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 11),
    legend.position = "none"
  ) +
  scale_fill_manual(values = colors)

# Saved to Images as "AvgEntropyDiffs_KMA.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#                 ENTROPY ESTIMATIONS FOR SEWAGE COUNT FROM KMA HITS
#                            WITHOUT VARIANCE TOWERS
# ==============================================================================
# Load the KMA abundance data without variance towers and estimate entropies
load("../RData Objects/kma_after_var.RData")
kma_after_var_entropies <- estimate_individual_entropies(kma_after_var)

# Save estimated entropies as RData object
save(kma_after_var_entropies, file = "../RData Objects/kma_after_var_entropies.RData")

# Load the object entropies object for KMA without variance towers data when required
load("../RData Objects/kma_after_var_entropies.RData")
# ------------------------------------------------------------------------------
# ==============================================================================
#  DEVIATION OF SAMPLE-WISE ENTROPY VALUES FROM AVG ENTROPY FOR EACH METHOD
#                            WITHOUT VARIANCE TOWERS
# ==============================================================================
# Deviation of Sample-wise Entropy Values from Average Entropy for Each Method
avg_entropies <- colMeans(kma_after_var_entropies)

# Subtract each value in the columns from the corresponding column mean and take the absolute difference
diff_entropies <- abs(sweep(kma_after_var_entropies, 2, avg_entropies, "-"))
diff_entropies$Sample <- as.numeric(rownames(diff_entropies))

# Reshape the data to long format
long_diff_entropies <- diff_entropies %>%
  pivot_longer(cols = -Sample, names_to = "Method", values_to = "Difference")

long_diff_entropies$Sample <- as.factor(long_diff_entropies$Sample)
long_diff_entropies$Method <- as.factor(long_diff_entropies$Method)

# Define distinct colors for the methods
colors <- colorspace::sequential_hcl(6, "Blues", rev = TRUE)

# Reshape the data to long format
long_diff_entropies <- diff_entropies %>%
  pivot_longer(cols = -Sample, names_to = "Method", values_to = "Difference")

# Create the stacked bar plot using ggplot2
ggplot(long_diff_entropies, aes(x = factor(Sample), y = Difference, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Samples", y = "Absolute Differences between Entropy and Avg Entropy") +
  ggtitle(label = "Deviation of Sample-wise Entropy Values from Average Entropy for Each Method",
          subtitle = "Average entropy for each method is calculated as the mean entropy value across all samples.\n\nBased on KMA counts without variance towers") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face="italic"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(vjust=-2),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  scale_fill_manual(values = colors)

# Saved to Images as "AvgEntropyDiffs_KMA_After_Var.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#      RELATIONSHIP BETWEEN AVERAGE ENTROPIES AND DISTANCE MEASURES
#           Only for Original KMA hits with variance towers
# ==============================================================================
# Calculate average distance of each sample from all other samples
df <- clr(kma_com)
dist_matrix <- vegan::vegdist(df, method = "euclidean")
dist_matrix_no_diag <- as.matrix(dist_matrix)
diag(dist_matrix_no_diag) <- NA # Putting NA for the diagonals so that zero values do not affect the avg calculation
avg_distance <- apply(dist_matrix_no_diag, 1, function(x) mean(x, na.rm = TRUE))

# Combine entropy and average distance into a data frame
modeling_data <- kma_entropies
modeling_data$Sample <- rownames(kma_entropies)
modeling_data$AvgDistance <- avg_distance

# Sort observations in the order of sample numbers
modeling_data <- modeling_data[order(as.numeric(modeling_data$Sample)), ]

# LINEAR MODEL
# Fit a linear model of average distance with every entropy method values
model_naive <- lm(Naive.entropy ~ AvgDistance, data = modeling_data)
model_chao <- lm(Chao.Shen ~ AvgDistance, data = modeling_data)
model_hausser <- lm(Hausser.Strimmer ~ AvgDistance, data = modeling_data)
model_wolpert <- lm(Wolpert.Wolf ~ AvgDistance, data = modeling_data)
model_nsb <- lm(NSB ~ AvgDistance, data = modeling_data)
model_piga <- lm(Piga ~ AvgDistance, data = modeling_data)

# Print model summaries
print(summary(model_naive))
print(summary(model_chao))
print(summary(model_hausser))
print(summary(model_wolpert))
print(summary(model_nsb))
print(summary(model_piga))

# Residuals vs Fitted Plots For All Methods
par(mfrow = c(3, 2))
plot(model_naive, which = 1, main = "Naive entropy ~ AvgDistance\n")
plot(model_chao, which = 1, main = "Chao-Shen entropy ~ AvgDistance\n")
plot(model_hausser, which = 1, main = "Hausser-Strimmer entropy ~ AvgDistance\n")
plot(model_wolpert, which = 1, main = "Wolpet-Wolf entropy ~ AvgDistance\n")
plot(model_nsb, which = 1, main = "NSB entropy ~ AvgDistance\n")
plot(model_piga, which = 1, main = "Piga's entropy ~ AvgDistance\n")

# Saved to Images as "residualVsFitted_KMA.pdf"
# ------------------------------------------------------------------------------
# CORRELATION ANALYSIS FOR KMA HITS

# Pearson correlation test
pearson_cor_df <- data.frame(Naive.entropy_avgDist = cor(modeling_data$AvgDistance, 
                                                         modeling_data$Naive.entropy, 
                                                         method="pearson"),
                             Chao.Shen_avgDist = cor(modeling_data$AvgDistance, 
                                                     modeling_data$Chao.Shen, 
                                                     method="pearson"),
                             Hausser.Strimmer_avgDist = cor(modeling_data$AvgDistance, 
                                                            modeling_data$Hausser.Strimmer, 
                                                            method="pearson"),
                             Wolpert.Wolf_avgDist = cor(modeling_data$AvgDistance, 
                                                        modeling_data$Wolpert.Wolf, 
                                                        method="pearson"),
                             NSB_avgDist = cor(modeling_data$AvgDistance, 
                                               modeling_data$NSB, 
                                               method="pearson"),
                             Piga_avgDist = cor(modeling_data$AvgDistance, 
                                                modeling_data$Piga, 
                                                method="pearson")
)

# Print Pearson correlation values
t(pearson_cor_df)


# Spearman correlation test
Spearman_cor_df <- data.frame(Naive.entropy_avgDist = cor(modeling_data$AvgDistance, 
                                                          modeling_data$Naive.entropy, 
                                                          method="spearman"),
                              Chao.Shen_avgDist = cor(modeling_data$AvgDistance, 
                                                      modeling_data$Chao.Shen, 
                                                      method="spearman"),
                              Hausser.Strimmer_avgDist = cor(modeling_data$AvgDistance, 
                                                             modeling_data$Hausser.Strimmer, 
                                                             method="spearman"),
                              Wolpert.Wolf_avgDist = cor(modeling_data$AvgDistance, 
                                                         modeling_data$Wolpert.Wolf, 
                                                         method="spearman"),
                              NSB_avgDist = cor(modeling_data$AvgDistance, 
                                                modeling_data$NSB, 
                                                method="spearman"),
                              Piga_avgDist = cor(modeling_data$AvgDistance, 
                                                 modeling_data$Piga, 
                                                 method="spearman")
)

# Print Spearman correlation values
t(Spearman_cor_df)
# ------------------------------------------------------------------------------
# ==============================================================================
