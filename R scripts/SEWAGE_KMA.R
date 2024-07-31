# ==============================================================================
# --------- LOADING REQUIRED LIBRARIES AND AESTHETICS --------------------------
# ==============================================================================
source("PACKAGES.R")
source("CLEANER.R")
source("ALPHA_ESTIMATORS.R")
source("ENTROPY_ESTIMATORS.R")
# ==============================================================================
#                  DATA LOADING, CLEANING AND PROCESSING
# ==============================================================================
# ------------------------------------------------------------------------------
# ------------------------- KMA MAPSTAT DATA  ---------------------------------
# ------------------------------------------------------------------------------
# Set working directory
dir <- "/Users/swati/Desktop/Master Thesis/Programming/sewage_mapstats"
setwd(dir)
# ------------------------------------------------------------------------------
# Reading all the mapstat files in the directory
mapstat_df <- mapstatHelpers::read_multiple_mapstats(dir)

# Collapsing counts by sample id and then reshaping the dataset
reshaped_df <- mapstat_df %>%
  dplyr::select(sample_id, refSequence, fragmentCount) %>%
  pivot_wider(names_from = refSequence, values_from = fragmentCount, values_fill = 0)

reshaped_df <- data.frame(reshaped_df)
dim(reshaped_df)

# Clean and order the sewage abundance data
df <- clean_and_order(reshaped_df)
dim(df)

# Remove the sample column for diversity estimation and analysis
com <- df[-1]
kma_com <- com

# Save the KMA abundance set as RData object
save(kma_com, file = "kma_com.RData")
# write.csv(com, "com_before.csv")
dim(df)
dim(kma_com)
# ------------------------------------
# OR SIMPLY LOAD FROM RData Objects
load("kma_com.RData")
# ------------------------------------------------------------------------------
# ==============================================================================
# ------------------------- SPECIES ACCUMULATION CURVES ------------------------
# ==============================================================================
# Function to run species accumulation
calculate_species_curve <- function(df) {
  spec_curve <- specaccum(comm = df, method = "random", permutations = 1000)
  data.frame(
    samples = spec_curve$sites,
    richness = spec_curve$richness,
    sd = spec_curve$sd
  )
}

# Specaccum for KMA results
spec_curve <- calculate_species_curve(kma_com)

# Identify the approximate points where the plots start to flatten
flatten_point_kma <- 16

# Plot for KMA data using ggplot2
p1 <- ggplot(spec_curve, aes(x = samples, y = richness)) +
  geom_line(color = "darkblue", size = 1) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), fill = "purple", alpha = 0.2) +
  geom_vline(xintercept = flatten_point_kma, color = "orange", linetype = "dashed") +
  geom_hline(yintercept = 3792, color = "orange", linetype = "dashed") +
  annotate("text", x = flatten_point_kma, y = max(spec_curve$richness) * 0.95, 
           label = "Potential stabilization point\nafter 16 samples", hjust = 1.1, color = "steelblue") +
  theme_minimal() +
  labs(title = "Mapped to BacATG database using KMA",
       x = "Sampling Effort",
       y = "Estimated Richness") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = spec_curve$samples)
# ------------------------------------------------------------------------------
# ==============================================================================
#               OBSERVED RICHNESS, CHAO1 ESTIMATES AND ALPHA ESTIMATES
# ==============================================================================
dim(kma_com)
cat("Observed number of taxa for KMA hits:",dim(kma_com)[2],  "\n")
kma_chao1 <- fossil::chao1(kma_com, taxa.row=FALSE)
cat("Chao1 diversity for KMA hits:", round(kma_chao1), "\n")
# ------------------------------------------------------------------------------
# Alpha estimates for KMA hits
load("/Users/swati/Desktop/RData Objects/kma_com.RData")
KMA_alphas <- compute_MoM_alpha(kma_com)
KMA_alpha <- round(median(KMA_alphas$Alphas), 2)
KMA_alpha # 7.43

load("/Users/swati/Desktop/RData Objects/kma_after_var.RData")
KMA_alphas_after_var <- compute_MoM_alpha(kma_after_var)
KMA_alpha_after_var <- round(median(KMA_alphas_after_var$Alphas), 2)
KMA_alpha_after_var # 2.93

# Make dataframes for plotting
KMA_alphas_df <- data.frame(Alphas = KMA_alphas$Alphas)

# Scatter plots for all alphas
alpha_plot_a <- ggplot(KMA_alphas_df, aes(x = 1:nrow(KMA_alphas_df), y = Alphas)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  theme_minimal() +
  labs(title = "Estimated alpha values for all taxa based on KMA hits",
       x = "Taxa categories",
       y = "Alpha values") +
  geom_hline(aes(yintercept = median(Alphas)), color = "red", linetype = "dashed") +
  annotate("text", x = nrow(KMA_alphas_df) * 0.9, y = median(KMA_alphas_df$Alphas), 
           label = paste("Min =", round(min(KMA_alphas_df$Alphas), 2), "\n",
                         "Median =", round(median(KMA_alphas_df$Alphas), 2), "\n",
                         "Max =", round(max(KMA_alphas_df$Alphas), 2)), 
           hjust = 1, vjust = -8, color = "red")

# Saved to Images as "alpha_plotKMA.pdf"

KMA_after_var_alphas_df <- data.frame(Alphas = KMA_alphas_after_var$Alphas)

ggplot(KMA_after_var_alphas_df, aes(x = 1:nrow(KMA_after_var_alphas_df), y = Alphas)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  theme_minimal() +
  labs(title = "Estimated alpha values for all taxa based on KMA hits without variance towers",
       x = "Taxa categories",
       y = "Alpha values") +
  geom_hline(aes(yintercept = median(Alphas)), color = "red", linetype = "dashed") +
  annotate("text", x = nrow(KMA_after_var_alphas_df) * 0.9, 
           y = median(KMA_after_var_alphas_df$Alphas), 
           label = paste("Min =", round(min(KMA_after_var_alphas_df$Alphas), 2), "\n",
                         "Median =", round(median(KMA_after_var_alphas_df$Alphas), 2), "\n",
                         "Max =", round(max(KMA_after_var_alphas_df$Alphas), 2)), 
           hjust = 1, vjust = -5, color = "red")

# Saved to Images as "KMA_alpha_after_var_plot.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#                         ACOMP BARPLOTS
# ==============================================================================
kma_com_acomp <- acomp(kma_com, total=1)
kma_after_var_acomp <- acomp(kma_after_var, total=1)
color_palette_kma <- colorspace::sequential_hcl(ncol(kma_com_acomp), "Blues", rev = TRUE)


# Create the barplots
# Adjust layout to add space for the legend
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(4, 4))

# # Acomp bar plot for KMA mapped hits
par(mar = c(5, 3, 2, 1))  # Set margins for the plots
barplot.acomp(kma_com_acomp, 
              legend.text = FALSE, 
              plotMissings = FALSE,
              missingColor = "white",
              col = color_palette_kma, 
              main = "KMA Mapped Hits, Median Alpha=7.43",
              xlab = "Samples", ylab = "Proportion")

# Saved to Images as "sewage_barplots.pdf" were saved for both KMA and Kraken2 together
# ------------------------------------------------------------------------------
# Acomp bar plot for KMA mapped hits without variance towers
color_palette_kma_after_var <- colorspace::sequential_hcl(ncol(kma_after_var_acomp), "Blues", rev = TRUE)
par(mar = c(5, 3, 2, 1))  # Set margins for the plots
barplot.acomp(kma_after_var_acomp, 
              legend.text = FALSE, 
              plotMissings = FALSE,
              missingColor = "white",
              col = color_palette_kma_after_var, 
              main = "\nKMA Mapped Hits Without Variance Towers, Median Alpha=2.93\n",
              xlab = "Samples", ylab = "Proportion")

# Saved to Images as "sewage_barplots_after_var.pdf"

# Add a gradient legend
# Create an empty plot to display only the legend
plot.new()  # Create a new empty plot
par(mar = c(5, 1, 5, 10))  # Set margins for the plots
# Add a gradient legend
image.plot(zlim = c(0, 1), legend.only = TRUE, col = color_palette_kma, 
           legend.shrink = 0.75, horizontal = TRUE, 
           legend.args = list(text = "Relative Abundances of Taxa From Low to High", 
                              side = 1, line = 2.5, cex = 1))
# ------------------------------------------------------------------------------
# ==============================================================================
#                        CLR-TRANSFORMATION
#         BOXPLOTS FOR SPECIES DISTRIBUTION ACROSS SAMPLES
# ==============================================================================
kma_clr_com <- clr(kma_com)

# Function for orfered boxplots
create_ordered_boxplots <- function(data, title) {
  # Calculate species-wise mean CLR
  species_means <- as.data.frame(data) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "Means")
  
  # Calculate species-wise variance in CLR values
  species_vars <- as.data.frame(data) %>%
    summarise(across(everything(), var)) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "Variance")
  
  
  # Order species by descending mean CLR abundance
  ordered_species_means <- species_means %>%
    arrange(Means) %>%
    pull(Species)
  
  ordered_species_vars <- species_vars %>%
    arrange(desc(Variance)) %>%
    pull(Species)
  
  # Convert the data to long format and reorder species
  plot_data_long <- as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "CLR_Abundance") %>%
    mutate(Species = factor(Species, levels = ordered_species_vars))
  
  # Boxplot using ggplot2
  ggplot(plot_data_long, aes(x = Species, y = CLR_Abundance)) +
    geom_boxplot(aes(group = Species), fill = "darkblue",
                 color = "steelblue", outlier.color = "red") +
    labs(title = title,
         subtitle = "Boxplots for CLR-transformed Abundances\nArranged by descending species-wise variances across samples",
         x = "Species",
         y = "CLR Abundance") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          legend.position = "none")
}


kma_clr_boxplot <- create_ordered_boxplots(kma_clr_com, 
                                           "KMA Mapped hits: Taxa-wise variation across samples")

# Print boxplot
kma_clr_boxplot
# ==============================================================================
#   CLUSTERING BASED ON EUCLIDEAN DISTANCE FOR SAMPLE HOMOGENEITY TESTING
#                     FOR KMA Counts With Variance Towers
# ==============================================================================
# Creating a distance matrix for KMA_com
load("kma_com.RData")
df <- clr(kma_com)
dist_matrix <- vegan::vegdist(df, method = "euclidean")
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
hc <- hclust(dist_matrix, method = "complete")
clusters <- cutree(hc, k = 2)

plot(hc, labels = rownames(df), 
     main = "Clustering of samples based on distance matrix", xlab = "", sub = "")

# PCoA plot with ggplot2
library(ggdendro)
library(ggalt)

pcoa_data <- as.data.frame(pcoa_result$points)
pcoa_data$Sample <- rownames(pcoa_data)
pcoa_data$Clusters <- as.factor(clusters)

# Create the plot with improved circles around clusters
ggplot(pcoa_data, aes(x = V1, y = V2, label = Sample, color = Clusters)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, color = "black") +
  geom_encircle(aes(group = Clusters), stat = "identity", color = "darkblue", 
                size = 2, expand = 0.04, linetype = "dotted") +
  theme_minimal() +
  labs(title="PCoA for testing samples' homogeneity using KMA counts",
       subtitle="Clustering of samples is based on complete linkage with Euclidean distances") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  scale_color_manual(values = c("blue", "orange")) +
  xlim(min(pcoa_data$V1) - 6, max(pcoa_data$V1) + 2) +
  ylim(min(pcoa_data$V1) - 4, max(pcoa_data$V1) + 2) # Adjust the x-axis limits

# PCoA plot saved as "sewage_KMA_PCoA.pdf"
# ------------------------------------------------------------------------------
# DISTRIBUTION OF PAIRWISE EUCLIDEAN DISTANCES FOR EACH SAMPLE

# Create a long format data frame for all pairwise distances
dist_matrix_long <- as.data.frame(as.table(as.matrix(dist_matrix)))
colnames(dist_matrix_long) <- c("Sample1", "Sample2", "Distance")
dist_matrix_long <- dist_matrix_long[dist_matrix_long$Sample1 != dist_matrix_long$Sample2, ]

# Violin plot
ggplot(dist_matrix_long, aes(x = Sample1, y = Distance)) +
  geom_violin(alpha=0.5, fill = "darkblue") +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Pairwise Distances for Each Sample",
       subtitle="Based on KMA hits",
       x = "Sample",
       y = "Pairwise Distance") +
  theme(axis.text.x = element_text(size=11, angle = 0, hjust = 1, vjust = 4),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(face="italic")
  )

# Saved to Images as "violinPlotsKMA.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#   CLUSTERING BASED ON EUCLIDEAN DISTANCE FOR SAMPLE HOMOGENEITY TESTING
#                     FOR KMA Counts Without Variance Towers
# ==============================================================================
# Creating a distance matrix for KMA_com after removing variance towers
df <- clr(kma_after_var)
dist_matrix <- vegan::vegdist(df, method = "euclidean")
# ------------------------------------------------------------------------------
# Applying PCoA on distance matrix
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
hc <- hclust(dist_matrix, method = "complete")
clusters <- cutree(hc, k = 2)

plot(hc, labels = rownames(df), 
     main = "Clustering of samples based on distance matrix", xlab = "", sub = "")

# PCoA plot with ggplot2
library(ggdendro)
library(ggalt)

pcoa_data <- as.data.frame(pcoa_result$points)
pcoa_data$Sample <- rownames(pcoa_data)
pcoa_data$Clusters <- as.factor(clusters)

# Create the plot with improved circles around clusters
ggplot(pcoa_data, aes(x = V1, y = V2, label = Sample, color = Clusters)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, color = "black") +
  geom_encircle(aes(group = Clusters), stat = "identity", color = "darkblue", 
                size = 2, expand = 0.04, linetype = "dotted") +
  theme_minimal() +
  labs(title="PCoA for testing samples' homogeneity using KMA counts without variance towers",
       subtitle="Clustering of samples is based on complete linkage with Euclidean distances") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  scale_color_manual(values = c("blue", "orange")) +
  xlim(min(pcoa_data$V1) - 6, max(pcoa_data$V1) + 2) +
  ylim(min(pcoa_data$V1) - 4, max(pcoa_data$V1) + 2) # Adjust the x-axis limits

# PCoA plot saved as "sewage_KMA_After_Var_PCoA.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================



