# ==============================================================================
# --------- LOADING REQUIRED LIBRARIES AND AESTHETICS --------------------------
# ==============================================================================
source("PACKAGES.R")
source("CLEANER.R")
source("ALPHA_ESTIMATORS.R")
source("ENTROPY_ESTIMATORS.R")
# ==============================================================================
#                  DATA LOAD, CLEANING AND PROCESSING
# ==============================================================================
# ------------------------------------------------------------------------------
# -------------------------------- KRAKEN2  DATA  ------------------------------
# ------------------------------------------------------------------------------
dir <- "/Users/swati/Desktop/Master Thesis/Programming/abundance_tables/sewage/kraken_abundances"
setwd(dir)
combined_kraken2_output <- read.csv("combined_kraken2_output.csv", 
                                    header=TRUE, stringsAsFactors=TRUE)

# Data wrangling to convert Kraken2 results into an abundance table
combined_kraken2_output <- combined_kraken2_output[-1] # Remove first column of index
colnames(combined_kraken2_output) # Check remaining column names
combined_kraken2_output <- t(combined_kraken2_output) # Transpose the df
colnames(combined_kraken2_output) <- combined_kraken2_output[1,] # Convert colnames to species names
combined_kraken2_output <- as.data.frame(combined_kraken2_output[-1,]) # Convert to df, remove first row
class(combined_kraken2_output)

# Converting chars to numeric
combined_kraken2_output <- combined_kraken2_output %>%
  mutate(across(everything(), ~as.numeric(gsub("\\s+", "", .))))
str(combined_kraken2_output)

# Check dimensions
dim(combined_kraken2_output)
# View(combined_kraken2_output)

# Extract the sample numbers from the row names
sample_numbers <- unique(gsub(".*_(\\d+)$", "\\1", rownames(combined_kraken2_output)))
# sample_numbers <- as.numeric(unique(gsub(".*_(\\d+)$", "\\1", rownames(combined_kraken2_output))))

# Initialize a new data frame to store the combined results
combined_kraken2_results <- data.frame(matrix(ncol = ncol(combined_kraken2_output), 
                                              nrow = length(sample_numbers)))
rownames(combined_kraken2_results) <- sample_numbers
colnames(combined_kraken2_results) <- colnames(combined_kraken2_output)

# Combine the paired and singleton observations for each sample
for (sample in sample_numbers) {
  paired_rows <- grepl(paste0("report_paired_", sample), rownames(combined_kraken2_output))
  singleton_rows <- grepl(paste0("report_singletons_", sample), rownames(combined_kraken2_output))
  
  combined_kraken2_results[sample, ] <- (combined_kraken2_output[paired_rows,]) + 
    (combined_kraken2_output[singleton_rows,])
}

# Sort the observations by row numbers
combined_kraken2_results <- combined_kraken2_results[order(as.numeric(row.names(combined_kraken2_results))), ]

# Check the combined results
dim(combined_kraken2_results)

# Save the Kraken2 abundance set as RData object
kraken2_com <- combined_kraken2_results
dim(kraken2_com)
save(kraken2_com, file = "kraken2_com.RData")
# ------------------------------------
# OR SIMPLY LOAD FROM RData Objects
load("kraken2_com.RData")
# ==============================================================================
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

# Specaccum for Kraken2 results
Kraken2_speccurve <- calculate_species_curve(kraken2_com)

# Identify the approximate points where the plots start to flatten
flatten_point_kraken2 <- 17

# Plot for Kraken2 data using ggplot2
p2 <- ggplot(Kraken2_speccurve, aes(x = samples, y = richness)) +
  geom_line(color = "purple", size = 1) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), fill = "pink", alpha = 0.2) +
  geom_vline(xintercept = flatten_point_kraken2, color = "orange", linetype = "dashed") +
  geom_hline(yintercept = 9909, color = "orange", linetype = "dashed") +
  annotate("text", x = flatten_point_kraken2, y = max(Kraken2_speccurve$richness) * 0.95, 
           label = "Potential stabilization point\nbut still not quite there", hjust = 1.1, color = "steelblue") +
  theme_minimal() +
  labs(title = "\nMapped to NCBI's Bacterial and Plasmid Libraries using Kraken2",
       x = "Sampling Effort",
       y = "Estimated Richness") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks = spec_curve$samples)

# Combine the plots
grid.arrange(p1, p2, ncol = 1) # Combined with plot "p1" from SEWAGE_KMA.R
# ------------------------------------------------------------------------------
# ==============================================================================
#               OBSERVED RICHNESS, CHAO1 ESTIMATES AND ALPHA ESTIMATES
# ==============================================================================
dim(kraken2_com)
cat("Observed number of taxa for Kraken2 hits:",dim(kraken2_com)[2], "\n")
kraken2_chao1 <- fossil::chao1(kraken2_com, taxa.row=FALSE)
cat("Chao1 diversity for Kraken2 hits:", round(kraken2_chao1), "\n")
# ------------------------------------------------------------------------------
# Alpha estimates for Kraken hits
load("kraken2_com.RData")
Kraken2_alphas <- compute_MoM_alpha(kraken2_com)
Kraken2_alpha <- round(median(Kraken2_alphas$Alphas), 2)
Kraken2_alpha # 15.88

# Make dataframe for plotting
Kraken2_alphas_df <- data.frame(Alphas = Kraken2_alphas$Alphas)

# Scatter plot for Kraken2_alphas
alpha_plot_b <- ggplot(Kraken2_alphas_df, aes(x = 1:nrow(Kraken2_alphas_df), y = Alphas)) +
  geom_point(alpha = 0.5, color = "darkblue") +
  theme_minimal() +
  labs(title = "Estimated alpha values for all taxa based on Kraken2 hits",
       x = "Taxa categories",
       y = "Alpha values") +
  geom_hline(aes(yintercept = median(Alphas)), color = "red", linetype = "dashed") +
  annotate("text", x = nrow(Kraken2_alphas_df) * 0.9, y = median(Kraken2_alphas_df$Alphas), 
           label = paste("Min =", round(min(Kraken2_alphas_df$Alphas), 2), "\n",
                         "Median =", round(median(Kraken2_alphas_df$Alphas), 2), "\n",
                         "Max =", round(max(Kraken2_alphas_df$Alphas), 2)), 
           hjust = 1, vjust = -8, color = "red")

alpha_plot_b

alpha_plot_a / alpha_plot_b # Combined with alpha plots for KMA from SEWAGE_KMA.R

# Saved to Images as "alpha_plotKraken.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#                         ACOMP BARPLOTS
# ==============================================================================
# Acomp bar plot for Kraken2 mapped hits
kraken2_com_acomp <- acomp(kraken2_com, total=1)
color_palette_kraken2 <- colorspace::sequential_hcl(ncol(kraken2_com_acomp), "Blues", rev = TRUE)

par(mar = c(5, 3, 2, 1))  # Set margins for the plots
barplot.acomp(kraken2_com_acomp, 
              legend.text = FALSE, 
              col = color_palette_kraken2, 
              plotMissings = FALSE,
              missingColor = "white",
              main = "Kraken2 Mapped Hits, Median Alpha=15.88",
              xlab = "Samples", ylab = "Proportion")

# Saved to Images as "sewage_barplots.pdf" were saved for both KMA and Kraken2 together

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
kraken2_clr_com <- clr(kraken2_com)

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

# Creating boxplots
kraken2_clr_boxplot <- create_ordered_boxplots(kraken2_clr_com, 
                                               "Kraken2 Mapped hits: Taxa-wise variation across samples")

# Print boxplot
kraken2_clr_boxplot
# ------------------------------------------------------------------------------
# ==============================================================================
#   CLUSTERING BASED ON EUCLIDEAN DISTANCE FOR SAMPLE HOMOGENEITY TESTING
#                               FOR KRAKEN2 Counts
# ==============================================================================
# Creating a distance matrix for Kraken2_com
load("kraken2_com.RData")
df <- clr(kraken2_com)
dist_matrix <- vegan::vegdist(df, method = "euclidean")
pcoa_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
hc <- hclust(dist_matrix, method = "complete")
clusters <- cutree(hc, k = 2)
# ------------------------------------------------------------------------------
# Create a dendrogram object from the hclust result
plot(hc, labels = rownames(df), 
     main = "Clustering of samples based on distance matrix", xlab = "", sub = "")
# ------------------------------------------------------------------------------
# PCoA plot with ggplot2
pcoa_data <- as.data.frame(pcoa_result$points)
pcoa_data$Sample <- rownames(pcoa_data)
pcoa_data$Clusters <- as.factor(clusters)

# Create the plot with improved circles around clusters
ggplot(pcoa_data, aes(x = V1, y = V2, label = Sample, color = Clusters)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, color = "black") +
  geom_encircle(aes(group = Clusters), stat = "identity", color = "darkblue", 
                size = 2, expand = 0.05, linetype = "dotted") +
  theme_minimal() +
  labs(title="PCoA for testing samples' homogeneity using Kraken2 counts",
       subtitle="Clustering of samples is based on complete linkage with Euclidean distances") +
  xlab("PCoA1") +
  ylab("PCoA2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  scale_color_manual(values = c("blue", "orange")) +
  xlim(min(pcoa_data$V1) - 2, max(pcoa_data$V1) + 2) +
  ylim(min(pcoa_data$V1) - 2, max(pcoa_data$V1) + 2) # Adjust the x-axis limits

# PCoA plot saved as "sewage_Kraken2_PCoA.pdf"
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
       subtitle="Based on Kraken2 hits",
       x = "Sample",
       y = "Pairwise Distance") +
  theme(axis.text.x = element_text(size=11, angle = 0, hjust = 1, vjust = 4),
        #plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.title = element_blank(),
        plot.subtitle = element_text(face="italic")
  )

# Saved to Images as "violinPlotsKraken2.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
