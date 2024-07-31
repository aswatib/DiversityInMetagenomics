# ==============================================================================
#               SNR ANALYSIS AND PLOTS FOR KMA AND KRAKEN HITS
# ==============================================================================
# Function to calculate SNR for synthetic abundance counts
calculate_snr <- function(df) {
  
  # Perform CLR transformation
  clr_data <- as.data.frame(clr(df))
  
  # Estimate the signal (mean) and noise (standard deviation)
  signals <- colMeans(clr_data)
  noises <- colSds(as.matrix(clr_data))
  
  # Calculate SNR
  snr <- signals / noises
  
  # Collect SNR results as a data frame
  snr_results <- data.frame(
    "Species" = colnames(clr_data),
    "Signal" = signals,
    "Noise" = noises,
    "SNR" = snr)
  
  # Calculate tolerable noise
  total_noise <- snr_results$Noise
  
  # Evaluate conditions separately
  is_snr_below_0 <- snr_results$SNR < 0
  is_noise_below_median <- snr_results$Noise < median(total_noise)
  snr_results$Color <- ifelse(!is_snr_below_0, "darkblue",
                              ifelse(is_snr_below_0 & !is_noise_below_median, "grey50", "orange"))
  
  # Category column for the legend
  snr_results$Category <- factor(snr_results$Color, levels = c("darkblue", "grey50", "orange"),
                                 labels = c("High-quality signal and low noise (major part of the composition)",
                                            "Low-quality signal and high noise (potentially spurious counts)",
                                            "Low-quality signal but relatively low noise (potentially rare taxa)"))
  
  return(snr_results)
}

# Function to plot SNR results
plot_snr <- function(snr_results, title, subtitle) {
  ggplot(snr_results, aes(x = Signal, y = SNR, color=Category)) +
    geom_point() + 
    labs(title = title,
         subtitle = subtitle,
         x = "Signal",
         y = "SNR",
         color = "Taxa categories indicated by solid colored points:") +
    theme_minimal() +
    scale_color_manual(values=c("darkblue", "grey50", "orange"))
}
# ------------------------------------------------------------------------------
# Load KMA and KRAKEN2 abundance datasets
load("../RData Objects/kma_com.RData")
load("../RData Objects/kraken2_com.RData")

# Apply SNR on KMA and KRAKEN2 abundance datasets
kma_snr <- calculate_snr(kma_com)
kraken2_snr <- calculate_snr(kraken2_com)

# Plot SNR
a <- plot_snr(kma_snr, "SNR Analysis for Sewage Data", "For KMA mapped hits")
b <- plot_snr(kraken2_snr, "", "For Kraken2 mapped hits")

# Combine plots and apply common theme settings with a single legend
combined_plot <- a + b + plot_layout(guides = "collect") &
  theme(plot.title = element_text(size=11, vjust=1.5, face="bold"),
        plot.subtitle = element_text(size=10, face="italic"),
        legend.position = "bottom",
        legend.direction = "vertical")

# Print the combined plot for KMA and Kraken2 SNR
print(combined_plot)

# Saved to Images as "snr.pdf"
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SNR based filtering for KMA com
# Identifying species that have higher than median noise
kma_species_to_filter <- kma_snr[kma_snr$Color != "grey50", ]
kma_species_to_filter <- kma_species_to_filter$Species

# Filter out these species from the original DataFrame
kma_filtered_df <- kma_com[, colnames(kma_com) %in% kma_species_to_filter] # Taxa that are orange and darkblue
dim(kma_filtered_df) # Taxa that are orange and darkblue
dim(kma_com)
# ------------------------------------------------------------------------------
# SNR based filtering for Kraken2 com
# Identifying species that have higher than median noise
kraken2_species_to_filter <- kraken2_snr[kraken2_snr$Color != "grey50", ] # Taxa that are orange and darkblue
kraken2_species_to_filter <- kraken2_species_to_filter$Species # Taxa that are orange and darkblue

# Filter out these species from the original DataFrame
kraken2_filtered_df <- kraken2_com[, colnames(kraken2_com) %in% kraken2_species_to_filter]

dim(kraken2_filtered_df) # Taxa that are orange and darkblue
dim(kraken2_com)

kraken2_filtered_df

# List of rare taxa
kma_rare_taxa <- kma_snr[kma_snr$Color == "orange", ]
kraken2_rare_taxa <- kraken2_snr[kraken2_snr$Color == "orange", ]

dim(kma_rare_taxa)
dim(kraken2_rare_taxa)

common_species <- intersect(kma_rare_taxa$Species, kraken2_rare_taxa$Species)
length(common_species)

dim(kraken2_com)
# ==============================================================================
#     ENTROPY AND EFFECTIVE RICHNESS AFTER REMOVING NOISE BASED ON SNR
# ==============================================================================
# Estimate KMA entropies after removing noise and save as RData object
kma_after_noise_entropies <- estimate_individual_entropies(kma_filtered_df)
save(kma_after_noise_entropies, file="../RData Objects/kma_after_noise_entropies.RData")

# Load again when required
load("../RData Objects/kma_after_noise_entropies.RData")
# ------------------------------------------------------------------------------
# Kraken2 entropies after removing noise and save as RData object
kraken2_after_noise_entropies <- estimate_individual_entropies(kraken2_filtered_df)
save(kraken2_after_noise_entropies, file="../RData Objects/kraken2_after_noise_entropies.RData")

# Load again when required
load("../RData Objects/kraken2_after_noise_entropies.RData")
# ------------------------------------------------------------------------------
# ==============================================================================