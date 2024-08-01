# ==============================================================================
#                ANALYSIS OF MAPPED COUNTS FOR SYNTHETIC READS
# ==============================================================================
#                            SOURCE SCRIPTS
# ==============================================================================
source("PACKAGES.R")
source("CLEANER.R")
source("ALPHA_ESTIMATORS.R")
source("ENTROPY_ESTIMATORS.R")
source("SYNTHETIC_UNMAPPED.R")
# ==============================================================================
#                      WORKING DIRECTORY AND DATA (MAPSTATs)
# ==============================================================================
dir <- "/Users/swati/Desktop/Master Thesis/Programming/synthetic_mapstats/"
setwd(dir)

# Load BacOriginal mapstat files
BacOriginal_Dirichlet_alpha10_mapstat <- read_multiple_mapstats("synthetic_bacOriginal_mapstats/Dirichlet_alpha10/")
BacOriginal_Dirichlet_alpha1_mapstat <- read_multiple_mapstats("synthetic_bacOriginal_mapstats/Dirichlet_alpha1/")
BacOriginal_sparseDirichlet_mapstat <- read_multiple_mapstats("synthetic_bacOriginal_mapstats/sparseDirichlet/")
BacOriginal_sparseDirichlet2_mapstat <- read_multiple_mapstats("synthetic_bacOriginal_mapstats/sparseDirichlet2/")

BacATG_Dirichlet_alpha1_mapstat <- read_multiple_mapstats("synthetic_bacdb_mapstats/Dirichlet_alpha1/")
BacATG_Dirichlet_alpha10_mapstat <- read_multiple_mapstats("synthetic_bacdb_mapstats/Dirichlet_alpha10/")
BacATG_sparseDirichlet_mapstat <- read_multiple_mapstats("synthetic_bacdb_mapstats/sparseDirichlet/")
BacATG_sparseDirichlet2_mapstat <- read_multiple_mapstats("synthetic_bacdb_mapstats/sparseDirichlet2/")
# ------------------------------------------------------------------------------
# ==============================================================================
#        VISUALIZING OUTLIERS IN DEPTH VARIANCE OF MAPSTAT RESULTS
# ==============================================================================
# mapstat_df <- BacATG_Dirichlet_alpha1_mapstat
mapstat_df <- BacATG_sparseDirichlet_mapstat

# TESTING FILTERING USING OUTLIERS
Q1 <- quantile(mapstat_df$depthVariance, 0.25)
Q3 <- quantile(mapstat_df$depthVariance, 0.75)
IQR <- Q3 - Q1
upper_threshold <- Q3 + 1.5 * IQR # for outliers

# Create a logical vector for outliers
is_outlier <- mapstat_df$depthVariance > upper_threshold
table(is_outlier)

# Plot the data
plot(mapstat_df$depthVariance, 
     main = "Depth Variance in KMA's Mapstat Results",
     sub = "Based on one of the mapped reads results, with a total of 63 outlier values",
     xlab = "Template or reference sequences (refSequence column in mapstat)",
     ylab = "Depth Variance (depthVariance column in mapstat)",
     col = ifelse(is_outlier, "red", "black"),
     pch = ifelse(is_outlier, 16, 1))

# Add a legend
legend("topright", legend = c("Outlier", "Normal"), 
       col = c("red", "black"), pch = c(16, 1))

# Saved to Images as "scatterplotVarianceTowers.pdf"

# Boxplot
ggplot(mapstat_df, aes(x = refSequence, y = depthVariance)) +
  geom_boxplot(fill="grey50") +
  geom_point(data = mapstat_df[is_outlier, ], 
             aes(x = refSequence, y = depthVariance), color = "red") +
  labs(title = "Boxplot of Depth Variance in KMA's Mapstat Results",
       subtitle = "Extreme values in depthVariance column have been termed as Variance towers",
       x = "Reference Sequences",
       y = "Depth Variance") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_blank(),  # Hide x-axis labels to avoid clutter
    axis.text.y = element_text(size = 10),
    axis.ticks.x = element_line(size = 0.5),
    axis.ticks.y = element_line(size = 0.5)
  )

# Saved to Images as "boxplotVarianceTowers.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#                          RESHAPE THE DATA
# ==============================================================================
# Function to process mapstats before and after removing the variance towers
process_dataset <- function(mapstat_df) {
  df <- mapstat_df %>%
    dplyr::select(sample_id, refSequence, fragmentCount) %>%
    pivot_wider(names_from = refSequence, values_from = fragmentCount, values_fill = 0)
  
  reshaped_df1 <- data.frame(df)
  before_var <- clean_and_order(as.data.frame(reshaped_df1)) # clean_and_order function from CLEANER.R
  
  # Calculate IQR and upper threshold for filtering
  Q1 <- quantile(mapstat_df$depthVariance, 0.25)
  Q3 <- quantile(mapstat_df$depthVariance, 0.75)
  IQR <- Q3 - Q1
  upper_threshold <- Q3 + 1.5 * IQR
  
  # Filter out outliers
  filtered_df <- mapstat_df %>%
    filter(depthVariance <= upper_threshold)
  
  reshaped_df2 <- filtered_df %>%
    dplyr::select(sample_id, refSequence, fragmentCount) %>%
    pivot_wider(names_from = refSequence, values_from = fragmentCount, values_fill = 0)
  
  reshaped_df2 <- data.frame(reshaped_df2)
  after_var <- clean_and_order(as.data.frame(reshaped_df2))
  
  list(before_var = before_var, after_var = after_var)
}

# Define dataset names and corresponding data frames
BacOriginal_datasets <- c("BacOriginal_Dirichlet_alpha10_mapstat", "BacOriginal_Dirichlet_alpha1_mapstat",
                          "BacOriginal_sparseDirichlet_mapstat", "BacOriginal_sparseDirichlet2_mapstat")

# Loop over each dataset and plot manually -- for BacOriginal
for (i in 1:length(BacOriginal_datasets)) {
  df <- get(BacOriginal_datasets[i])
  
  # Process mapstat to make abundance tables
  processed <- process_dataset(df)
  
  # Remove the "_mapstat" suffix from the dataset name
  base_name <- sub("_mapstat$", "", BacOriginal_datasets[i])
  
  # Save the before variance dataframe with the base name
  assign(base_name, processed$before_var, envir = .GlobalEnv)
  
  # Save the after variance dataframe with "_after_var" suffix
  assign(paste0(base_name, "_after_var"), processed$after_var, envir = .GlobalEnv)
}

# Define dataset names and corresponding data frames
BacATG_datasets <- c("BacATG_Dirichlet_alpha10_mapstat", "BacATG_Dirichlet_alpha1_mapstat", 
                     "BacATG_sparseDirichlet_mapstat", "BacATG_sparseDirichlet2_mapstat")


# Loop over each dataset and plot manually -- for BacATG
for (i in 1:length(BacATG_datasets)) {
  df <- get(BacATG_datasets[i])
  
  # Process mapstat to make abundance tables
  processed <- process_dataset(df)
  
  # Remove the "_mapstat" suffix from the dataset name
  base_name <- sub("_mapstat$", "", BacATG_datasets[i])
  
  # Save the before variance dataframe with the base name
  assign(base_name, processed$before_var, envir = .GlobalEnv)
  
  # Save the after variance dataframe with "_after_var" suffix
  assign(paste0(base_name, "_after_var"), processed$after_var, envir = .GlobalEnv)
}
# ------------------------------------------------------------------------------
# ==============================================================================
#           OBSERVED ABUNDANCE BASED ON MAPPING RESULTS
# ==============================================================================
# Define dataset names and corresponding data frames
df_names <- c("BacATG_Dirichlet_alpha10", "BacATG_Dirichlet_alpha1", 
              "BacATG_sparseDirichlet", "BacATG_sparseDirichlet2",
              "BacOriginal_Dirichlet_alpha10", "BacOriginal_Dirichlet_alpha1",
              "BacOriginal_sparseDirichlet", "BacOriginal_sparseDirichlet2")

# Initialize a list to store the dimensions
rm(df_dimensions)
df_dimensions <- list()

# Function to refine names
refine_name <- function(name) {
  name <- gsub("_", " ", name) # Replace underscores with spaces
  name <- tools::toTitleCase(name) # Convert to title case
  return(name)
}

# Loop through each data frame name, get the data frame, store its dimensions, and refine the name
for (df_name in df_names) {
  df <- get(df_name)
  refined_name <- refine_name(df_name)
  df_dimensions[[refined_name]] <- dim(df)[2]
}

# Create a data frame to hold the results
OSR <- data.frame(
  Observed_Species_Richness = sapply(df_dimensions, function(x) paste(x, collapse = " x "))
)

kable(OSR) %>% kable_styling()
# ------------------------------------------------------------------------------
# ==============================================================================
#                    SPECIES ACCUMULATION
# ==============================================================================
# Function for species accumulation
calculate_species_curve <- function(df) {
  spec_curve <- specaccum(comm = df, method = "random", permutations = 1000)
  data.frame(
    samples = spec_curve$sites,
    richness = spec_curve$richness,
    sd = spec_curve$sd
  )
}
# ------------------------------------------------------------------------------
#     SPECIES ACCUMULATION FOR KMA RESULTS USING BacOriginaldb -- SYNTHETIC DATA
# ------------------------------------------------------------------------------
# Define dataset names and corresponding data frames
BacOriginal_datasets <- c("BacOriginal_Dirichlet_alpha10", "BacOriginal_Dirichlet_alpha1",
                          "BacOriginal_sparseDirichlet", "BacOriginal_sparseDirichlet2")

BacOriginal_labels <- c("Set A\nDirichlet, targeted alpha=10", 
                        "Set B\nDirichlet, targeted alpha=1", 
                        "Set C\nSparse Dirichlet, targeted alpha=0.5", 
                        "Set D\nSparse Dirichlet, targeted alpha=0.01")

BacOriginal_y_limits <- list(c(38, 42), c(36, 42), c(31, 42), c(0, 42))

# Set up the plotting layout
par(mfrow = c(2, 2))

# Loop over each dataset and plot manually
for (i in 1:length(BacOriginal_datasets)) {
  df <- get(BacOriginal_datasets[i])
  df <- df[-1]
  spec_curve <- calculate_species_curve(df)
  
  # Plot the base plot
  plot(spec_curve$samples, spec_curve$richness, type="l", col = "darkblue", 
       main = BacOriginal_labels[i], 
       xlab = "Sampling Effort", 
       ylab = "Estimated Richness", 
       ylim = BacOriginal_y_limits[[i]],
       lwd=2,
       xaxt = "n")
  
  # Adding custom x-axis
  axis(1, at = seq(1, max(spec_curve$samples), by = 1))  # Customize the x-axis
  
  # Adding the standard deviation as a shaded area
  polygon(c(spec_curve$samples, rev(spec_curve$samples)),
          c(spec_curve$richness + spec_curve$sd, 
            rev(spec_curve$richness - spec_curve$sd)),
          col = adjustcolor("purple", alpha.f = 0.1), border = NA)
  
  # Adding horizontal abline to point out original species richness
  abline(h = 40, lty = 2, col = "red")
}
# Saved to images as "specaccum_BacOriginal.pdf"
# ------------------------------------------------------------------------------
#     SPECIES ACCUMULATION FOR KMA RESULTS USING BacATG -- SYNTHETIC DATA
# ------------------------------------------------------------------------------
# Define dataset names and corresponding data frames
BacATG_datasets <- c("BacATG_Dirichlet_alpha10", "BacATG_Dirichlet_alpha1", 
                     "BacATG_sparseDirichlet", "BacATG_sparseDirichlet2")

BacATG_labels <- c("Set A\nDirichlet, targeted alpha=10", "Set B\nDirichlet, targeted alpha=1", 
                   "Set C\nSparse Dirichlet, targeted alpha=0.5", 
                   "Set D\nSparse Dirichlet, targeted alpha=0.01")

# y_limits <- list(c(10, 80), c(20, 80), c(0, 60), c(5, 150))
BacATG_y_limits <- list(c(30, 70), c(30, 70), c(30, 70), c(0, 50))

# Set up the plotting layout
par(mfrow = c(2, 2))

# Loop over each dataset and plot manually
for (i in 1:length(BacATG_datasets)) {
  df <- get(BacATG_datasets[[i]])
  df <- df[-1]
  spec_curve <- calculate_species_curve(df)
  
  # Plot the base plot
  plot(spec_curve$samples, spec_curve$richness, type="l", col = "darkblue", 
       main = BacATG_labels[i], 
       xlab = "Sampling Effort", 
       ylab = "Estimated Richness", 
       ylim = BacATG_y_limits[[i]],
       lwd=2,
       xaxt = "n")
  
  # Adding custom x-axis
  axis(1, at = seq(1, max(spec_curve$samples), by = 1))  # Customize the x-axis
  
  # Adding the standard deviation as a shaded area
  polygon(c(spec_curve$samples, rev(spec_curve$samples)),
          c(spec_curve$richness + spec_curve$sd, 
            rev(spec_curve$richness - spec_curve$sd)),
          col = adjustcolor("purple", alpha.f = 0.1), border = NA)
  
  # Adding horizontal abline to point out original species richness
  abline(h = 40, lty = 2, col = "red")
}
# Saved to images as "specaccum_BacATG.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#                    MAPPING AND ABUNDANCE HEATMAPS
#                Heatmaps using normalized counts (proportions)
# ==============================================================================
# HEATMAPS FOR MAPPED HITS FOR BOTH BacOriginal AND BacATG DATABASES

# Create a list containing the names of abundance counts data frames
df_names <- c("BacATG_Dirichlet_alpha10", "BacATG_Dirichlet_alpha1", 
              "BacATG_sparseDirichlet", "BacATG_sparseDirichlet2",
              "BacOriginal_Dirichlet_alpha10", "BacOriginal_Dirichlet_alpha1",
              "BacOriginal_sparseDirichlet", "BacOriginal_sparseDirichlet2")

# Loop through each data frame and apply modifications
for (df_name in df_names) {
  # Get the data frame by its name
  df <- get(df_name)
  
  # Wrangle the dataframe
  df <- df[, -1] # Remove the sample ID column
  
  # Normalize the counts to get proportions
  df_freq <- apply(df, 2, function(x) {x / sum(x, na.rm = TRUE)})
  
  # Apply CLR -- Agreed upon with CB that normalized counts are fine for these heatmaps
  # df_freq <- clr(df)
  
  # Create long dataframe using reshape2::melt
  df_melted <- setNames(reshape2::melt(as.matrix(df_freq), na.rm = FALSE), 
                        c("Sample", "Species", "Value"))
  
  # Assign column data types
  df_melted <- df_melted %>%
    mutate(Sample = as.factor(Sample),
           Species = as.factor(Species),
           Value = as.numeric(Value))
  
  # Assign the modified long dataframe back to its original name
  assign(paste0(df_name, "_melted"), df_melted, envir = .GlobalEnv)
}

# PLOTTING
# Determine the overall minimum and maximum values for both datasets
overall_min <- 0
overall_max <- 1

par(mfrow=c(2,2))

# BacOriginal_db <- BacOriginal_Dirichlet_alpha1_melted
# BacATG_db <- BacATG_Dirichlet_alpha1_melted

BacOriginal_db <- BacOriginal_Dirichlet_alpha10_melted
BacATG_db <- BacATG_Dirichlet_alpha10_melted

# PLOTTING
# Determine the overall minimum and maximum values for both datasets
overall_min <- min(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)
overall_max <- max(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)

BacOriginal_heatmap <- ggplot(BacOriginal_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.title = element_text(size=14, 
                                  vjust=2),
        plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x=element_text(vjust=-1),
        legend.position = "none") +
  labs(title= "Heatmaps for the mapped reads' proportions in Set A (targeted alpha=10)",
       subtitle = "Mapping database: BacOriginal", x = "Samples", y = "Species", fill = "Frequency")

BacATG_heatmap <- ggplot(BacATG_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        text = element_text(size = 12), 
        axis.title.x = element_text(vjust=-1, size=12),
        legend.position = "none") +
  labs(title="",
       subtitle = "Mapping database: BacATG", x = "Samples", y = NULL, fill = "Frequency")

# BacOriginal_heatmap 
# BacATG_heatmap

combined_plots <- grid.arrange(BacOriginal_heatmap, BacATG_heatmap, nrow=1, ncol = 2)
combined_plots

# Saved to images as "combined_heatmap_Dirichlet10.pdf"
# Size 12 by 8
# ------------------------------------------------------------------------------
BacOriginal_db <- BacOriginal_Dirichlet_alpha1_melted
BacATG_db <- BacATG_Dirichlet_alpha1_melted

# PLOTTING
# Determine the overall minimum and maximum values for both datasets
overall_min <- min(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)
overall_max <- max(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)

BacOriginal_heatmap <- ggplot(BacOriginal_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.title = element_text(size=14, 
                                  vjust=2),
        plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x=element_text(vjust=-1),
        legend.position = "none") +
  labs(title= "Heatmaps for the mapped reads' proportions in Set B (targeted alpha=1)",
       subtitle = "Mapping database: BacOriginal", x = "Samples", y = "Species", fill = "Frequency")

BacATG_heatmap <- ggplot(BacATG_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x = element_text(vjust=-1, size=12),
        legend.position = "none") +
  labs(title="",
       subtitle = "Mapping database: BacATG", x = "Samples", y = NULL, fill = "Frequency")

# BacOriginal_heatmap 
# BacATG_heatmap

combined_plots <- grid.arrange(BacOriginal_heatmap, BacATG_heatmap, nrow=1, ncol = 2)
combined_plots

# Saved to images as "combined_heatmap_Dirichlet1.pdf"
# Size 12 by 8
# ------------------------------------------------------------------------------
BacOriginal_db <- BacOriginal_sparseDirichlet_melted
BacATG_db <- BacATG_sparseDirichlet_melted

# PLOTTING
# Determine the overall minimum and maximum values for both datasets
overall_min <- min(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)
overall_max <- max(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)

BacOriginal_heatmap <- ggplot(BacOriginal_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.title = element_text(size=14, 
                                  vjust=2),
        plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x=element_text(vjust=-1),
        legend.position = "none") +
  labs(title= "Heatmaps for the mapped reads' proportions in Set C (targeted alpha=0.5)",
       subtitle = "Mapping database: BacOriginal", x = "Samples", y = "Species", fill = "Frequency")

BacATG_heatmap <- ggplot(BacATG_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x = element_text(vjust=-1, size=12),
        legend.position = "none") +
  labs(title="",
       subtitle = "Mapping database: BacATG", x = "Samples", y = NULL, fill = "Frequency")

# BacOriginal_heatmap 
# BacATG_heatmap

combined_plots <- grid.arrange(BacOriginal_heatmap, BacATG_heatmap, nrow=1, ncol = 2)
combined_plots

# Saved to images as "combined_heatmap_sparseDirichlet.pdf"
# Size 12 by 8
# ------------------------------------------------------------------------------
BacOriginal_db <- BacOriginal_sparseDirichlet2_melted
BacATG_db <- BacATG_sparseDirichlet2_melted

# PLOTTING
# Determine the overall minimum and maximum values for both datasets
overall_min <- min(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)
overall_max <- max(c(BacOriginal_db$Value, BacATG_db$Value), na.rm = TRUE)

BacOriginal_heatmap <- ggplot(BacOriginal_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.title = element_text(size=14, 
                                  vjust=2),
        plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x=element_text(vjust=-1),
        legend.position = "none") +
  labs(title= "Heatmaps for the mapped reads' proportions in Set D (targeted alpha=0.01)",
       subtitle = "Mapping database: BacOriginal", x = "Samples", y = "Species", fill = "Frequency")

BacATG_heatmap <- ggplot(BacATG_db, aes(x = Sample, y = Species, fill = Value)) +
  geom_tile(position="identity", color = "white") +
  scale_fill_gradient2(low = "grey", mid= "lightgrey", high = "darkblue", limits = c(overall_min, overall_max)) +
  theme_minimal() +
  theme(plot.subtitle = element_text(size=11),
        axis.text.x = element_text(angle=0, hjust=0.75), 
        text = element_text(size = 12), 
        axis.title.x = element_text(vjust=-1, size=12),
        legend.position = "none") +
  labs(title="",
       subtitle = "Mapping database: BacATG", x = "Samples", y = NULL, fill = "Frequency")

# BacOriginal_heatmap 
# BacATG_heatmap

combined_plots <- grid.arrange(BacOriginal_heatmap, BacATG_heatmap, nrow=1, ncol = 2)
combined_plots

# Saved to images as "combined_heatmap_sparseDirichlet2.pdf"
# Size 12 by 8
# ------------------------------------------------------------------------------
# ==============================================================================
#           CLR-TRANSFORMATION OF ABUNDANCES FOR ALL 8 SETS
# ==============================================================================
# Apply CLR on the abundance data for of the 8 set
clr_BacOriginal_Dirichlet_alpha10 <- clr(BacOriginal_Dirichlet_alpha10[-1])
clr_BacOriginal_Dirichlet_alpha1 <- clr(BacOriginal_Dirichlet_alpha1[-1])
clr_BacOriginal_sparseDirichlet <- clr(BacOriginal_sparseDirichlet[-1])
clr_BacOriginal_sparseDirichlet2 <- clr(BacOriginal_sparseDirichlet2[-1])

clr_BacATG_Dirichlet_alpha10 <- clr(BacATG_Dirichlet_alpha10[-1])
clr_BacATG_Dirichlet_alpha1 <- clr(BacATG_Dirichlet_alpha1[-1])
clr_BacATG_sparseDirichlet <- clr(BacATG_sparseDirichlet[-1])
clr_BacATG_sparseDirichlet2 <- clr(BacATG_sparseDirichlet2[-1])
# ------------------------------------------------------------------------------
# ==============================================================================
#    BOXPLOTS FOR CLR-TRANSFORMED ABUNDANCES FOR EACH SPECIES IN EACH SET
#     TO IDENTIFY SAMPLES THAT HAVE SPECIES-SPECIFIC OUTLIERS AND VISUALISE 
#                         SPECIES-WISE DISTRIBUTION
# ==============================================================================
clr_dfs <- c("clr_BacOriginal_Dirichlet_alpha1",
             "clr_BacOriginal_Dirichlet_alpha10",
             "clr_BacOriginal_sparseDirichlet",
             "clr_BacOriginal_sparseDirichlet2",
             "clr_BacATG_Dirichlet_alpha1",
             "clr_BacATG_Dirichlet_alpha10",
             "clr_BacATG_sparseDirichlet",
             "clr_BacATG_sparseDirichlet2")

create_ordered_boxplots <- function(data, title, subtitle) {
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
         subtitle = subtitle,
         x = "Species",
         y = "CLR Abundance") +
    theme_minimal() +
    theme(plot.title=element_text(face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
          legend.position = "none")
}

# Create ordered plots for BacOriginal CLR abundances
BacOriginal_clr_boxplot1 <- create_ordered_boxplots(clr_BacOriginal_Dirichlet_alpha10, 
                                                    "Distribution of CLR-Transformed Abundances After Mapping", 
                                                    "Arranged by descending species-wise variances across samples\n\nSet A Dirichlet, targeted alpha=10\n\nBased on BacOriginal hits\n")

BacOriginal_clr_boxplot2 <- create_ordered_boxplots(clr_BacOriginal_Dirichlet_alpha1, 
                                                    "Distribution of CLR-Transformed Abundances After Mapping",
                                                    "Arranged by descending species-wise variances across samples\n\nSet B Dirichlet, targeted alpha=1\n\nBased on BacOriginal hits\n")

BacOriginal_clr_boxplot3 <- create_ordered_boxplots(clr_BacOriginal_sparseDirichlet, 
                                                    title="Distribution of CLR-Transformed Abundances After Mapping", 
                                                    subtitle="Arranged by descending species-wise variances across samples\n\nSet C Sparse Dirichlet, targeted alpha=0.5\n\nBased on BacOriginal hits\n")
BacOriginal_clr_boxplot4 <- create_ordered_boxplots(clr_BacOriginal_sparseDirichlet2, 
                                                    "Distribution of CLR-Transformed Abundances After Mapping", 
                                                    "Set D Sparse Dirichlet, targeted alpha=0.01\n\nBased on BacOriginal hits\n")

# Create ordered plots for BacATG CLR abundances
BacATG_clr_boxplot1 <- create_ordered_boxplots(clr_BacATG_Dirichlet_alpha10, "", "\nBased on BacATG hits\n")
BacATG_clr_boxplot2 <- create_ordered_boxplots(clr_BacATG_Dirichlet_alpha1, "", "\nBased on BacATG hits\n")
BacATG_clr_boxplot3 <- create_ordered_boxplots(clr_BacATG_sparseDirichlet, "", "\nBased on BacATG hits\n")
BacATG_clr_boxplot4 <- create_ordered_boxplots(clr_BacATG_sparseDirichlet2, "", "\nBased on BacATG hits\n")

# Print all the boxplots
BacOriginal_clr_boxplot3 # <---- Saved to Images as "synthetic_clr_boxplots1.pdf"
BacATG_clr_boxplot3 # <---- Saved to Images as "synthetic_clr_boxplots2.pdf"

BacOriginal_clr_boxplot3 / BacATG_clr_boxplot3 # <---- Saved to Images as "synthetic_clr_boxplots.pdf"
BacOriginal_clr_boxplot1 / BacATG_clr_boxplot1 # <---- Saved to Images as "synthetic_clr_boxplots1.pdf"
BacOriginal_clr_boxplot2 / BacATG_clr_boxplot2 # <---- Saved to Images as "synthetic_clr_boxplots2.pdf"
BacOriginal_clr_boxplot4 / BacATG_clr_boxplot4 # <---- Saved to Images as "synthetic_clr_boxplots4.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#           CHAO1 ESTIMATES BEFORE REMOVING VARIANCE TOWERS
# ==============================================================================
# Define dataset names and corresponding data frames
datasets <- c("BacOriginal_Dirichlet_alpha10", "BacOriginal_Dirichlet_alpha1",
              "BacOriginal_sparseDirichlet", "BacOriginal_sparseDirichlet2",
              "BacATG_Dirichlet_alpha10", "BacATG_Dirichlet_alpha1", 
              "BacATG_sparseDirichlet", "BacATG_sparseDirichlet2")

# Calculate chao1 for BacOriginal sets
for (i in 1:length(datasets)) {
  dataset_name <- datasets [i]
  df <- get(datasets [i])
  df <- df[-1]
  H_chao1 <- fossil::chao1(df, taxa.row=FALSE)
  print(paste("Chao1 diversity for:", dataset_name, ":", round(H_chao1)))
}
# ------------------------------------------------------------------------------
# ==============================================================================
#           ALPHA ESTIMATES BEFORE REMOVING VARIANCE TOWERS
# ==============================================================================
# Re-estimation of alpha values for mapped counts using compute_MoM_alpha from ALPHA_ESTIMATORS.R
for (i in 1:length(datasets)) {
  dataset_name <- datasets[i]
  df <- get(datasets[i])
  df <- df[-1]
  
  alpha <- compute_MoM_alpha(df)
  print(paste("Mean MoM Alpha for:", dataset_name, ":", round(mean(alpha$Alphas), 2)))
}
# 
# [1] "Mean MoM Alpha for: BacOriginal_Dirichlet_alpha10 : 8.1"
# [1] "Mean MoM Alpha for: BacOriginal_Dirichlet_alpha1 : 1.2"
# [1] "Mean MoM Alpha for: BacOriginal_sparseDirichlet : 0.69"
# [1] "Mean MoM Alpha for: BacOriginal_sparseDirichlet2 : 0.04"
# [1] "Mean MoM Alpha for: BacATG_Dirichlet_alpha10 : 4.49"
# [1] "Mean MoM Alpha for: BacATG_Dirichlet_alpha1 : 0.69"
# [1] "Mean MoM Alpha for: BacATG_sparseDirichlet : 0.43"
# [1] "Mean MoM Alpha for: BacATG_sparseDirichlet2 : 0.04"
# ------------------------------------------------------------------------------
# ==============================================================================
#           CHAO1 ESTIMATES AFTER REMOVING VARIANCE TOWERS
# ==============================================================================
datasets_after_var <- c("BacATG_Dirichlet_alpha10_after_var", 
                        "BacATG_Dirichlet_alpha1_after_var", 
                        "BacATG_sparseDirichlet_after_var", 
                        "BacATG_sparseDirichlet2_after_var",
                        "BacOriginal_Dirichlet_alpha10_after_var", 
                        "BacOriginal_Dirichlet_alpha1_after_var",
                        "BacOriginal_sparseDirichlet_after_var", 
                        "BacOriginal_sparseDirichlet2_after_var")

# Calculate chao1 for all sets
for (i in 1:length(datasets_after_var)) {
  dataset_name <- datasets_after_var[i]
  df <- get(datasets_after_var[i])
  df <- df[-1]
  H_chao1 <- fossil::chao1(df, taxa.row=FALSE)
  print(paste("Chao1 diversity for:", dataset_name, ":", round(H_chao1)))
}
# ------------------------------------------------------------------------------
# ==============================================================================
#           ALPHA ESTIMATES AFTER REMOVING VARIANCE TOWERS
# ==============================================================================
# Calculate alpha for all sets
for (i in 1:length(datasets_after_var)) {
  dataset_name <- datasets_after_var[i]
  df <- get(datasets_after_var[i])
  df <- df[-1]
  
  alpha <- compute_MoM_alpha(df)
  print(paste("Mean MoM Alpha for:", dataset_name, ":", round(mean(alpha$Alphas), 2)))
}
# ------------------------------------------------------------------------------
# ==============================================================================
#         BacATG ALIGNMENT COUNTS -- FOR COMPARISON WITH MAPPED COUNTS
# ==============================================================================
# Function to process mapstats before and after removing the variance towers for alignment counts
process_dataset_aln <- function(mapstat_df) {
  reshaped_df <- mapstat_df %>%
    dplyr::select(sample_id, refSequence, fragmentCountAln) %>%
    pivot_wider(names_from = refSequence, values_from = fragmentCountAln, values_fill = 0)
  
  reshaped_df <- data.frame(reshaped_df)
  before_var <- clean_and_order(as.data.frame(reshaped_df))
  
  median_depthVar <- median(mapstat_df$depthVariance)
  high_threshold <- median_depthVar * 1
  
  varbelow10k_df <- mapstat_df %>%
    filter(depthVariance <= high_threshold)
  
  reshaped_df <- varbelow10k_df %>%
    dplyr::select(sample_id, refSequence, fragmentCountAln) %>%
    pivot_wider(names_from = refSequence, values_from = fragmentCountAln, values_fill = 0)
  
  reshaped_df <- data.frame(reshaped_df)
  after_var <- clean_and_order(as.data.frame(reshaped_df))
  
  list(before_var = before_var, after_var = after_var)
}

# Define dataset names and corresponding data frames
datasets <- c("BacOriginal_Dirichlet_alpha1_mapstat", "BacOriginal_Dirichlet_alpha10_mapstat", 
              "BacOriginal_sparseDirichlet_mapstat", "BacOriginal_sparseDirichlet2_mapstat",
              "BacATG_Dirichlet_alpha1_mapstat", "BacATG_Dirichlet_alpha10_mapstat", 
              "BacATG_sparseDirichlet_mapstat", "BacATG_sparseDirichlet2_mapstat")


# Loop over each dataset and plot manually -- for BacOriginal
for (i in 1:length(datasets)) {
  df <- get(datasets[i])
  
  # Process mapstat to make abundance tables
  processed <- process_dataset_aln(df)
  
  # Remove the "_mapstat" suffix from the dataset name
  base_name <- sub("_mapstat$", "", datasets[i])
  
  # Save the before variance dataframe with the base name
  assign(paste0(base_name, "_before_var_aln"), processed$before_var, envir = .GlobalEnv)
  
  # Save the after variance dataframe with "_after_var" suffix
  assign(paste0(base_name, "_after_var_aln"), processed$after_var, envir = .GlobalEnv)
}

BacOriginal_Dirichlet_alpha1_before_var_aln
BacOriginal_Dirichlet_alpha1_after_var_aln
# ==============================================================================