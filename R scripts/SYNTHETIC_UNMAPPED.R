# ------------------------------------------------------------------------------
# ==============================================================================
#                            SOURCE SCRIPTS
# ==============================================================================
source("PACKAGES.R")
source("CLEANER.R")
source("ENTROPY_ESTIMATORS.R")
source("ALPHA_ESTIMATORS.R")
# ------------------------------------------------------------------------------
# ==============================================================================
#                      WORKING DIRECTORY AND DATA
# ==============================================================================
# Set working directory
dir <- "/Users/swati/Desktop/Master Thesis/Programming/sample_stats/"
setwd(dir)
# ==============================================================================
# Read the sample stats CSV files for all synthetic reads
Dirichlet_alpha10 <- read.csv("stats_Dirichlet_alpha10.csv", 
                              header=TRUE, sep=",", stringsAsFactors=TRUE)

Dirichlet_alpha1 <- read.csv("stats_Dirichlet_alpha1.csv", 
                             header=TRUE, sep=",", stringsAsFactors=TRUE)


sparseDirichlet <- read.csv("stats_sparseDirichlet.csv", 
                            header=TRUE, sep=",", stringsAsFactors=TRUE)


sparseDirichlet2 <- read.csv("stats_sparseDirichlet2.csv", 
                             header=TRUE, sep=",", stringsAsFactors=TRUE)

# ------------------------------------------------------------------------------
# ==============================================================================
# # ALPHA ESTIMATION -- based on method of moments from Wikipedia
# ==============================================================================
# Define the count sets and their names
count_sets <- list(Dirichlet_alpha10 = Dirichlet_alpha10, 
                   Dirichlet_alpha1 = Dirichlet_alpha1, 
                   sparseDirichlet = sparseDirichlet, 
                   sparseDirichlet2 = sparseDirichlet2)

# Initialize a data frame to store results
results_df <- data.frame(
  Count_Set = character(),
  Min_Alpha_J = numeric(),
  Max_Alpha_J = numeric(),
  Mean_Alpha_J = numeric(),
  Median_Alpha_J = numeric(),
  Alpha_0 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each count set
for (count_set_name in names(count_sets)) {
  df <- count_sets[[count_set_name]]
  
  # Convert to matrix and remove the first column (assuming it is an index or non-numeric column)
  data <- as.matrix(df[-1])
  
  # Calculate n, mean proportions, and variance proportions
  n <- rowSums(data, na.rm = TRUE)
  mean_proportions <- colMeans(data / n, na.rm = TRUE)
  variance_proportions <- apply(data / n, 2, function(x) var(x, na.rm = TRUE))
  
  # Calculate alpha_j and alpha_0
  alpha_j <- mean_proportions * (((mean_proportions * (1 - mean_proportions)) / variance_proportions) - 1)
  alpha_0 <- sum(alpha_j, na.rm = TRUE)
  
  # Store results in the data frame
  results_df <- results_df %>%
    add_row(
      Count_Set = count_set_name,
      Min_Alpha_J = min(alpha_j, na.rm = TRUE),
      Max_Alpha_J = max(alpha_j, na.rm = TRUE),
      Mean_Alpha_J = mean(alpha_j, na.rm = TRUE),
      Median_Alpha_J = median(alpha_j, na.rm = TRUE),
      Alpha_0 = alpha_0
    )
}

# Print the results in a table
print(results_df)
kable(results_df)
# ==============================================================================
# ALPHA ESTIMATION USING THE FUNCTION FROM ALPHA_ESTIMATORS.R
# ==============================================================================
# Calculate summaries for each dataset
summary_alpha10 <- compute_MoM_alpha(Dirichlet_alpha10[-1])
summary_alpha10_dirmult <- compute_alpha_dirmult(Dirichlet_alpha10[-1])


summary_alpha1 <- compute_MoM_alpha(Dirichlet_alpha1[-1])
summary_alpha1_dirmult <- compute_alpha_dirmult(Dirichlet_alpha1[-1])

summary_alpha05 <- compute_MoM_alpha(sparseDirichlet[-1])
summary_alpha05_dirmult <- compute_alpha_dirmult(sparseDirichlet[-1])

summary_alpha001 <- compute_MoM_alpha(sparseDirichlet2[-1])
summary_alpha001_dirmult <- compute_alpha_dirmult(sparseDirichlet2[-1])

# Combine all summaries into one table
summary_table <- rbind(
  "Set A" = summary_alpha10[[1]],
  "Set B" = summary_alpha1[[1]],
  "Set C" = summary_alpha05[[1]],
  "Set D" = summary_alpha001[[1]]
)

summary_table
# ------------------------------------------------------------------------------
# ==============================================================================
#             TARGETED DISTRIBUTIONS OF UNMAPPED SYNTHETIC READS
#                                 Histograms
# ==============================================================================
# Do a colSums for all species (K categories) and add as a new row named "Total"
# Stats set: Dirichlet alpha 1
Total <- (colSums(Dirichlet_alpha1[,-1]))
Total <- as.data.frame(t(Total))
normalize_Dirichlet_alpha1 <- Total[-1]
norm1 <- normalize_Dirichlet_alpha1 / sum(normalize_Dirichlet_alpha1)
norm1$Sample <- "Dirichlet_alpha1"
normalized_Dirichlet_alpha1 <- norm1 %>%
  dplyr::select(Sample, everything())

# Stats set: Dirichlet alpha 10
Total <- (colSums(Dirichlet_alpha10[,-1]))
Total <- as.data.frame(t(Total))
normalize_Dirichlet_alpha10 <- Total[-1]
norm1 <- normalize_Dirichlet_alpha10 / sum(normalize_Dirichlet_alpha10)
norm1$Sample <- "Dirichlet_alpha10"
normalized_Dirichlet_alpha10 <- norm1 %>%
  dplyr::select(Sample, everything())


# Stats set: Sparse Dirichlet
Total <- (colSums(sparseDirichlet[,-1]))
Total <- as.data.frame(t(Total))
normalize_sparseDirichlet <- Total[-1]
norm1 <- normalize_sparseDirichlet / sum(normalize_sparseDirichlet)
norm1$Sample <- "sparseDirichlet"
normalized_sparseDirichlet <- norm1 %>%
  dplyr::select(Sample, everything())

# Stats set: Sparse Dirichlet2
Total <- (colSums(sparseDirichlet2[,-1]))
Total <- as.data.frame(t(Total))
normalize_sparseDirichlet2 <- Total[-1]
norm1 <- normalize_sparseDirichlet2 / sum(normalize_sparseDirichlet2)
norm1$Sample <- "sparseDirichlet2"
normalized_sparseDirichlet2 <- norm1 %>%
  dplyr::select(Sample, everything())


# Combine the last rows into a new dataframe
counts <- rbind(normalized_Dirichlet_alpha10, normalized_Dirichlet_alpha1,
                normalized_sparseDirichlet, normalized_sparseDirichlet2)
rownames(counts) <- counts$Sample
counts <- counts[,-1]

# checking entropy and ESR
com <- counts
ShannonIndex <- vegan::diversity(com, index = "shannon")
exp(ShannonIndex)

# Plot all stats with targeted distributions
par(mfrow=c(2,2))

stat1 <- t(counts[1,])
stat2 <- t(counts[2,])
stat3 <- t(counts[3,])
stat4 <- t(counts[4,])

# Define statistic names and their corresponding titles
stat_names <- c("stat1", "stat2", "stat3", "stat4")
stat_titles <- c("Dirichlet alpha=10", 
                 "Dirichlet alpha=1", 
                 "Sparse Dirichlet alpha=0.5", 
                 "Sparse Dirichlet alpha=0.01")

stat_color <- c("steelblue", "blue", "darkblue", "darkturquoise")
x_points <- 1:34  # Creates a sequence from 1 to 34

# Iterate over each statistic
for (i in seq_along(stat_names)) {
  # Get the statistic and its title
  stat <- get(stat_names[i])
  title <- stat_titles[i]
  #  color <- stat_color[i]
  
  # Plot the statistic
  stat <- sort(stat, decreasing = TRUE)
  
  # Create the plots
  targeted_distributions <- plot(stat, 
                                 type = "h", 
                                 xlab = "Species categories (i)", 
                                 ylab = "Proportions (p_i)", 
                                 # xaxt="n",
                                 main = title, 
                                 col = "grey50",
                                 lwd=3,
                                 cex.main=1.2, cex.lab=1.2)
  # Add custom x-axis ticks
  #axis(1, at = x_points, labels = x_points)
  
  # Print the plots
  targeted_distributions
}

# Saved to Images as "targeted_distributions.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#     TARGETED DISTRIBUTIONS USING HMP SIMULATED COUNTS (NOT FROM THE READS)
#                               Histograms
# ==============================================================================
library(HMP)
set.seed(123)

stat1 <- colSums(Dirichlet.multinomial(Nrs=rep(200, 10), shape=rep(10,40)))
stat2 <- colSums(Dirichlet.multinomial(Nrs=rep(200, 10), shape=rep(1,40)))
stat3 <- colSums(Dirichlet.multinomial(Nrs=rep(200, 10), shape=rep(0.5,40)))
stat4 <- colSums(Dirichlet.multinomial(Nrs=rep(200, 10), shape=rep(0.01,40)))

# Define statistic names and their corresponding titles
stat_names <- c("stat1", "stat2", "stat3", "stat4")
stat_titles <- c("Dirichlet alpha=10", 
                 "Dirichlet alpha=1", 
                 "Sparse Dirichlet alpha=0.5", 
                 "Sparse Dirichlet alpha=0.01")

stat_color <- c("steelblue", "blue", "darkblue", "darkturquoise")
x_points <- 1:40  # Creates a sequence from 1 to 34

# Plot all stats with targeted distributions
par(mfrow=c(2,2))

# Iterate over each statistic
for (i in seq_along(stat_names)) {
  # Get the statistic and its title
  stat <- get(stat_names[i])
  title <- stat_titles[i]
  #  color <- stat_color[i]
  
  # Plot the statistic
  stat <- stat / sum(stat)
  stat <- sort(stat, decreasing = TRUE)
  
  # Create the plots
  targeted_distributions <- plot(stat,
                                 type = "h",
                                 xlab = "Species categories (i)",
                                 ylab = "Proportions (p_i)",
                                 xaxt="n",
                                 main = title,
                                 col = "grey50",
                                 lwd=3,
                                 cex.main=1.2, cex.lab=1.2)
  lines(stat, col= "darkblue", lwd=0.5)
  
  # Add custom x-axis ticks
  axis(1, at = x_points, labels = x_points)
  
  # Print the plots
  targeted_distributions
}

# Saved to Images as "targeted_distributions.pdf"
# ------------------------------------------------------------------------------
# ==============================================================================
#         TARGETED DISTRIBUTIONS OF UNMAPPED SYNTHETIC READS
# SHOWN AS SPECIES-WISE VARIANCE THROUGH BOXPLOTS WITH CLR TRANSFORMED COUNTS
# ==============================================================================
# Apply CLR after removing the sample ID column
clr_Dirichlet_alpha1 <- clr(Dirichlet_alpha1[-1])
clr_Dirichlet_alpha10 <- clr(Dirichlet_alpha10[-1])
clr_sparseDirichlet <- clr(sparseDirichlet[-1])
clr_sparseDirichlet2 <- clr(sparseDirichlet2[-1])

# Function to create a plot data frame with species reordered
create_plot_data <- function(clr_data) {
  plot_data <- as.data.frame(clr_data) %>%
    pivot_longer(cols = everything(), names_to = "Species", values_to = "CLR_Abundance")
  
  # Calculate variance or mean for each species
  species_stats <- plot_data %>%
    group_by(Species) %>%
    summarize(Variance = var(CLR_Abundance, na.rm = TRUE),
              Mean = mean(CLR_Abundance, na.rm = TRUE))
  
  # Reorder species by variance or mean
  plot_data <- plot_data %>%
    mutate(Species = factor(Species, levels = species_stats$Species[order(species_stats$Variance, decreasing = TRUE)]),
           SpeciesIndex = as.numeric(Species))
  
  return(plot_data)
}

# Create plot data for each dataset
plot_data_long10 <- create_plot_data(clr_Dirichlet_alpha10)
plot_data_long1 <- create_plot_data(clr_Dirichlet_alpha1)
plot_data_long05 <- create_plot_data(clr_sparseDirichlet)
plot_data_long001 <- create_plot_data(clr_sparseDirichlet2)

# Define plot aesthetics
my_plot_style <- list(
  scale_x_continuous(breaks = 1:length(unique(plot_data_long05$Species)),
                     labels = 1:length(unique(plot_data_long05$Species))),
  theme_minimal(), 
  theme(plot.title = element_text(face="bold"),
    plot.subtitle = element_text(size=10),
    axis.text.x = element_text(size = 10),
        axis.title.x = element_text(vjust=-2),
        legend.position = "none")
)

# Generate boxplots
a <- ggplot(plot_data_long10, aes(x = SpeciesIndex, y = CLR_Abundance)) +
  geom_boxplot(aes(group = Species), fill = "darkblue",
               color = "grey", outliers = TRUE, outlier.size = 0.5, 
               outlier.color = "blue") +
  labs(title="Distribution of normalized species counts from the raw unmapped synthetic reads",
    subtitle = "Arranged in desceding order of species-wise variances across samples\n\nTargeted alpha = 10\nMean of estimated alpha = 8.17",
       x = "Species",
       y = "Normalized counts") + my_plot_style


b <- ggplot(plot_data_long1, aes(x = SpeciesIndex, y = CLR_Abundance)) +
  geom_boxplot(aes(group = Species), fill = "steelblue",
               color = "grey", outliers = TRUE, outlier.size = 0.5, 
               outlier.color = "blue") +
  labs(subtitle = "Targeted alpha = 1\nMean of estimated alpha = 1.16",
       x = "Species",
       y = "Normalized counts") + my_plot_style

c <- ggplot(plot_data_long05, aes(x = SpeciesIndex, y = CLR_Abundance)) +
  geom_boxplot(aes(group = Species), fill = "skyblue",
               color = "grey50", outliers = TRUE, outlier.size = 0.5, 
               outlier.color = "blue") +
  labs(subtitle = "Targeted alpha = 0.5\nMean of estimated alpha = 0.69",
       x = "Species",
       y = "Normalized counts") + my_plot_style


d <- ggplot(plot_data_long001, aes(x = SpeciesIndex, y = CLR_Abundance)) +
  geom_boxplot(aes(group = Species), fill = "lightblue",
               color = "grey50", 
               outliers = TRUE, outlier.size = 0.5, 
               outlier.color = "blue") +
  labs(subtitle = "Targeted alpha = 0.01\nMean of estimated alpha = 0.039",
       x = "Species",
       y = "Normalized counts") + my_plot_style

# Arrange the plots
library(patchwork)
a / b / c / d

# Saved to Images as "targeted_distributions_boxplots.pdf"
# ==============================================================================