# Entropy estimate function for bootstrapping
estimate_entropy_boot <- function(data, indices) {
  data <- data[indices, , drop = FALSE]
  H_piga <- piga_entropy(data)
  return(H_piga)
}

# Tweaked entropy estimate function for bootstrapping
estimate_entropy_boot_Singles <- function(data, indices) {
  data <- data[indices, , drop = FALSE]
  H_piga <- apply(data, 1, function(row) piga_entropy_single(as.numeric(row)))
  return(H_piga)
}

# ------------------------------------------------------------------------------
# Set up parallel processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
# ------------------------------------------------------------------------------
# BOOTSTRAPPING FOR KMA COUNT AFTER REMOVING VARIANCE TOWERS
dim(kma_com)
dim(kma_after_var)

set.seed(123)  # for reproducibility
kma_after_var_data_matrix <- as.matrix(kma_after_var)
kma_after_var_bootstrap_resultsCombined <- boot(kma_after_var_data_matrix, 
                                                     statistic = estimate_entropy_boot, 
                                                     R = 50, stype="i")

save(kma_after_var_bootstrap_resultsCombined, 
     file="../RData Objects/kma_after_var_bootstrap_resultsCombined.RData")

load("../RData Objects/kma_after_var_bootstrap_resultsCombined.RData")

# Calculate the mean and 95% confidence intervals for each entropy method
bootstrap_means <- colMeans(kma_after_var_bootstrap_resultsCombined$t)
bootstrap_ci <- apply(kma_after_var_bootstrap_resultsCombined$t, 2, function(x) quantile(x, c(0.025, 0.975)))
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
set.seed(123)  # for reproducibility
kma_after_var_bootstrap_resultsSingles <- boot(kma_after_var_data_matrix, 
                                                    statistic = estimate_entropy_boot_Singles, 
                                                    R = 50, stype="i")

save(kma_after_var_bootstrap_resultsSingles, 
     file="../RData Objects/kma_after_var_bootstrap_resultsSingles.RData")

load("../RData Objects/kma_after_var_bootstrap_resultsSingles.RData")

# Calculate the mean and 95% confidence intervals for each entropy method
bootstrap_means <- colMeans(kma_after_var_bootstrap_resultsSingles$t)
bootstrap_ci <- apply(kma_after_var_bootstrap_resultsSingles$t, 2, function(x) quantile(x, c(0.025, 0.975)))
# ------------------------------------------------------------------------------
# BEST ESTIMATE OF SPECIES RICHNESS (BASED ON KMA AFTER VAR) AND PIGA's METHOD
# BASED ON SINGLES (INDIVIDUAL SAMPLE ENTROPIES AFTER BOOTSTRAPPING)

mean(bootstrap_means)
round(mean(exp(bootstrap_means)))

# BEST ESTIMATES OF CIs
bootstrap_ci_mat <- matrix(bootstrap_ci, ncol = 18, byrow = TRUE)

# Extract the lower and upper bounds
lower_bounds <- apply(bootstrap_ci, 2, function(x) x[1])
upper_bounds <- apply(bootstrap_ci, 2, function(x) x[2])

# Calculate the mean of the lower and upper bounds
mean_lower_bound <- mean(lower_bounds)
mean_upper_bound <- mean(upper_bounds)

# Print the mean CI values
cat("Mean 2.5% CI:", mean_lower_bound, "\n")
cat("Mean 97.5% CI:", mean_upper_bound, "\n")

# CI for effective richness
exp(mean_lower_bound)
exp(mean_upper_bound)

# Create a density plot for all bootstrapped values
density_values <- density(kma_after_var_bootstrap_resultsSingles$t)

# Plot with customizations
plot(density_values,
     main = "Density of bootstrapped entropy estimates",
     sub="Higher peaks represent entropy estimates that are more common among bootstrapped samples",
     xlab = "Bootstrapped entropy values",
     ylab = "Density",
     col = "darkblue",
     lwd = 2,
     cex.lab = 1,
     cex.axis = 1,
     bty = "n")

# Fill the density plot
polygon(density_values, col = "darkblue", border = "darkblue")

# Density_values is the density object created from your data
density_data <- data.frame(x = density_values$x, y = density_values$y)

# Create the density plot using ggplot2
ggplot(density_data, aes(x = x, y = y)) +
  geom_area(fill = "darkblue", color = "darkblue") +
  geom_vline(xintercept = 6.531496, color = "skyblue", linetype = "dashed", size = 0.5) +
  annotate("text", x = 6.5, y = max(density_data$y) * 0.9, label = "Mean Entropy = 6.531496", color = "darkblue", angle = 0, lwd=2) +
  labs(title = "Density of bootstrapped entropy estimates",
       subtitle = "Higher peaks represent entropy estimates that are more common among bootstrapped samples",
       x = "Bootstrapped entropy values",
       y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Saved to Images as "PDFBootstrappedEntropy.pdf"
# ------------------------------------------------------------------------------
# Check for NaNs or extreme values in bootstrap results
summary(kma_after_var_bootstrap_resultsSingles$t)
hist(kma_after_var_bootstrap_resultsSingles$t,
     main="Distribution of bootstrapped entropy estimes",
     xlab="Bootstrapped values of entropy (Piga's method)")
# ------------------------------------------------------------------------------
# Stop the clusters after bootstrapping
stopCluster(cl)
registerDoSEQ()
