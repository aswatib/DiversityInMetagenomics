# LIBRARIES AND PACKAGES 
# ----------------------- Custom Packages from Github --------------------------
# install.packages("devtools")
# devtools::install_github("roeder/mapstatHelpers")
# --------------------------- CRAN Packages ------------------------------------
# List of required packages
required_packages <- c(
  "rmarkdown",
  "pagedown",
  "vegan", 
  "boot",
  "knitr", 
  "readr", 
  "kableExtra", 
  "tidyr", 
  "dplyr", 
  "purrr",
  "ggplot2",
  "reshape2", 
  "gridExtra", 
  "cowplot", 
  "viridis", 
  "viridisLite", 
  "RColorBrewer", 
  "patchwork", 
  "gsl",
  "compositions",
  "DiagrammeR",
  "pracma",
  "torch",
  "pscl",
  "rentrez",
  "networkD3",
  "car",
  "stats",
  "data.table",
  "mapstatHelpers",
  "entropy",
  "MASS",
  "matrixStats",
  "fossil",
  "data.table"
)


# Check if each package is installed, install if not and then load
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ------------------------ BiocManager Packages --------------------------------
# Ensure that BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of required Bioconductor packages
required_bioc_packages <- c("DESeq2", "DirichletMultinomial", "HilbertCurve")

# Check if each Bioconductor package is installed, install if not, and then load
for (pkg in required_bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ------------------------------- AESTHETICS -----------------------------------
colors_viridis <- viridis(37, option = "C")
set1 <- brewer.pal(9, "Blues")
set2 <- brewer.pal(9, "Greens")
set3 <- brewer.pal(9, "YlGnBu")
set4 <- brewer.pal(9, "PuOr")
pastel1 <- brewer.pal(9, "Pastel1")
colors <- c(set1, set2, set3, set4, set1, set2)
colors2 <- c(pastel1, colors_viridis, pastel1)
