# Gsub options
#cleaned_names <- gsub("^X\\.|\\.$", "", columns) # Remove any leading X's
#colnames(y1) <- gsub("\\.", " ", cleaned_names) # Replace the period between genus and species by space

# Function to clean the species names in the columns except the first column of sample numbers
clean_names <- function(df) {
  samples <- gsub('([[:alpha:]]+)([[:digit:]]+)', '\\1 \\2', df[,1]) # Clean sample names
  
  # Retain the column name "Sample" for the first column
  colnames(df)[1] <- "Sample"
  
  # Order rows based on sample numbers
  df <- df[order(df$Sample), ]
  
  # Extract column names from the dataframe excluding the first column
  species <- colnames(df[-1])
  
  # Define the input string
  input_string <- as.character(species)
  
  # Replace periods with spaces
  input_string <- gsub("\\.", " ", input_string)
  
  # Remove trailing numbers or periods from species names
  input_string <- gsub("\\.\\d*?$", "", input_string)
  
  # Initialize an empty vector to store the species names
  species_names <- character(length(input_string))
  
  # Extract the second and third English words from each string
  for (i in 1:length(input_string)) {
    words <- unlist(strsplit(input_string[i], "\\s+"))
    if (length(words) >= 3) {
      species_names[i] <- paste(words[3:4], collapse = " ")
    } else {
      species_names[i] <- NA
    }
  }
  
  # Replace original column names with cleaned species names
  non_missing_indices <- !is.na(species_names)
  colnames(df)[-1][non_missing_indices] <- species_names[non_missing_indices]
  # colnames(df)[-1] <- species_names
  
  # Add back the sample column
  df$Sample <- samples
  
  # Order rows based on sample numbers
  df <- df[order(df$Sample), ]
  
  # Change row names back to just numbers without the word "sample"
  rownames(df) <- as.numeric(gsub("sample", "", rownames(df)))
  
  # Return the modified dataframe
  return(df)
}

# ------------------------------------------------------------------------------
# Function to collapse all records up to the species level
collapse_species <- function(df) {
  # Save the Sample column as a vector
  samples <- (df[[1]])
  
  # Select all columns except the first one
  df_excluding_sample <- df[, -1]
  rownames(df_excluding_sample) <- NULL
  
  # Select only numeric columns
  df_numeric <- df_excluding_sample[, sapply(df_excluding_sample, is.numeric)]
  
  # Identify unique species names
  colnames(df_numeric) <- gsub("\\d+", "", colnames(df_numeric))
  colnames(df_numeric) <- gsub("\\d+\\.*$", "", colnames(df_numeric))
  colnames(df_numeric) <- gsub("\\.$", "", colnames(df_numeric))
  
  unique_species1 <- unique(gsub("\\d+", "", colnames(df_numeric)))
  unique_species2 <- unique(gsub("\\d+\\.*$", "", unique_species1))
  unique_species <- gsub("\\.$", "", unique_species2)
  
  # Create a new dataframe to store aggregated counts
  new_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  
  # Iterate over each unique species
  for (species in unique_species) {
    # Identify columns for the current species
    species_columns_to_sum <- grep(paste0("^", species), colnames(df_numeric))
    
    # Debugging statement
    #   print(species_columns_to_sum)
    
    # Sum up the counts for the current species
    new_df[[species]] <- rowSums(df_numeric[, species_columns_to_sum, drop = FALSE], na.rm = TRUE)
    
  }
  
  # Replace periods at the end of column names with an empty string
  colnames(new_df) <- gsub("\\.$", "", colnames(new_df))
  
  return(new_df)
}

# ------------------------------------------------------------------------------
# Function to collapse all records up to the species level after the species names are cleaned
collapse_genus <- function(data, species_col = "Species") {
  genus_data <- data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    pivot_longer(-Sample, names_to = species_col, values_to = "Counts") %>%
    mutate(Genus = sapply(strsplit(get(species_col), " "), `[`, 1)) %>%
    group_by(Sample, Genus) %>%
    summarise(TotalCounts = sum(Counts), .groups = 'drop') %>%
    pivot_wider(names_from = Genus, values_from = TotalCounts, values_fill = list(TotalCounts = 0))
  
  # Convert back to a matrix if needed
  genus_matrix <- as.matrix(genus_data[-1])
  rownames(genus_matrix) <- genus_data$Sample
  
  return(genus_matrix)
}
# ------------------------------------------------------------------------------
# Function to identify columns that have family and order names instead of genus
identify_family_order <- function(species_names) {
  # Calculate the condition
  condition = grepl("ceae$", species_names) | grepl("ales$", species_names)
  
  # Loop through each species name and print it if condition is TRUE
  if (any(condition)) {
    cat("Species with family or order suffixes:\n")
    for (i in seq_along(species_names)) {
      if (condition[i]) {
        cat(species_names[i], "\n")
      }
    }
  }
}
# ------------------------------------------------------------------------------
# Function to identify family or order names and return matching species names
identify_family_order <- function(species_names) {
  # Calculate the condition
  condition = grepl("eae$", species_names) | grepl("ales$", species_names) | grepl(" bacterium$", species_names)
  
  # Collect and return species names that meet the condition
  matching_species <- species_names[condition]
  
  # Print matching species names if any
  if (length(matching_species) > 0) {
    cat("Species with family or order suffixes:\n")
    #cat(matching_species, sep="\n")
  }
  
  return(matching_species)
}
# ------------------------------------------------------------------------------
# Function to order the row names in ascending order of sample numbers
clean_and_order <- function(df) {
  
  # Apply clean_names function to clean the names of species
  df <- clean_names(df)
  
  # Apply collapse up to species level
  df <- collapse_species(df)
  
  # Set row names and change casing
  rownames(df) <- gsub("sample", "Sample", df$Sample, ignore.case=TRUE)
  
  # Extract only numeric part of the Sample column
  sample_numbers <- as.numeric(gsub("\\D", "", df$Sample))
  
  # Order rows
  df <- df[order(sample_numbers), ]
  
  # Change the sample column with ordered row names
  df$Sample <- rownames(df)
  
  # Change rownames to numbers only
  rownames(df) <- gsub("Sample", "", df$Sample, ignore.case=TRUE)
  
  # Assign the result back to the original data frame
  return(df)
}

# ------------------------------------------------------------------------------
# check_column_names <- function(df) {
#   # Extract column names
#   column_names <- colnames(df)
#   
#   # Check each column name for the presence of numbers
#   columns_with_numbers <- grep("[0-9]", column_names)
#   
#   # Print column names with numbers
#   if (length(columns_with_numbers) > 0) {
#     cat("Column names with numbers:\n")
#     for (col_index in columns_with_numbers) {
#       cat(column_names[col_index], "\n")
#     }
#   } else {
#     cat("No column names with numbers found.\n")
#   }
# }