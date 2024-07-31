#!/bin/bash

# Base directory where the raw read folders are located
BASE_DIR="/home/people/s220868/project/synthetic_data/read_generation/"

# Function to process sample files in a given folder
process_samples() {
  local folder_name=$1

  # Full path to the folder
  local full_path="${BASE_DIR}${folder_name}"

  # Loop over all files in the sample folder
  for file in "${full_path}"/*.fastq; do
    # Extract the filename
    filename=$(basename "$file")
    # Remove the extension from the filename
    filename_no_ext="${filename%.fastq}"

    output_folder="./bacOriginal.map/$folder_name"
    mkdir -p "$output_folder"
    # Run KMA with the appropriate input and output filenames
    echo "RUNNING KMA WITH BacOriginal DATABASE ON: $filename_no_ext"
    kma -i "$file" -o "$output_folder/${filename_no_ext}" -t_db "/home/people/s220868/project/synthetic_data/kma_runs/bacOriginal.db/bacOriginal" -ts 0 -oa -ef -nf -nc
  done
}

# Remove any existing directories
rm -rf bacOriginal.db
rm -rf bacOriginal.map

# Make fresh directories for kma
mkdir bacOriginal.db
mkdir bacOriginal.map

# Index the genome files from Entrez (already done once so no need to repeat unless the original genome files are changed)
kma index -i /home/people/s220868/project/synthetic_data/read_generation/bacOriginal.genomes/*.fna -o bacOriginal.db/bacOriginal

# Process samples in each sample folder
process_samples "Dirichlet_alpha10"
process_samples "Dirichlet_alpha1"
process_samples "sparseDirichlet"
process_samples "sparseDirichlet2"
