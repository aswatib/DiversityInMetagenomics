#!/bin/bash

# Base directory where raw read folders are located
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

    output_folder="./bacdb.map/$folder_name"
    mkdir -p "$output_folder"
    # Run KMA with the appropriate input and output filenames
    echo "RUNNING KMA WITH BACDB ATG DATABASE ON: $filename_no_ext"
    kma -i "$file" -o "$output_folder/${filename_no_ext}" -t_db "/home/people/s220868/project/databases/bacdb/bacteria.ATG" -ts 0 -oa -ef -nf -nc -tmp
  done
}

# Remove existing kma directories
rm -rf "bacdb.map"

# Make fresh directories for kma
mkdir "bacdb.map"

# Process samples in each sample folder
process_samples "Dirichlet_alpha10"
process_samples "Dirichlet_alpha1"
process_samples "sparseDirichlet"
process_samples "sparseDirichlet2"
