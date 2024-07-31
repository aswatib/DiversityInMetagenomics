#!/bin/bash

# Directory containing your reports
REPORT_DIR="./individual_reports/"

# Bracken database directory
DB_DIR="/home/people/s220868/project/databases/Bacteria_Kraken"

# Output directory
OUTPUT_DIR="./bracken_temp_out"

# Final combined output file
FINAL_OUTPUT="./combined_bracken.report"


# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"
mkdir -p "(dirname "$FINAL_OUTPUT")"

# Loop through each report file in the report directory
for report in "$REPORT_DIR"/*.out; do

    # Extract the base name of the report file
    base_name=$(basename "$report" .out)

    # Define the output file name
    output_file="$OUTPUT_DIR/${base_name}_bracken.report"

    # Run the bracken command
    bracken -d "$DB_DIR" -i "$report" -o "$output_file" -l S

done


# Combine all individual Bracken output files into one
cat "$OUTPUT_DIR"/*.report > "$FINAL_OUTPUT"

# Clean up temporary files
# rm -r "$OUTPUT_DIR"
