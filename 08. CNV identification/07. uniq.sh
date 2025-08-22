#!/bin/bash

# Define input and output paths
input_dir="/path/blastp_chr_same"
output_dir="/path/blastp_chr_same_uniq"

# Create output directory
mkdir -p "$output_dir"

# Process each sample file
for input_file in "$input_dir"/*.txt; do
    # Get sample filename
    sample_name=$(basename "$input_file")
    output_file="$output_dir/$sample_name"
    
    # Remove duplicates and sort by first column
    awk '!seen[$1]++' "$input_file" | sort -k1,1 > "$output_file"
done