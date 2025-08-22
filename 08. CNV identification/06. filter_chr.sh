#!/bin/bash

# Define input and output paths
input_dir="/path/blastp_chr_gene_id"
output_dir="/path/blastp_chr_same"

# Create output directory
mkdir -p "$output_dir"

# Process each sample file
for input_file in "$input_dir"/*.txt; do
    # Get sample filename
    sample_name=$(basename "$input_file")
    output_file="$output_dir/$sample_name"
    
    # Filter rows where chromosome numbers match
    awk '
    {
        # Extract chromosome number from first column 
        split($1, col1_parts, "_");
        chr_num1 = substr(col1_parts[1], 4);  
        
        # Extract chromosome number from second column 
        split($2, col2_parts, ".");
        ah_part = col2_parts[5];  # Get Ah01g000200 part
        chr_num2 = substr(ah_part, 3, 2);  
        
        # Only keep rows where chromosome numbers match
        if (chr_num1 == chr_num2) {
            print $0;
        }
    }
    ' "$input_file" > "$output_file"
done