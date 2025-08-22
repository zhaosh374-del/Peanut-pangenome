#!/bin/bash

# Define input and output paths
input_dir="/path/gff3_clean_fin"
output_dir="/path/chr_gene_id"

# Create output directory
mkdir -p "$output_dir"

# Process each sample file
for input_file in "$input_dir"/*.gff3; do
    # Get sample filename 
    sample_name=$(basename "$input_file" .gff3)
    output_file="$output_dir/$sample_name.txt"
    
    # Extract chromosome, MP ID and Target from mRNA lines
    awk '
    BEGIN {OFS="\t"}
    $3 == "mRNA" {
        # Extract chromosome (first column)
        chr = $1;
        
        # Extract MP ID
        split($9, id_parts, /ID=|;/);
        mp_id = id_parts[2];
        
        # Extract Target 
        split($9, target_parts, /Target=| /);
        target = target_parts[2];
        
        # Output three columns: chromosome, MP ID, Target
        print chr, mp_id, target;
    }' "$input_file" > "$output_file"
done