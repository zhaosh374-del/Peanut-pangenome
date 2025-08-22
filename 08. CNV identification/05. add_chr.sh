#!/bin/bash

# Define paths
blastp_dir="/path/blastp"
chr_gene_dir="/path/chr_gene_id"
output_dir="//path/blastp_chr_gene_id"

# Create output directory
mkdir -p "$output_dir"

# Process each sample
for blastp_file in "$blastp_dir"/*.txt; do
    # Get sample name
    sample_name=$(basename "$blastp_file")
    chr_gene_file="$chr_gene_dir/$sample_name"
    output_file="$output_dir/$sample_name"
    
    # Add chromosome information
    awk -F'\t' '
    # Read chr_gene file to build mapping
    NR==FNR {
        # $2 is MP ID, $1 is chromosome
        chr_map[$2] = $1
        next
    }
    # Process blastp file
    {
        # If first column exists in mapping, add chromosome prefix
        if ($1 in chr_map) {
            $1 = chr_map[$1] "_" $1
        }
        # Print entire line (tab-separated)
        print $0
    }
    ' "$chr_gene_file" "$blastp_file" > "$output_file"
done