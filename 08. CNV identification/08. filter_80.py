#!/usr/bin/env python3
import os
from Bio import SeqIO
from collections import defaultdict

# Define paths
protein_fasta = "/ref/arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa"
input_dir = "/path/blastp_chr_same_uniq"
output_dir = "/path/blastp_chr_same_uniq_80_80"

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def load_protein_lengths(fasta_file):
    """Read protein FASTA file and return {protein_id: length} dictionary"""
    protein_lengths = {}
    with open(fasta_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            protein_id = record.id  # Keep full ID
            protein_lengths[protein_id] = len(record.seq)
    return protein_lengths

def process_blastp_file(input_path, output_path, protein_lengths):
    """Process a single BLASTP result file"""
    # Store all records
    all_records = []
    # Store Tifrunner ID to records mapping
    tifrunner_to_records = defaultdict(list)
    
    with open(input_path, "r") as infile:
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 12:
                continue
                
            query_id = parts[0]
            subject_id = parts[1]
            identity = float(parts[2])
            alignment_length = int(parts[3])
            
            record = {
                'line': line,
                'query_id': query_id,
                'subject_id': subject_id,
                'identity': identity,
                'alignment_length': alignment_length,
                'protein_length': protein_lengths.get(subject_id, 0)
            }
            
            all_records.append(record)
            tifrunner_to_records[subject_id].append(record)
    
    # Identify Tifrunner IDs with multiple MP ID mappings
    multi_mapped = {tid for tid, recs in tifrunner_to_records.items() if len(recs) > 1}
    
    # Filter records
    filtered_records = []
    for record in all_records:
        subject_id = record['subject_id']
        
        # If 1:1 mapping, keep directly
        if subject_id not in multi_mapped:
            filtered_records.append(record)
            continue
            
        # If many-to-one mapping, apply filtering conditions
        identity = record['identity']
        alignment_length = record['alignment_length']
        protein_length = record['protein_length']
        
        # Condition 1: identity > 80
        if identity <= 80:
            continue
            
        # Condition 2: alignment_length > 80% of protein length
        if protein_length > 0 and alignment_length >= 0.8 * protein_length:
            filtered_records.append(record)
    
    # Write results
    with open(output_path, "w") as outfile:
        for record in filtered_records:
            outfile.write(record['line'])

def main():
    # Step 1: Load protein length information
    print("Loading protein lengths...")
    protein_lengths = load_protein_lengths(protein_fasta)
    print(f"Loaded lengths for {len(protein_lengths)} proteins")
    
    # Step 2: Process each BLASTP result file
    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            
            print(f"Processing {filename}...")
            process_blastp_file(input_path, output_path, protein_lengths)
            print(f"Finished processing {filename}")
    
    print("All files processed successfully!")

if __name__ == "__main__":
    main()
