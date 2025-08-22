#!/usr/bin/env python3
import os
import re
from collections import defaultdict

input_dir = "/path/gff3"
output_dir = "/pathpath/gff3_clean_fin"

os.makedirs(output_dir, exist_ok=True)

def process_file(input_path, output_path):
    # Store best mRNA records by start position
    best_records_st = {}
    all_entries = []  # Store all lines in order
    gene_blocks = defaultdict(list)  # Store mRNA and child features
    
    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            # Handle comment lines
            if line.startswith('##'):
                if line.startswith('##PAF'):
                    continue
                all_entries.append(('COMMENT', line))
                continue
                
            parts = line.split('\t')
            if len(parts) < 9:
                all_entries.append(('OTHER', line))
                continue
                
            if parts[2] == 'mRNA':
                attrs = {k:v for k,v in (item.split('=') for item in parts[8].split(';') if '=' in item)}
                pos_key_st = (parts[0], parts[3], parts[6])  # chr, start, strand
                pos_key_end = (parts[0], parts[4], parts[6])  # chr, end, strand
                mrna_id = attrs['ID']
                
                # Get comparison attributes
                rank = float(attrs.get('Rank', float('inf')))
                identity = float(attrs.get('Identity', 0))
                
                # Determine best record by start position
                current_best_st = best_records_st.get(pos_key_st)
                if (current_best_st is None or 
                    rank < current_best_st['rank'] or 
                    (rank == current_best_st['rank'] and identity > current_best_st['identity'])):
                    
                    best_records_st[pos_key_st] = {
                        'rank': rank,
                        'identity': identity,
                        'id': mrna_id
                    }
                
                gene_blocks[mrna_id].append(line)
                all_entries.append(('GENE_BLOCK', mrna_id))
                
            # Handle other features (CDS/exon etc)
            else:
                attrs = {k:v for k,v in (item.split('=') for item in parts[8].split(';') if '=' in item)}
                if 'Parent' in attrs:
                    parent_id = attrs['Parent']
                    gene_blocks[parent_id].append(line)
                else:
                    all_entries.append(('OTHER', line))
    
    # Get mRNA IDs to keep (by start position)
    keep_ids_st = {rec['id'] for rec in best_records_st.values()}
    
    # Second phase: filter by end position
    best_records_end = {}
    filtered_entries = []
    temp_gene_blocks = defaultdict(list)
    
    for entry_type, content in all_entries:
        if entry_type in ('COMMENT', 'OTHER'):
            filtered_entries.append((entry_type, content))
        elif entry_type == 'GENE_BLOCK':
            mrna_id = content
            if mrna_id in keep_ids_st:
                for line in gene_blocks[mrna_id]:
                    if line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if parts[2] == 'mRNA':
                        attrs = {k:v for k,v in (item.split('=') for item in parts[8].split(';') if '=' in item)}
                        pos_key_end = (parts[0], parts[4], parts[6])  # chr, end, strand
                        mrna_id = attrs['ID']
                        
                        rank = float(attrs.get('Rank', float('inf')))
                        identity = float(attrs.get('Identity', 0))
                        
                        current_best_end = best_records_end.get(pos_key_end)
                        if (current_best_end is None or 
                            rank < current_best_end['rank'] or 
                            (rank == current_best_end['rank'] and identity > current_best_end['identity'])):
                            
                            if current_best_end:
                                temp_gene_blocks[current_best_end['id']].append('TO_DELETE')
                            
                            best_records_end[pos_key_end] = {
                                'rank': rank,
                                'identity': identity,
                                'id': mrna_id
                            }
                        
                        temp_gene_blocks[mrna_id].append(line)
                        filtered_entries.append(('GENE_BLOCK_END', mrna_id))
                    else:
                        temp_gene_blocks[mrna_id].append(line)
    
    # Get mRNA IDs to keep (by end position)
    keep_ids_end = {rec['id'] for rec in best_records_end.values()}
    
    # Write final output file (PAF lines removed)
    with open(output_path, 'w') as f:
        for entry_type, content in filtered_entries:
            if entry_type == 'COMMENT':
                f.write(content + '\n')
            elif entry_type == 'GENE_BLOCK_END':
                mrna_id = content
                if mrna_id in keep_ids_end:
                    for line in temp_gene_blocks[mrna_id]:
                        if line != 'TO_DELETE':
                            f.write(line + '\n')
            elif entry_type == 'OTHER':
                f.write(content + '\n')

if __name__ == '__main__':
    for filename in os.listdir(input_dir):
        if filename.endswith('.gff3'):
            input_file = os.path.join(input_dir, filename)
            output_file = os.path.join(output_dir, filename)
            print(f"Processing {filename}...")
            process_file(input_file, output_file)
    print("All files processed.")