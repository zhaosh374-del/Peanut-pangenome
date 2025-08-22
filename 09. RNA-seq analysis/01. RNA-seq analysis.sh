# Data filtering
fastp -i ${raw_1} -I ${raw_2} -o ${clean_1} -O ${clean_2}

# Build index
extract_splice_sites.py ref.gtf > ref.ss
extract_exons.py ref.gtf > ref.exon
hisat2-build --ss ref.ss \
               --exon ref.exon \
               ref.fa \
               ref

# Alignment
hisat2 -p 128 --dta \
    -x ref \
    --summary-file summary.txt \
    -1 ${clean_1} \
    -2 ${clean_2} \
    -S ${sam}

# Convert SAM to BAM
samtools view -S ${sam} -b > ${bam}
samtools sort ${bam} -o ${sorted_bam}
samtools index ${sorted_bam}

# Transcript assembly
stringtie -p 128 -G ref.gtf \
            -l sample \
            -o ${assembl_data}/sample.gtf \
            ${sorted_bam}

# Create mergelist file
ls -l ${assembl_data}/*.gtf | awk '{print $9}' > ${assembl_data}/mergelist.txt

# Transcript merging
stringtie --merge -p 10 -G ref.gtf \
                    -o ${assembl_data}/merged.gtf \
                    ${assembl_data}/mergelist.txt

# Quantification
stringtie -p 128 -e -G ${assembl_data}/merged.gtf \
             -o ${quantity_data}/sample.gtf \
             -A ${quantity_data}/gene_abundances.tsv \
             ${sorted_bam}

# Extract merged count results (The sample_list.txt file specifies the paths to the GTF files for each sample)
python prepDE.py -i sample_list.txt

# Extract TPM values
perl stringtie_expression_matrix.pl --expression_metric=TPM \
                                    --result_dirs='results' \
                                    --transcript_matrix_file=transcript_tpms_all_samples.tsv \
                                    --gene_matrix_file=gene_tpms_all_samples.tsv

# prepDE.py is available from http://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
# stringtie_expression_matrix.pl is available from https://github.com/griffithlab/rnaseq_tutorial/blob/master/scripts/stringtie_expression_matrix.pl