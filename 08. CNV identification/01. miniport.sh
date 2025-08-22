# Create index
miniprot -t 128 -d ${sample_name}.mpi ${sample_name}.fasta

# Alignment
miniprot -t 128 --gff ${sample_name}.mpi ${protein_dir} > ${sample_name}.gff3

# gffread
gffread ${sample_name}.gff3 -g ${sample_name}.fasta -y ${sample_name}.protein.fasta

# Create BLAST database
makeblastdb -in ${Tifrunner_protein} -dbtype prot -out ${Tifrunner_protein}

# Alignment
blastp -query ${sample_name}.protein.fasta -db ${Tifrunner_protein} -out ${sample_name}.txt -outfmt 6 -evalue 1e-5 -num_threads 128

