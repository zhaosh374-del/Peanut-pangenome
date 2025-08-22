python quartet_centrominer.py -i genome.fasta --TE genome.TEanno.gff3 -n 200 -m 400 -t 20

TRASH_run.sh genome.fasta --par=52 --o /TRASH/

mafft --auto Monomer sequences.fasta > mafft-Monomer sequences.fasta
iqtree3  --threads-max 52 -s mafft-Monomer sequences.fasta -bb 1000

makeblastdb -in Monomer sequences.fasta -dbtype nucl -out Monomer sequences -parse_seqids

blastn -query Monomer sequences.fasta -out Monomer sequences.blast -db Monomer sequences -outfmt 6 -evalue 1e-5 -num_threads 10

cd-hit-est -i Monomer sequences.fasta -o Monomer sequences.fasta.cdhit.fa -c 0.8 -n 10 -d 0 -M 16000 -T 10