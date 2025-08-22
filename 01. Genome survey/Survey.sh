# Count k-mers
jellyfish count -g generators -G 2 -m 21 -s 5G -t 128 -C -o reads.jf <(zcat "$fq1") <(zcat "$fq2")
# Generate histogram
jellyfish histo -h 65535 -t 128 reads.jf > reads.histo
# Genome stats
genomescope reads.histo 21 150 genomescope_output
