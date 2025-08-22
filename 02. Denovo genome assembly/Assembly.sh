1. hifiasm -o ${Sample}.asm -t 40 ${Sample}.ccs.fasta

2. hifiasm -o ${Sample}.asm -t 128 --ul UL.fq.gz ${Sample}.ccs.fasta

3. HiC-Pro -c config-hicpro.txt -o Analysis -i HIC_read

4. python2 HiCPlotter/HiCPlotter.py -f iced/1000000/sample1_1000000_iced.matrix -o ${Sample} -r 1000000 -tri 1 -bed sample1/raw/1000000/sample1_1000000_abs.bed -n ${Sample} -wg 1 -chr Chr11

5. ##get contig length
   perl fastaDeal.pl -attr id:len ${Sample}.contigs.fa > ${Sample} .contigs.fa.len 
   ##draw contig Hi-C heatmaps with 10*100000 (1-Mb) resolution
   matrix2heatmap.py ${Sample}_HiC_100000_abs.bed ${Sample}_HiC_100000.matrix 10
   ##Run one round, when the contig assembly is quite good
   perl endhic.pl ${Sample}.contigs.fa.len ${Sample}_HiC_100000_abs.bed ${Sample}_HiC_100000.matrix ${Sample}_HiC_100000_iced.matrix
   ln Round_A.04.summary_and_merging_results/z.EndHiC.A.results.summary.cluster* ./
   ##convert cluster file to agp file
   perl cluster2agp.pl Round_A.04.summary_and_merging_results/z.EndHiC.A.results.summary.cluster ${Sample}.contigs.fa.len > ${Sample}.scaffolds.agp
   ##get final scaffold sequence file
   perl agp2fasta.pl ${Sample}.scaffolds.agp ${Sample}.contigs.fa > ${Sample}.scaffolds.fa

6. ragtag.py scaffold $REF $SCF -o ./ragtag_default -t 26
   ragtag.py scaffold $REF $SCF -o ./ragtag_nucmer -t 26 --aligner 'nucmer'


