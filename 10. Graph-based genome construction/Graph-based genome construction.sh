jobstore=../1.cactus_jobstore
outdir=../1.cactus_outdir
workdir=../1.cactus_workdir
cactus-pangenome ${jobstore} genomes.list --binariesMode local --workDir $workdir \
--outDir $outdir --outName graph_pan --reference ref_name --mgCores 128 --logFile cactus.log --vcf --giraffe --gfa --gbz --viz

vg minimizer -t 52 -d graph_pan.d2.dist -o graph_pan.d2.shortread.withzip.min -z graph_pan.d2.shortread.zipcodes graph_pan.d2.gbz
vg giraffe -m graph_pan.d2.shortread.withzip.min -d graph_pan.d2.dist -Z graph_pan.d2.gbz \
-f ${sample}_1_clean.fq.gz -f ${sample}_2_clean.fq.gz -N ${sample} -t 128 > ${sample}.d2.gam
vg pack -x graph_pan.d2.gbz -g ${sample}.d2.gam -Q 5 -o ${sample}.d2.pack -t 128
vg snarls -t 128 graph_pan.d2.gbz > graph_pan.d2.snarls
vg call graph_pan.d2.gbz -r graph_pan.d2.snarls -k  ${sample}.d2.pack -a -A -t 128 -s ${sample} > ${sample}.d2.vcf

vcf=`ls *vcf.gz | tr "\n" " "`
bcftools merge -o graph_pan.d2.vcf.gz $vcf

vcftools --gzvcf graph_pan.d2.vcf.gz --mac 5 --minDP 4 --minQ 100 --max-missing 0.8 --recode --stdout | bgzip -c > graph_pan_mr0.8_vcf.gz
tabix graph_pan_mr0.8_vcf.gz
