#!/usr/bin/env bash
set +o posix
set -eo pipefail

genome=Luidia_sarsii.fa
bam=ERR10934081.sort.bam
prefix=ERR10934081

angsd \
  -i ${bam} -ref ${genome} \
  -GL 1 -doSaf 1 -anc ${genome} \
  -out ${prefix}

realSFS ${prefix}.saf.idx > ${prefix}.sfs

realSFS saf2theta ${prefix}.saf.idx -sfs ${prefix}.sfs -outname ${prefix}

thetaStat do_stat ${prefix}.thetas.idx
sed '1d' ${prefix}.thetas.idx.pestPG | awk '$NF!=0 {print $0"\t"$5/$NF}' | cat <(head -n 1 ${prefix}.thetas.idx.pestPG) - > pi.deal.txt

thetaStat do_stat ${prefix}.thetas.idx -win 10000 -step 10000 -outnames ${prefix}.pi.window
sed '1d' ${prefix}.pi.window.pestPG | awk '$NF!=0 {print $0"\t"$5/$NF}' | cat <(head -n 1 ${prefix}.pi.window.pestPG) - > pi.window.deal.txt

angsd \
  -i ${bam} -ref ${genome} \
  -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 \
  -out ${prefix}

zcat ${prefix}.geno.gz | awk 'BEGIN{het=0;tot=0} {tot++; if($3==1) het++} END{print "Heterozygosity =", het/tot, " (het=",het,", total=",tot,")"}' > Heterozygosity.txt
sed '1d' pi.window.deal.txt | awk '{print $NF}' | datamash mean 1 > pi.window.deal.mean.txt
sed '1d' pi.deal.txt | awk '{print $NF}' | datamash mean 1 > pi.deal.mean.txt

# rm *geno.mafs.gz *geno.geno.gz *.thetas.gz *.saf.pos.gz *.saf.gz
