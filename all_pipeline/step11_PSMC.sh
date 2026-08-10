#!/bin/bash
set +o posix

REF=Ophiura_sarsii.fa
VCF=SRR27344560.vcf.gz

[ -e ${REF}.fai ] || ln -s ../${REF}.fai ./
[ -e ${REF} ] || ln -s ../${REF} ./
[ -e merge.keep.bed ] || ln -s ../merge.keep.bed .

bcftools view -m2 -M2 -v snps $VCF -Oz -o SNP.vcf.gz
bcftools index SNP.vcf.gz

bcftools view -i 'GQ>=20 && DP>30 && DP<110' SNP.vcf.gz -O z -o SNP.filted.vcf.gz

bcftools index SNP.filted.vcf.gz

bcftools view -T merge.keep.bed SNP.filted.vcf.gz | bgzip > keep.vcf.gz

bcftools index keep.vcf.gz

cat $REF | bcftools consensus --iupac-codes keep.vcf.gz > consensus.fa

/01_software/psmc/utils/fq2psmcfa -q20 consensus.fa > 1.psmcfa

/01_software/psmc/psmc -N25 -t15 -r5 -p "6+20*2+4+6" -o 1.psmc 1.psmcfa
