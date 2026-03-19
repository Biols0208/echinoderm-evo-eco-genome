#!/usr/bin/env bash
set -e

## polish by craq
genome=YZ.asm.hic.p_ctg.fasta               # genome fasta file from step1
sms_fq=yz.fastq.gz                          # Hi-Fi data
ngs_fq=YZ-C-3_R1.fq.gz,YZ-C-3_R2.fq.gz      # WGS data
cpu=20
outdir=Result_CRAQ

craq -g ${genome} -sms ${sms_fq} -ngs ${ngs_fq} -x map-hifi --plot T --break T --thread ${cpu} --output_dir ${outdir}

### Get user-specified regional(i.e. window=50000) AQI score.
cat ${outdir}/runAQI_out/strER_out/out_final.CSE.bed ${outdir}/runAQI_out/locER_out/out_final.CRE.bed > ${outdir}/runAQI_out/CRE_CSE.bed
window=50000
perl /CRAQ/src/regional_AQI.pl ${outdir}/seq.size ${window} ${window} ${outdir}/runAQI_out/CRE_CSE.bed > ${outdir}/runAQI_out/plot_AQI.out
### plot.  the scaffolds ids you want to present is ok. see CRAQcircos.py --help
python /CRAQ/src/CRAQcircos.py --genome_size ${outdir}/seq.size --genome_error_loc ${outdir}/runAQI_out/CRE_CSE.bed --genome_score ${outdir}/runAQI_out/plot_AQI.out --output ${outdir}/runAQI_out/plot_AQI.out.pdf