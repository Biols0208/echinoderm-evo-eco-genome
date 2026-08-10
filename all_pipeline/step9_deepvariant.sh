#!/bin/bash
TF_ENABLE_ONEDNN_OPTS=0
set -euo pipefail

sampleID=SRR28413286
OUTPUT_DIR=$PWD
ref=$PWD/Acanthaster_planci.fa
bam=${OUTPUT_DIR}/SRR28413286.bam
out_dir=${OUTPUT_DIR}/deepvariant_vcf
out_gvcf=${out_dir}/${sampleID}.g.vcf.gz
out_vcf=${out_dir}/${sampleID}.vcf.gz
cpu=20
sif=deepvariant.1.8.0_saved_model.sif
TMPdir=${OUTPUT_DIR}/tmp_${sampleID}

[ -d ${TMPdir} ] || mkdir -p ${TMPdir}
[ -d ${out_dir} ] || mkdir -p ${out_dir}

singularity exec -B /ldfsqd1/:/ldfsqd1/ -B /usr/lib/locale/:/usr/lib/locale/ ${sif} \
     /opt/deepvariant/bin/run_deepvariant \
     --model_type=WGS \
    --vcf_stats_report=true \
    --ref=${ref} \
    --reads=${bam} \
    --output_vcf=${out_vcf} \
    --output_gvcf=${out_gvcf} \
    --intermediate_results_dir=${TMPdir} \
    --num_shards=${cpu}

rm -rf ${TMPdir}
