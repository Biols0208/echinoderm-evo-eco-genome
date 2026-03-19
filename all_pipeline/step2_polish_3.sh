#!/usr/bin/env bash
set -e

## polish by nextpolish

nextPolish nextpolish.cfg

<<!!
nextpolish.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 24
multithread_jobs = 2
genome = ./HSC.asm.hic.p_ctg.fasta #genome file
genome_size = auto
workdir = ./rundir
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa

[hifi_option]
hifi_fofn = ./hifi.fofn
hifi_options = -min_read_len 1k -max_depth 100
hifi_minimap2_options = -x map-pb

sgs.fofn
/Full_PATH/HCS-ch-TYD20223773_R1.fq.gz
/Full_PATH/HCS-ch-TYD20223773_R2.fq.gz

hifi.fofn
/Full_PATH/HCS.fastq.gz