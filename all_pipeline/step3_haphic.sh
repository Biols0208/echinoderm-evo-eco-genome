#!/usr/bin/env bash
set -e

## for genome scaffolding

genome=$PWD/nextpolish.fa
hic1=$PWD/HS_R1.clean.fq.gz
hic2=$PWD/HS_R2.clean.fq.gz
cpu=48

### index
bwa index ${genome}

### map
bwa mem -5SP -t ${cpu} ${genome} ${hic1} ${hic2} | samblaster | samtools view - -@ ${cpu} -S -h -b -F 3340 -o HiC.bam

### filter
filter_bam HiC.bam 1 --nm 3 --threads ${cpu} | samtools view - -b -@ ${cpu} -o HiC.filtered.bam

### scaffolding
haphic pipeline ${genome} HiC.filtered.bam ${nchrs} --threads ${cpu} --processes ${cpu} --max_inflation 50 --min_inflation 2 --inflation_step 0.5