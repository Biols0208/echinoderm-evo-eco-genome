#!/usr/bin/env bash
set -e

cwd=$PWD
indir=FastOMA2/proteome
outdir=$cwd/Orthomm_out
cpu=1
evalue=0.01
single_copy_threshold=0.5   ## taxon occupancy threshold for single-copy orthologs
inflation_value=1.5         ## Lower values are more permissive resulting in larger orthogroups. Higher values are stricter resulting in smaller orthogroups.
substitution_matrix=BLOSUM45    ## BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, PAM120, and PAM240. The default is BLOSUM62.
phmmer=phmmer
mcl=mcl

[ -d ${outdir} ] || mkdir ${outdir}

orthohmm ${indir} --phmmer ${phmmer} --evalue ${evalue} --substitution_matrix ${substitution_matrix} --cpu ${cpu} --single_copy_threshold ${single_copy_threshold} --mcl  ${mcl} --inflation_value ${inflation_value} --output_directory ${outdir}