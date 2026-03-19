#!/usr/bin/env bash
set -e

## for genome assembly, by hifiasm-0.19.9
# *.clean.fq.gz are Hi-c data
# HCS.fastq.gz is Hi-Fi data
hifiasm -o HSC.asm -t48 -z20 --h1 HS_R1.clean.fq.gz --h2 HS_R2.clean.fq.gz HCS.fastq.gz
