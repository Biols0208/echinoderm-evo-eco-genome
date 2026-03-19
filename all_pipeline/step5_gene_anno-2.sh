#!/usr/bin/env bash
set -e

genome=HCS_chr.fa
pep=pep.for_anno.fa
prefix=EFGHABCD_                  ### prefix for IDs in GFF3
cpu=15
max_intron_size=20k            ### max intron size [200k]
splice_model=1                  ### splice model: 2=mammal, 1=general, 0=none (see Detail) [1]
weight_of_splice_penalty=1      ### weight of splice penalty; 0 to ignore splice signals [1]

miniprot -G ${max_intron_size} -j ${splice_model} -t ${cpu} --gff -P ${prefix} -C ${weight_of_splice_penalty} ${genome} ${pep} --outs=0.99 > EFGHABCD.gff

grep -A 1 "##PAF" EFGHABCD.gff | awk '$1!~/--/' | paste - - | awk '($5-$4)/$3>0.9' | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$(NF-10)"\t"$(NF-6)"\t"$(NF-5)"\t"$(NF-4)"\t"$(NF-2)}' | sed 's/;/\t/g;s/ID=//g' > filter.ABCD.info
cat EFGHABCD.gff | grep -v -e "^#" -e stop_codon | gffread -C -G -K -Q -Y -M --cset -d dup -H -V -P -N -Z - -g ${genome} > EFGHABCD.gff.gffread
cut -f 9 filter.ABCD.info | fishInWinter.pl -bf table -ff gff - EFGHABCD.gff.gffread | awk '{print $0";"}' | sed "s/miniprot/miniprot_EFGHABCD/1" > EFGHABCD.gff.gffread.gff
Covert_for_evm.pl EFGHABCD.gff.gffread.gff miniprot | awk '!a[$1"\t"$4"\t"$5]++{print $0}' > EFGHABCD.gff.gffread.gff.forevm.gff3