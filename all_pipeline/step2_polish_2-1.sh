#!/usr/bin/env bash
set -e

## polish by kmerdup

prefix=hcs
genome=genome.fa
mer_len=19
cpu=60
counter_len=8
NGS_fq1=HCS-ch-TYD20223773_R1.fq.gz
NGS_fq2=HCS-ch-TYD20223773_R2.fq.gz
ccs_fa=HCS.fasta.gz

## step1 count k-mers
jellyfish count -m ${mer_len} -s 100M -t ${cpu} -c ${counter_len} -C -o ${prefix}.count <(pigz -p 5 -d -c ${NGS_fq1}) <(pigz -p 5 -d -c ${NGS_fq2}) <(zcat ${ccs_fa})

## if jellyfish gives you more than one count file, you need merge them firs
# jellyfish merge -v -o ${prefix}.count.jf ${prefix}.count_*

mv ${prefix}.count ${prefix}.count.jf

## step2 stat and histo (you can skip)
#jellyfish stats -o ${prefix}.stats ${prefix}.count.jf
#jellyfish histo -t ${cpu} ${prefix}.count.jf | perl -lane 'my ($dpt, $cnt) = split(/\s+/, $_); my $nn = $dpt * $cnt;print "$dpt\t$cnt\t$nn"' > ${prefix}.histo

## step3 dump k-mers
jellyfish dump -c -t -o ${prefix}.dump.all ${prefix}.count.jf
perl /kmerDedup/kmerFilter.pl -d ${prefix}.dump.all -o ${prefix}.filt.fa -l 3 -u 100000
perl /kmerDedup/splitFasta.pl -f ${prefix}.filt.fa -o split -k ${prefix}.kmer

## step4 mapping k-mers
/dellfsqd2/ST_OCEAN/USER/lishuo11/01_soft/mambaforge/bin/perl /kmerDedup/kmerDedup/fa2fa.pl -f ${genome} -o ${prefix}.format.fa -c F -n F -l 1000

bowtie2-build --threads ${cpu} ${prefix}.format.fa ${prefix}.format
ls split/ > split.id
[ -d shell ] || mkdir shell

## vim bowtie2_demo.sh
<<!!
prefix=HCS
cpu=5
work_dir=$PWD/../
input_dir=${work_dir}/split
output_dir=${work_dir}/mapping

[ -d ${output_dir} ] || mkdir ${output_dir}
bowtie2 --very-sensitive -k 1000 --score-min L,-0.6,-0.2 --end-to-end --reorder -L 21 --rg-id ${prefix} --rg SM:${prefix} -p ${cpu} -f ${input_dir}/ABCD -x ${work_dir}/${prefix}.format | samtools view -@ ${cpu} -F 4 -bS - > ${output_dir}/ABCD.bam
!!

## make shells
for i in $(cat split.id); do sed "s/ABCD/${i}/g" bowtie2_demo.sh > shell/bowtie2_${i}.sh; done