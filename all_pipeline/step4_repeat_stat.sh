#!/usr/bin/env bash
set -e

ln -s ../Asterina_phylactica.fa ./
faToTwoBit Asterina_phylactica.fa Asterina_phylactica.fa.2bit
/RepeatMasker/util/buildSummary.pl -species Asterina_phylactica -genome Asterina_phylactica.fa.2bit -useAbsoluteGenomeSize HiTE.update.out > HiTE.update.out.stat
/RepeatMasker/util/rmOutToGFF3.pl HiTE.update.out > HiTE.update.gff
repeat_to_gff.pl HiTE.update.out
grep -e Satellite -e Simple_repeat HiTE.update.out.gff | sed 's/Transposon/TandemRepeat/g' > HiTE.update.out.TR.gff
grep -v -e Satellite -e Simple_repeat HiTE.update.out.gff | sed 's/Transposon/TandemRepeat/g' > HiTE.update.out.noTR.gff
grep -v -e Satellite -e Simple_repeat HiTE.update.out > repeatmasker.out
cat HiTE.update.out.TR.gff HiTE.update.out.noTR.gff > all.gff
ln -s HiTE.update.out.noTR.gff all_without_trf.gff
ln -s HiTE.update.out.noTR.gff repeatmasker.gff
ln -s HiTE.update.out.TR.gff trf.gff

perl stat.pl -repeatmasker -trf Asterina_phylactica.fa
awk '{(tot+=$2)};END{print tot}' Asterina_phylactica.fa.len

genome_len=$(awk '{(tot+=$2)};END{print tot}' Asterina_phylactica.fa.len)
distribution_TE_divergence.pl repeatmasker.out $genome_len -X_step 10 -X_end 40 -Y_step 0.5 -Y_end 3

perl /Parsing-RepeatMasker-Outputs/parseRM_simple.pl -genfile Asterina_phylactica.fa -RMout repeatmasker.out

awk -F "[;\t]" '{print $6"\t"$10"\t"$11"\t"$12}' all_without_trf.gff | sed 's/Class=//g' | sed 's/PercDiv=//g;s/Target=//g' | awk '{print $5"\t"$2"\t"$1"\t"$NF}' | awk '{if ($1~/DNA/)print "DNA\t"$0;else if ($1~/LTR/)print "LTR\t"$0;else if ($1~/LINE/)print "LINE\t"$0;else if ($1~/SINE/)print "SINE\t"$0 }' | grep -v pValue | csvtk -t add-header -n class,type,fam,len,div > Sum.txt

awk '{print $5"\t"$4}' Sum.txt > Sum.huoli.txt
awk '{print $1"\t"$5}' Sum.txt | sed '1d' > Sum.div.density.txt

huoli.r Sum.huoli.txt
density.r Sum.div.density.txt class div

sed '1d' Sum.huoli.txt > Sum.huoli2.txt
huoli2.r Sum.huoli2.txt

rm tmp1 tmp2 repeatmasker.gff.type.tmp denovo_repeat.statistics.xls Asterina_phylactica.fa.Nrem Asterina_phylactica.fa.Nrem.fa Asterina_phylactica.fa.Nrem.fa.index Asterina_phylactica.fa.Nrem.fa.length repeatmasker.out.parseRM.1.log Sum.huoli2.txt Asterina_phylactica.fa.2bit Asterina_phylactica.fa.len
repeatmasker.out.parsed1/repeatmasker.out.parseRM.Summary.tab
pigz --best -p 1 HiTE.update.gff HiTE.update.out HiTE.update.out.gff HiTE.update.out.TR.gff HiTE.update.out.noTR.gff repeatmasker.out Sum.huoli.txt Sum.txt Sum.div.density.txt repeatmasker.out.parsed1/*tab
sed '1d' TE_repeat.statistics.xls | awk '{print $1"\t"$NF}' | csvtk transpose -t -H | cut -f 2-7 | paste - <(zcat repeatmasker.out.parsed1/repeatmasker.out.parseRM.Summary.tab.gz | grep -A 1 "%masked_minus_%in-double" | awk '{print $NF}') > tongji.sum