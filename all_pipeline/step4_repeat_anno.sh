#!/usr/bin/env bash
set -e

## genome annotation: repeat

genome=$PWD/Anseropoda_placenta.fa	### absolute paths
out_dir=$PWD/Hite_out		### absolute paths
cpu=48
curated_lib=/Full_path/animal.rmdup.lib		## Provide a fully trusted curated library, which will be used to pre-mask highly homologous sequences in the genome.
isplant=0		## Is it a plant genome, 1: true, 0: false.
isremove_nested=1	## Whether to remove nested TE, 1: true, 0: false.
isrecover=1		## whether to enable recovery mode to avoid starting from the beginning, 1: true, 0: false.
isdomain=1		## Whether to obtain TE domains, HiTE uses RepeatPeps.lib from RepeatMasker to obtain TE domains, 1: true, 0: false.
isannotate=1		## Whether to annotate the genome using the TE library generated, 1: true, 0: false.
isintact_anno=1		## Whether to generate annotation of full-length TEs, 1: true, 0: false.
isBM_RM2=1		## Whether to conduct benchmarking of RepeatModeler2, 1: true, 0: false.
isBM_HiTE=1		## Whether to conduct benchmarking of HiTE, 1: true, 0: false.


/usr/bin/singularity exec --bind HiTE.sif \
	python /HiTEv3.2/main.py \
	--genome ${genome} \
	--thread ${cpu} \
	--outdir ${out_dir} \
	--chunk_size 200 \
	--plant ${isplant} \
 	--remove_nested ${isremove_nested} \
	--domain ${isdomain} \
	--recover ${isrecover} \
	--annotate ${isannotate} \
	--intact_anno ${isintact_anno} \
	--BM_RM2 ${isBM_RM2} \
	--BM_HiTE ${isBM_HiTE} \
	--curated_lib ${curated_lib}