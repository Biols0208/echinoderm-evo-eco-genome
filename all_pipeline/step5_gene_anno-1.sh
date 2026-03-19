#!/usr/bin/env bash
set -e

export NXF_SINGULARITY_CACHEDIR=/egapx-0.3.1-alpha/NXF_SINGULARITY_CACHEDIR
export JAVA_HOME=/jdk-11.0.1
export TMPDIR=$PWD

file_path=local.yaml
outdir=Output
main_script=/egapx-0.3.1-alpha/ui/egapx.py
workdir=$PWD/workdir

[ -d egapx_config ] || mkdir -p egapx_config && cp /egapx-0.3.1-alpha/egapx_config/singularity.config egapx_config/singularity.config

python3 ${main_script} ${file_path} -e singularity -w ${workdir} -o ${outdir} -lc /egapx-0.3.0-alpha/support_data --func_name seastar

<<!!
local.yaml
genome: Anseropoda_placenta.fa
taxid: 290518
proteins: pep.fa
annotation_provider: Biols
annotation_name_prefix: A_placenta
locus_tag_prefix: A_placenta
!!