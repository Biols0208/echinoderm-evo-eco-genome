#!/usr/bin/env bash
set -e

[ -d TMPdir ] || mkdir TMPdir && chmod -R 777 TMPdir
export JAVA_HOME=/jdk-11.0.1
export TMPDIR=$PWD/TMPdir

inputdir=FastOMA2
outdir=FastOMA2/result
databse=/OMArk_database/LUCA-v2.0.0.h5
config=/FastOMA/nextflow.config

nextflow -C ${config} run FastOMA.nf -profile standard --input_folder ${inputdir} --output_folder ${outdir} --omamer_db ${databse}
rm -rf TMPdir .nextflow* work