#!/bin/bash

infile=$1
outdir=$2
batch=$3

if [ $# -eq 0 ]
  then
    echo "parameters to be provided: infile outdir batch"
fi

script_base="/home/icb/michaela.mueller/workspace/Benchmarking_data_integration/scripts/run_"
tools=(bbknn scanorama seurat)

for tool in ${tools[*]}
do
    outfile="$outdir$tool.h5ad"
    script="$script_base/$tool.py"
    echo "submitting to queue: $script -i $infile -o $outfile -b $batch"
    #qsub -q long_fed25 -hard -l job_mem=50G ./runqueue.sh python $script -i $infile -o $outfile -b $batch
done
