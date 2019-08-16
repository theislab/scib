#!/bin/bash

indir=$1
outdir=$2
ls $indir
for fname in $indir/*
do
    base="$(basename $fname)"
    outfile="$outdir/$base"
    echo $outfile
    qsub -q long_fed25 -hard -l job_mem=50G ./runqueue.sh python /home/icb/daniel.strobl/Benchmarking_data_integration/scripts/preprocessing.py -i $fname -o $outfile 
done
