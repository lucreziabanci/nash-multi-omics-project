#!/bin/bash

# ATAC-seq peakcalling using MACS3
# Dataset:PRJNA725028

# 1. Peak calling
conda activate MACS3

for i in $(ls mapped | grep -E "\_sorted\.bam$" | perl -ne 's/_sorted\.bam//; print'); do
    outpath="try"
    echo -e "\nProcessing sample $i"
    macs3 callpeak -t mapped/$i"_sorted.bam" -f BAMPE -g hs --keep-dup all --outdir $outpath/$i"/" -n $i -q 0.01 --verbose 2 -B --broad --broad-cutoff 0.05
    # LEAVE BAMPE for paired-end data
done


