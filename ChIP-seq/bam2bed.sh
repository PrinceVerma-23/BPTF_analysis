#!/bin/bash

# Script was ran using SLURM scheduler.
# Convert all bam files to bed files.
# singularity container for bedtools
bedtools_program="singularity run --app bedtools2300 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf"

for bams in *.bam
do

  echo "Processing $bams"
  name=$(basename "$bams" .sorted.nodups.bam)
  $bedtools_program bedtools bamtobed -i $bams > $name.bed

done
