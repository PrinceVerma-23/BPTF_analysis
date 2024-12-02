#!/bin/bash

macs2_program="singularity run --app macs22271 /share/singularity/images/ccs/conda/amd-conda8-rocky8.sinf"

# Effective genome_size=840196558 total non-N bases

# Output directories for MACS2
macs_out="/home/pve232/scratch/bptf_analysis/atac_seq/out_macs2"
bedpe_dir="/home/pve232/scratch/bptf_analysis/atac_seq/out_BEDPE"

# Run Macs2 to call peak for ATAC samples
for i in $bedpe_dir/*.minimal.bedpe
do
  echo $i
  name=$(basename "$i" | cut -d '.' -f 1)
  echo $name

  # Run macs2 callpeak for all samples
  $macs2_program macs2 callpeak -t $i \
    -f BEDPE \
    -g 840196558 \
    -n $name 
    -q 0.05 \
    --keep-dup all \
    --outdir $macs_out \
done
