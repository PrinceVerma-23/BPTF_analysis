#!/bin/bash
# Script was ran using SLURM scheduler.

#file with length of each chromosome in genome
chr_length="/home/pve232/scratch/bptf_analyses/chr_length.txt"

# diffreps program loc
diffrep_path="/home/pve232/scratch/programs/diffreps"

# output directory for diffreps results
diffrep_out="/home/pve232/scratch/bptf_analyses/diffrep_output"

# perl module is required to load to run diffrep script
module load perl-5.34.0-gcc-8.4.1-4mv6qmh

# window sizes
window_sizes=("1000" "500" "1500")

for window in "${window_sizes[@]}"; do
  perl $diffrep_path/diffrep.pl --treatment bptf_X1_k4me3_rep1.bed bptf_X1_k4me3_rep2.bed \
  --control unc_X1_k4me3_rep1.bed unc_X1_k4me3_rep2.bed \
  --btr bptf_X1_input_rep1.bed bptf_X1_input_rep2.bed \
  --bco unc_X1_input_rep1.bed unc_X1_input_rep2.bed \
  --nsd sharp --noanno --window $window --nohs --frag 250 \
  --report $diffrep_out/out_bptf_unc_diffrep_${window}.tsv --nproc 12 --chrlen $chr_length
done
