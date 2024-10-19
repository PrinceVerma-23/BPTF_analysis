#!/bin/bash

#singularity container with MACS2 program
macs2_program="singularity run --app macs22271 /share/singularity/images/ccs/conda/amd-conda8-rocky8.sinf"
macs2_out="/home/pve232/scratch/bptf_analyses/h3k4me3_bptf_analyses/macs2_out/"
genome_size=840196558

# Run MACS2 for each sample and replicate
samples=("bptf_X1" "unc_X1")
reps=("rep1" "rep2")

# Loop over each sample and replicate
for sample in ${samples[@]}; do
    for rep in ${reps[@]}; do
        # Define the treatment and control BAM file names
        treatment_bam=${sample}_k4me3_${rep}.sorted.nodups.bam
        control_bam=${sample}_input_${rep}.sorted.nodups.bam
        output_name=${sample}_${rep}
        
        # Run MACS2 peak calling
        $macs2_program macs2 callpeak -t $treatment_bam \
            -c $control_bam \
            -f BAM -g $genome_size \
            -n $output_name -q 0.01 --nomodel \
            --outdir $macs2_out
    done
done
