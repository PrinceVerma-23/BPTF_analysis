#!/bin/bash

# Script was ran using SLURM scheduler.
# singularity container for bowtie2 and samtools
bowtie2_program="singularity run --app bowtie2244 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf"
samtools_program="singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"

# Bowtie2 index files
smed_h1_idx="/home/pve232/scratch/new_chr_genome/schmed_h1_bowtie2_idx/schmed_h1"


##bptf chip align

for i in *.fastq.gz;
do
    name=$(basename $i .fastq.gz)
    echo $i

    echo "processing $name"

    #Run bowtie2 and capture the standard output
    align_out=$($bowtie2_program bowtie2 -p 4 -x $smed_h1_idx -q $i -S $name.sam 2>&1)

    #Print sample name
    echo "Sample: $name"

    #Print alignment output
    echo "$align_out"

done
