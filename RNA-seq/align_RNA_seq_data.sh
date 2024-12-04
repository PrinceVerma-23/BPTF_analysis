#!/bin/bash

# Script was ran using SLURM scheduler.
# singularity container for hisat2, samtools, and stringtie
hisat2_program="singularity run --app hisat2221 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf"
samtools_program="singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
stringtie_program="singularity run --app stringtie217  /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf"

#Path to files
chr_genome_fasta="/home/pve232/scratch/new_chr_genome/schMedS3_h1.fa"
smed_chr_hisat2_idx="/home/pve232/scratch/new_chr_genome/schmed_h1_hisat2_idx/schmed_h1"
gtf_chr_file="/home/pve232/scratch/new_chr_genome/schMedS3_h1.gtf"

#build Hisat2 index for schMedS3_h1 chr level genome
$hisat2_program hisat2-build -f $chr_genome_fasta $smed_chr_hisat2_idx

#Day9 files
day9_fastqs="/home/pve232/scratch/bptf_analyses/bulk_RNA_X1_unc_bptf/day9"
day9_stringtie_out="/home/pve232/scratch/bptf_analyses/bulk_RNA_X1_unc_bptf/day9/"

#day12 files
day12_fastqs="/home/pve232/scratch/bptf_analyses/bulk_RNA_X1_unc_bptf/day12"
day12_stringtie_out="/home/pve232/scratch/bptf_analyses/bulk_RNA_X1_unc_bptf/day12/"


# function to process RNA-seq data
process_rnaseq() {
    fastqs_dir=$1
    output_dir=$2

    for i in $fastqs_dir/*.fastq.gz; 
    do
        echo "Processing $i"
        name=$(basename "$i" .fastq.gz)

        # Align
        $hisat2_program hisat2 -q -x $smed_chr_hisat2_idx -U $i -S $name.sam

        # Convert SAM to BAM
        $samtools_program samtools sort -@ 64 -o $name.bam $name.sam
        rm $name.sam

        # Quantification with StringTie
        $stringtie_program stringtie $name.bam -p 64 -e -G $gtf_chr_file \
              -A $output_dir/$name.gene_abundances.tsv \
              -o $output_dir/$name.gtf

    done
}

# Process Day 9 and Day 12
process_rnaseq $day9_fastqs $day9_stringtie_out
process_rnaseq $day12_fastqs $day12_stringtie_out
