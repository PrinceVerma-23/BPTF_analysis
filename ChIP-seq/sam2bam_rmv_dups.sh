#!/bin/bash

# Script was ran using SLURM scheduler.
# singularity container for bowtie2 and samtools
samtools_program="singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
picard="singularity run --app picard2264 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf java -jar /usr/local/picard-2.26.4/picard.jar"

# Convert sam to bam, sort bam,run flagstat before and after duplicate removal, and remove duplicates

for sam in *.sam; 
do
	echo $sam
	name=$(basename $sam .sam)
	
	echo " Processing: $name"

	#convert sam to bam
	$samtools_program samtools view -@ 64 -bh -o $name.unsorted.bam $sam

  #Sort bam
  $samtools_program samtools sort -@ 64 $name.unsorted.bam -o $name.sorted.bam

  #Run flagstat before duplicates removal
	$samtools_program samtools flagstat $name.sorted.bam > $name.flagstat

  #Run Picard to mark duplicates
	$picard MarkDuplicates -I $name.sorted.bam -O $name.sorted.nodups.bam -M $name.sorted.nodups.txt -REMOVE_DUPLICATES true

	#Create index for bam files
	$samtools_program samtools index $name.sorted.nodups.bam

	#Run flagstat
  $samtools_program samtools flagstat $name.sorted.nodups.bam > $name.sorted.nodups.bam.flagstat

  rm $name.unsorted.bam

  echo "Processing for $name completed."

done
