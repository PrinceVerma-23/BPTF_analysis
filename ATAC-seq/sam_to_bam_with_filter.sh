#!/bin/bash

samtools_program="singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
picard="singularity run --app picard2264 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf java -jar /usr/local/picard-2.26.4/picard.jar"

for i in *.sam
do
   echo $i
   name=$(basename "$i" .sam)

   #Make bam files with header included
   $samtools_program samtools view -b -h $i > $name.unsorted.bam

   #Create a sorted bam file
   $samtools_program samtools sort -@ 64 $name.unsorted.bam -o $name.sorted.bam

   #Run flagstat before doing filter step
   $samtools_program samtools flagstat $name.sorted.bam > $name.before.filter.flagstat

   #restrict to properly pair read
   $samtools_program samtools view -h -b -f 3 $name.sorted.bam > $name.filt.bam

   #Create a sorted bam file
   $samtools_program samtools sort -@ 64 $name.filt.bam -o $name.sorted.filt.bam

   #Create index file from bam
   $samtools_program samtools index $name.sorted.filt.bam

   #Run flagstat after doing filter step
   $samtools_program samtools flagstat $name.sorted.filt.bam > $name.after.filter.flagstat

   rm $name.unsorted.bam $name.sorted.bam $name.filt.bam 

done
