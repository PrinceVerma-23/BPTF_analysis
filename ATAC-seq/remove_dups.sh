samtools_program="singularity run --app samtools113 /share/singularity/images/ccs/ngstools/samtools-1.13+matplotlib-bcftoools-1.13.sinf"
picard="singularity run --app picard2264 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf java -jar /usr/local/picard-2.26.4/picard.jar"

## remove duplicates from BAM files ##

for i in *filt.sub.sorted.bam
do 
   echo "Processing sample $i"
   name=$(basename "$i" .bam)
   
   # Remove duplicates
   $picard MarkDuplicates -I $i -O $name.NODUPS.bam -M $name.NODUPS.txt -REMOVE_DUPLICATES true

   #Create a coordinate sorted bam file
   $samtools_program samtools sort -@ 64 $name.NODUPS.bam -o $name.NODUPS.sorted.bam

   #Create index file from bam
   $samtools_program samtools index $name.NODUPS.sorted.bam
done
