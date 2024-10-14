#deeptools singularity container
deeptools_program="singularity run --app deeptools351 /share/singularity/images/ccs/conda/amd-conda8-rocky8.sinf"

bam_files=("unc_combined_PE_reps_sorted.bam" "bptf_combined_PE_reps_sorted.bam")

for i in ${bam_files[@]}
do
   echo "Processing $i" 
   name=$(basename "$i" _reps_sorted.bam)

   #Create bigwig file from bam files with bin size 5 and smooth length 60
   $deeptools_program bamCoverage -b $i \
   --binSize 5 --smoothLength 60 \
   --normalizeUsing RPKM \
   -o ${name}_bs5_SL60_RPKM.bw -p max
done
