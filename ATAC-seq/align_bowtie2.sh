#bowtie2 singularity container loc
bowtie2_program="singularity run --app bowtie2244 /share/singularity/images/ccs/conda/amd-conda2-centos8.sinf"

#Bowtie2 index
smed_h1_idx=/home/pve232/scratch/new_chr_genome/schmed_h1_bowtie2_idx/schmed_h1


for i in {unc_rep1,unc_rep2,bptf_rep1,bptf_rep2}
do
    # input file names
    forward_read="${i}_R1.fastq.gz"
    reverse_read="${i}_R2.fastq.gz"
    
    echo "Processing sample: $i"

    # output filename
    sam_files="${i}.sam"

    #Run bowtie2 alignment
    align_out=$($bowtie2_program bowtie2 -p 64 --very-sensitive -X 1000 \
        -x $smed_h1_idx \
        -1 $forward_read -2 $reverse_read \
        -S $sam_files 2>&1)

    #Print alignment output
    echo "ATAC-seq sample: $i alignment output"
    echo "$align_out"
done
