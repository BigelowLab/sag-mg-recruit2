

#decimal_numbers=(100 97.5 95 92.5 90 85 80 70)
decimal_numbers=(80)

#targets=(sag mag)
targets=(sag mag)



for target in "${targets[@]}"; do
    for mg1 in /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/metag/1M_subsample/*1_1Msub.fasta.gz; do

    mg2="${mg1/1_paired.fasta.gz/2_1Msub.fasta.gz}"
    file_name=$(basename "$mg2")

    # Extract the part before the first underscore
    first_part="${file_name%%_*}"

    
    for decimal in "${decimal_numbers[@]}"; do

    outlog=/mnt/s1/labs/brown/sag-mag/out_logs/${target}${decimal}${first_part}.out
    qsub_cmd="qsub -q route -j oe -N ${target}${decimal}${first_part} -o ${outlog} -l ncpus=50,mem=100GB,walltime=50:00:00"
    #echo $qsub_cmd
    cmd="bash /mnt/s1/labs/brown/sag-mag/sag-mag_coverm.sh $target $decimal $mg1 $mg2"
    #echo $cmd
    echo $cmd | $qsub_cmd
    done
done
done



