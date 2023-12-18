outdir=/mnt/stepanauskas_nfs/projects/sag_mag/coverm/manual_recruit2/
logdir=/mnt/s1/labs/brown/sag-mag/out_logs/

for reference in /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/merge*.fna; do
    for mg in /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/metag/1M_subsample/*interleaved*1Msub.fasta.gz; do
        file_name=$(basename "$mg")
        # Extract the part before the first underscore
        mg_id="${file_name%%_*}"
        
        refname=$(basename "$reference")
        refid="${refname%%.*}"
        
        outlog=${logdir}${mg_id}_${refid}_3.out
        
        
        cmd="bash /mnt/s1/labs/brown/sag-mag/scripts/recruitment_workflow.sh $reference $mg $outdir"
        qsub_cmd="qsub -q route -j oe -N ${mg_id}_${refid} -o ${outlog} -l ncpus=32,mem=100GB,walltime=50:00:00"
        
        echo $cmd | $qsub_cmd
        done
    done
    
