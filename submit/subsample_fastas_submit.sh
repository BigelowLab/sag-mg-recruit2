

outdir=/mnt/stepanauskas_nfs/projects/sag_mag/seqs/metag/1M_subsample/

for mg1 in /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/metag/*paired.fasta.gz; do
    mg1base=$(basename "$mg1" | cut -d'_' -f1-2)
    mg1_out="${outdir}${mg1base}_1Msub.fasta"    
    cmd="module load anaconda3; source activate julia; seqtk sample -s 123 ${mg1} 1000000 > $mg1_out; gzip $mg1_out"
    qsub_cmd="qsub -q route -j oe -N ${mg1base} -o ${outlog} -l ncpus=1,mem=10GB,walltime=20:00"
    echo $cmd | $qsub_cmd
done

