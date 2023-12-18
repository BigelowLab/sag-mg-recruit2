outdir=/mnt/stepanauskas_nfs/projects/sag_mag/seqs/metag/1M_subsample/

for mg1 in /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/metag/*1_paired.fasta.gz; do
    mg1base=$(basename "$mg1" | cut -d'_' -f1-2)
    mg1_out="${outdir}${mg1base}_1Msub.fasta"
    
    
    mg2="${mg1/1_paired.fasta.gz/2_paired.fasta.gz}"
    mg2base=$(basename "$mg2" | cut -d'_' -f1-2)
    mg2_out="${outdir}${mg2base}_1Msub.fasta"
    
    seqtk sample -s 123 ${mg1} 1000000 > $mg1_out
    gzip $mg1_out
    
    seqtk sample -s 123 ${mg2} 1000000 > $mg2_out
    gzip $mg2_out
done

# interleave subsampled fastas
for reads_1 in /mnt/stepanauskas_nfs/projects/sag_mag/seqs/metag/1M_subsample/*_1_1Msub.fasta.gz; do
    reads_2="${reads_1/_1_/_2_}"
    file_name=$(basename "$reads_2")
    # Extract the part before the first underscore
    mg_name="${file_name%%_*}"
    interleaved_out="${reads_1/_1_/_interleaved_}"
    
    seqtk mergepe $reads_1 $reads_2 | gzip > $interleaved_out
done

