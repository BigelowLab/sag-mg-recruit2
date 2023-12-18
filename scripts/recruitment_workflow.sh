#!/usr/bin/env bash

# run coverm on "sag" or "mag" directory within /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/ at an input threshold
# Usage: bash sag-mag_coverm.sh target_genomes identity_threshold

module purge
module load anaconda3
source activate /mnt/s1/labs/brown/envs/coverm

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <reference> <interleaved_mg> <outdir>"
  exit 1
fi


reference="$1"
mg="$2"
outdir="$3"

cores=32

refname=$(basename "$reference")
refid="${refname%%.*}"
file_name=$(basename "$mg")
mg_name="${file_name%%_*}"

# index reference file

if test -e "${reference}.bwt"; then
    echo "reference file already indexed"
else
    bwa index $reference
fi

# prepping reads and running bwa mem

result="${outdir}${refid}_vs_${mg_name}.bam"
indexed_result="${result/.bam/.bai}"

if test -e "${indexed_result}"; then
    echo "alignment file exists"
else
# full piped version to sorted bam file
bwa mem -t ${cores} ${reference} -p ${mg} | samtools view \
    -ShuF4q2 --threads ${cores} | samtools sort -m 8G --threads ${cores} -o ${result}       
samtools index $result $indexed_result
fi


python /mnt/s1/labs/brown/sag-mag/scripts/filter_bam.py count-pctid-range $result --minlen 0 --overlap 0


