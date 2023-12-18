#!/usr/bin/env bash

# run coverm on "sag" or "mag" directory within /mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/ at an input threshold
# Usage: bash sag-mag_coverm.sh target_genomes identity_threshold

module purge
module load anaconda3
source activate /mnt/s1/labs/brown/envs/coverm

if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <sag_or_mag> <identity_threshold> <mg1> <mg2>"
  exit 1
fi


genome_target="$1"
pctid="$2"
mg1="$3"
mg2="$4"

cd /mnt/stepanauskas_nfs/projects/sag_mag/coverm/1Msub

#in_genomes=/mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/${genome_target}/
in_genomes=/mnt/scgc/stepanauskas_nfs/projects/sag_mag/seqs/merge_${genome_target}s.fna


file_name=$(basename "$mg2")

# Extract the part before the first underscore
first_part="${file_name%%_*}"

#outtsv=./${genome_target}_coverm_${pctid}_${first_part}.tsv
outtsv=./${genome_target}-concat_coverm-exclude_secondary_${pctid}_${first_part}.tsv

#if [ -e "$outtsv" ]; then
#echo "File $outtsv exists. exiting script."

#else

echo "Running coverm on $genome_target at $pctid percent identity in order to create $outtsv"



#cmd="coverm genome -1 $mg1 \
#-2 $mg2 \
#-d  ${in_genomes} \
#-p bwa-mem \
#--min-read-percent-identity $pctid \
#--min-covered-fraction 0 \
#--methods length count reads_per_base covered_bases covered_fraction \
#--output-file ${outtsv}"

cmd="coverm genome -1 $mg1 \
-2 $mg2 \
-f  ${in_genomes} \
-p bwa-mem \
--exclude-supplementary \
--min-read-percent-identity $pctid \
--min-covered-fraction 0 \
--methods length count reads_per_base covered_bases covered_fraction \
--output-file ${outtsv}"

echo $cmd

TMPDIR=./ $cmd
