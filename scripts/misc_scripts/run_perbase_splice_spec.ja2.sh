#!/bin/bash
#$ -V
#$ -t 1-24
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_rt=2:00:00,h_data=32G 
#$ -r n # job is NOT rerunable
#$ -m a # Email on abort
#$ -o run_perbase_splice_spec.out

. /u/local/Modules/default/init/modules.sh

ja=run_perbase_splice_spec.job 

PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
in_bam=${PARMS[0]}
out=${PARMS[1]}

log=log/perbase.$SGE_TASK_ID.log
exec 1>$log 2>&1

module load mamba
conda activate perbase

perbase base-depth $in_bam -b asrs_variants.bed -o $out -D 100000000 
