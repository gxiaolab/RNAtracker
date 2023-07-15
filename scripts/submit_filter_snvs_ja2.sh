#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -t 1-48
#$ -l h_data=32G,h_rt=3:00:00
#$ -m n
#$ -j y
#$ -r n
#$ -o filter_snvs.out


ja=filter_snvs.job
PARMS=($(awk "NR==$SGE_TASK_ID" $ja))
rdd1=${PARMS[0]}
rdd2=${PARMS[1]}
count_dat=${PARMS[2]}
dir_ase=${PARMS[3]}


mkdir -p $dir_ase

### THESE LINES NEED TO BE MODIFIED ACCORDING TO YOUR CONDA SETUP
source /u/local/Modules/default/init/modules.sh
source /u/project/gxxiao4/huange7/miniconda3/etc/profile.d/conda.sh
###

conda activate RNAtracker
python RNAtracker_filter_snvs.py  --rdd ${rdd1} ${rdd2} -c $count_dat --prefix ${dir_ase} > log/$ja.$SGE_TASK_ID 2>&1
conda deactivate
