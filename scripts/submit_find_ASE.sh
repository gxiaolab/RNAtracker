#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_data=16G,h_rt=3:00:00 #,highp
#$ -m n
#$ -j y
#$ -r n
#$ -o find_ASE.out



### THESE LINES NEED TO BE MODIFIED ACCORDING TO YOUR CONDA SETUP
source /u/local/Modules/default/init/modules.sh
source /u/project/gxxiao4/huange7/miniconda3/etc/profile.d/conda.sh
###

conda activate RNAtracker
cell_lines=../input_files/ENCODE_cell_line_names  # each row is a different cell line
input_dir=../input_files

while read p; do
	cell=$p
for stage in zero_hr two_hr six_hr
do
	#echo $stage
	dir_ase=$input_dir/$cell/$stage
	out_dir=$input_dir/$cell/$stage
	python RNAtracker_ase.py --in-prefix ${dir_ase} --prefix ${out_dir} -p 0.05 --nsigma 1.95996 --filter-test-set 2 10 --exclude-intronic --adjustP fdr_bh  > ./log/ase_run_${cell}_${stage}_strict 2>&1
done
done <$cell_lines


conda deactivate
