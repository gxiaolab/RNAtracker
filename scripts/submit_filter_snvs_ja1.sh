#!/bin/bash

ja=filter_snvs.job # this is the file to be created

cell_lines=../input_files/ENCODE_cell_line_names  # each row is a different cell line

rm $ja

count_dat=../input_files/ALL_COUNT.dat
input_dir=../input_files

while read p; do
	cell=$p
	for stage in zero_hr two_hr six_hr
	do
		#echo $stage
		dir_ase=$input_dir/${cell}/${stage}

		pickle_file=$dir_ase/data_all.pickle
		files=$(find $dir_ase/*.rdd.final.all_snvs.EN*.rdd |  tr '\n' ' ') 
		count_dat=$count_dat
		echo $files $count_dat $dir_ase >> $ja
	done
done <$cell_lines



