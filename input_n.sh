#!/bin/bash

#will run input_temps.sh for different lattice sizes
#best to comment out final line of input_temps.sh before running this, to prevent automatic graph plotting


#Array of sizes to pipe to input_temps.sh
N_vals=(18 14 10 6)

#files to store output data into
file_names=(data18.txt data14.txt data10.txt data6.txt)

for ((i=0;i<${#N_vals[@]};++i));
do
	date +"%T"
	echo "${N_vals[i]} ${file_names[i]}"|./input_temps.sh
	echo "${N_vals[i]} ${file_names[i]}"
	date +"%T"
#print when program finishes working

done





