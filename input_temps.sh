#!/bin/bash

#bash script to run ising.py for different values of T

#This will allow use module ising_functions.py
#next line may cause error since the path wont be valid on other computer, but works on mine.
export PYTHONPATH="${PYTHONPATH}:/home/vlad/Documents/ising"

#temp_vals contains temperatures that will be piped into .py program as command line arguments
temp_vals=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2)

 
echo print n=system size, .txt file_name to store output data into, name of .py file_to_run, space separated

read N file_name file_to_run
#e.g. entering '10 output.txt ising.py' will run ising.py for different T for a 10*10 lattice and store output in output.txt


for num in "${temp_vals[@]}"
do
	echo $N $num|ipython $file_to_run>>$file_name 
	echo $N $num
done


#file_name should now contain T, M, E, X, C, S values at equilibrium as columns

#This should run plotting routine:
echo "$file_name"|ipython plot_data.py


