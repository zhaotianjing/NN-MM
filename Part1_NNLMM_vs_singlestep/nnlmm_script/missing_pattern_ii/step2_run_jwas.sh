#! /bin/bash

#this file is to run JWAS
nRep=20 #rep1,...,rep20
methodAll=("omics_gblup")  
pctAll=("0.1" "0.3" "0.5" "0.7" "0.8" "0.9" "0.95" "0.99") 


for pct in "${pctAll[@]}"
do
	cd "$pct"

	for method in "${methodAll[@]}"
	do
		cd run_"$method"

		for rep in $( eval echo {1..$nRep} )  #rep1,...,rep20
		do
		    sbatch /group/qtlchenggrp/tianjing/omics_fix_all/random_missing/jwas.sbatch $method $rep $pct
		done

		cd .. 
	done
	cd ..
done