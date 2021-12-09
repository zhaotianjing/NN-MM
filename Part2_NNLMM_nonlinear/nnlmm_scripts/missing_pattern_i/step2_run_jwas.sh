#! /bin/bash

nRep=10
methodAll=("omics_gblup")  
pctAll=("0.1" "0.3" "0.5" "0.7" "0.9")  


for pct in "${pctAll[@]}"
do
	cd "$pct"

	for method in "${methodAll[@]}"
	do
		cd run_"$method"

		for rep in $( eval echo {1..$nRep} )  #rep1,...,rep10
		do
		    sbatch /group/qtlchenggrp/tianjing/omics_sample_all/nonlinear/varmi_2/missing_pattern_legarra/jwas.sbatch $method $rep $pct
		done

		cd .. 
	done
	cd ..
done