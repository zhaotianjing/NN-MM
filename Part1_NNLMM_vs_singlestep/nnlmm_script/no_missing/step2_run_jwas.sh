#! /bin/bash

#this file is to run JWAS

nRep=20  #data1,...,data20
methodAll=( "omics_gblup_1200")  


for method in "${methodAll[@]}"
do
	cd run_"$method"

	for rep in $( eval echo {1..$nRep} )  #data1,...,data20
	do
	    sbatch /group/qtlchenggrp/tianjing/omics_fix_all/full_omics/jwas.sbatch $method $rep
	done

	cd .. 

done