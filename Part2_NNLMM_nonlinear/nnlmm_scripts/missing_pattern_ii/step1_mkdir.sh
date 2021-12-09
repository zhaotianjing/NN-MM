#!/bin/bash

#this file is to make folders
nRep=10  #rep1,...,rep10
methodAll=("omics_gblup")  
pctAll=("0.5" "0.7" "0.9") 

for pct in "${pctAll[@]}"
do
	mkdir -p "$pct"
	cd "$pct"

	for method in "${methodAll[@]}"
	do
	    mkdir -p "$method"
		  cd "$method"
	    for j in $( eval echo {1..$nRep} ) 
	    do
	      mkdir -p rep"$j"
	    done
		  cd ..

			mkdir -p run_"$method"
			mkdir -p accuracy_"$method"
	done
	cd ..
	
done

