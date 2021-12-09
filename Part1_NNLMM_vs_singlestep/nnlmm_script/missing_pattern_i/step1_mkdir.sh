#!/bin/bash

#this file is to make folders
nRep=20  #rep1,...,rep20
methodAll=("omics_gblup")  
pctAll=("0.1" "0.3" "0.5" "0.7" "0.8" "0.9" "0.95" "0.99") 

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

