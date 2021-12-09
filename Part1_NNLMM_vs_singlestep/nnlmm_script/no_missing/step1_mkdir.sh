#!/bin/bash

#this file is to make folders

nRep=20  #data1,...,data20
methodAll=("omics_gblup_1200")  #NN-GBLUP with 1200 intermediate omics features

for method in "${methodAll[@]}"
do
    mkdir -p "$method"
	  cd "$method"
    for j in $( eval echo {1..$nRep} ) 
    do
      mkdir -p data"$j"
    done
	  cd ..

	mkdir -p run_"$method"
	mkdir -p accuracy_"$method"
done

