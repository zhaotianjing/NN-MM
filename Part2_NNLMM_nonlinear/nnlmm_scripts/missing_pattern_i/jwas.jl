# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

method            = ARGS[1]
rep               = parse(Int, ARGS[2])  #rep1,...,rep10
pct               = parse(Float64, ARGS[3]) 

@show method,rep,pct

mainpath="/group/qtlchenggrp/tianjing/omics_sample_all/nonlinear/varmi_2/missing_pattern_legarra/"
cd(mainpath*"$pct/$method/rep$rep/")

############ READ DATA ##########
#read missing index
missingIndex=Int.(vec(readdlm(mainpath*"missing_index/missingIndex.pct$pct.rep$rep.txt")))
@show missingIndex[1:3]
@show length(missingIndex)

#read phenotype
phenofile="/group/qtlchenggrp/tianjing/omics_data/nonlinear/varmi_2/data1/y.csv"
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])

#read omics
omicsfile  = "/group/qtlchenggrp/tianjing/omics_data/nonlinear/varmi_2/data1/omics.csv"
omics  = CSV.read(omicsfile,DataFrame,delim = ',',header=true,missingstrings=["NA"]);
insertcols!(omics,2,:y => phenotypes[:,:y]);
for i in 2:size(omics,2) #skip ID
    omics[!,i]= convert(Vector{Union{Missing,Float64}}, omics[!,i]);
end
omics[956:1055,:y] .= missing  #testing dataset
omics[missingIndex,3:end].=missing; #set missing patterns
@show [count(ismissing,col) for col in eachcol(omics)]

omics_name=["gene$i" for i in 1:1200]

@show size(omics)
@show omics[1:5,1:5];
@show omics[end-4:end,1:5];

#read genotype 
if occursin("gblup",method)  #NN-GBLUP
	genofile   = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data1/GRM_geno.csv"
	genotypes  = get_genotypes(genofile;separator=',',method="GBLUP",header=false);
end

model_equation = "y = intercept + genotypes";  


model = build_model(model_equation;
                    num_hidden_nodes=1200,
                    latent_traits=omics_name,
                    nonlinear_function="sigmoid"); #"linear" for linear activation function
out   = runMCMC(model,omics,chain_length=5000,printout_model_info=false);


