# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

method            = ARGS[1]
rep               = parse(Int, ARGS[2])  #rep1,...,rep10
pct               = parse(Float64, ARGS[3]) 

@show method,rep,pct

mainpath="/group/qtlchenggrp/tianjing/omics_sample_all/nonlinear/varmi_2/random_missing/"
cd(mainpath*"$pct/$method/rep$rep/")

############ READ DATA ##########
#read phenotype
phenofile="/group/qtlchenggrp/tianjing/omics_data/nonlinear/varmi_2/data1/y.csv"
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])


#read omics
omicsfile  = "/group/qtlchenggrp/tianjing/omics_data/nonlinear/varmi_2/data1/partial_missing_data/omics_partialmissing.pct$pct.rep$rep.csv"
omics  = CSV.read(omicsfile,DataFrame,delim = ',',header=true,missingstrings=["NA"]);

insertcols!(omics,2,:y => phenotypes[:,:y]);
omics[!,:y]= convert(Vector{Union{Missing,Float64}}, omics[!,:y]);
omics[956:1055,:y] .= missing  #testing dataset

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

model = build_model(model_equation,
                    num_hidden_nodes=1200,
                    latent_traits=omics_name,
                    nonlinear_function="sigmoid");  #"linear" for linear activation function
out   = runMCMC(model,omics,chain_length=5000,printout_model_info=false);




