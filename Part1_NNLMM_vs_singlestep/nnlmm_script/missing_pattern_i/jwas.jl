# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

method            = ARGS[1]
rep               = parse(Int, ARGS[2])  #rep1,...,rep20
pct               = parse(Float64, ARGS[3]) 

@show method,rep,pct

mainpath="/group/qtlchenggrp/tianjing/omics_fix_all/missing_pattern_legarra/"
cd(mainpath*"$pct/$method/rep$rep/")

############ READ DATA ##########
#read missing index
#some individuals have no omics data
missingIndex=Int.(vec(readdlm(mainpath*"missing_index/missingIndex.pct$pct.rep$rep.txt")))
@show missingIndex[1:3]
@show length(missingIndex)

#read phenotype
phenofile="/group/qtlchenggrp/tianjing/omics_data/andrew1055/data1/y_no_ar.csv"
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])

#read omics
omicsfile  = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data1/omics.csv"
omics  = CSV.read(omicsfile,DataFrame,delim = ',',header=true,missingstrings=["NA"]);
insertcols!(omics,2,:y => phenotypes[:,:y]); #insert phenotype into 2nd column of omics
for i in 2:size(omics,2) #skip ID
    omics[!,i]= convert(Vector{Union{Missing,Float64}}, omics[!,i]); #convert the data type to allow missing data
end
omics[956:1055,:y] .= missing  #missing phenotype for individuals in testing datasett
omics[missingIndex,3:end].=missing; #some individuals have no omics data
@show [count(ismissing,col) for col in eachcol(omics)]

omics_name=["gene$i" for i in 1:1200]

@show size(omics)
@show omics[1:5,1:5];
@show omics[end-4:end,1:5];

#fixed variances components required in Christensenet al.(2021)
if occursin("gblup",method) #NN-GBLUP
	G=Matrix(Diagonal(repeat([1.22],1200)))  #genetic variance for each omics is 1.22 in Christensenet al.(2021)
	R=Matrix(Diagonal(repeat([0.78],1200)))  #residual variance for omics is 0.78 in Christensenet al.(2021)
	genofile   = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data1/GRM_geno.csv"
	genotypes  = get_genotypes(genofile,G;separator=',',method="GBLUP",header=false,G_is_marker_variance = false,estimateVariance=false);
end

# build model
model_equation = "y = intercept + genotypes";  
#fixed variances components required in Christensenet al.(2021)   
user_σ2_yobs=19.37777777777778 #residual variance for observed phenotype
user_σ2_weightsNN=0.01         #variance for neural network weight between omics and observed phenotype

if occursin("gblup",method)   #NN-GBLUP
	model = build_model(model_equation,R;
                        num_hidden_nodes=1200,
                        latent_traits=omics_name,
                        nonlinear_function="linear");
	out   = runMCMC(model,omics,chain_length=5000,printout_model_info=false,user_σ2_yobs=user_σ2_yobs,user_σ2_weightsNN=user_σ2_weightsNN,estimate_variance=false);
end


# Step 5: Check Accuruacy
#read true bv
bv=vec(readdlm("/group/qtlchenggrp/tianjing/omics_data/andrew1055/data1/tbv_m.txt"));
accuruacy_all  = cor(out["EBV_NonLinear"][!,:EBV],bv)
accuruacy_test  = cor(out["EBV_NonLinear"][!,:EBV][956:1055],bv[956:1055])

@show pct,method,rep
@show accuruacy_all
@show accuruacy_test


res=[rep,accuruacy_test, accuruacy_all]
writedlm("res$rep.txt",res)





