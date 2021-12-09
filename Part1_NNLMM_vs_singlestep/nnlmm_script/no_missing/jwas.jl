# this file is to run JWAS
using DataFrames,CSV,Random,Distributions,JWAS,DelimitedFiles, LinearAlgebra
Random.seed!(123)

method            = ARGS[1]
rep               = parse(Int, ARGS[2])  #data1,...,data20
@show method,rep

mainpath="/group/qtlchenggrp/tianjing/omics_fix_all/full_omics/"
cd(mainpath*"$method/data$rep/")

############ READ DATA ##########
#read phenotype
phenofile="/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$rep/y_no_ar.csv"
phenotypes = CSV.read(phenofile,DataFrame,delim = ',',header=true,missingstrings=["NA"])

#read omics
omicsfile  = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$rep/omics.csv"
omics  = CSV.read(omicsfile,DataFrame,delim = ',',header=true,missingstrings=["NA"]);
insertcols!(omics,2,:y => phenotypes[:,:y]); #insert phenotype into 2nd column of omics 
for i in 2:size(omics,2) #skip ID
    omics[!,i]= convert(Vector{Union{Missing,Float64}}, omics[!,i]);  #convert the data type to allow missing data
end
omics[956:1055,:y] .= missing  #missing phenotype for individuals in testing dataset
omics_name=["gene$i" for i in 1:1200]

@show size(omics)
@show omics[1:5,1:5];
@show omics[end-4:end,1:5];


#fixed variances components required in Christensenet al.(2021)
G=Matrix(Diagonal(repeat([1.22],1200))) #genetic variance for each omics is 1.22 in Christensenet al.(2021)
R=Matrix(Diagonal(repeat([0.78],1200))) #residual variance for omics is 0.78 in Christensenet al.(2021)
#read genotype
if occursin("gblup",method)  #NN-GBLUP
	genofile   = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$rep/GRM_geno.csv"
	genotypes  = get_genotypes(genofile,G;separator=',',method="GBLUP",header=false,G_is_marker_variance = false,estimateVariance=false);
else #NN-BayesC
	genofile   = "/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$rep/geno.csv"
	if occursin("bayescpi",method)
		genotypes  = get_genotypes(genofile,G;separator=',',method="BayesC",estimatePi=false,G_is_marker_variance = false,estimateVariance=false);
	end
end


# build model
model_equation = "y = intercept + genotypes";
#fixed variances components required in Christensenet al.(2021)                            
user_σ2_yobs=19.37777777777778  #residual variance for observed phenotype
user_σ2_weightsNN=0.01          #variance for neural network weight between omics and observed phenotype
if occursin("gblup",method)   #NN-GBLUP
	model = build_model(model_equation,R;
                        num_hidden_nodes=1200,
                        latent_traits=omics_name,
                        nonlinear_function="linear");
	out   = runMCMC(model,omics,chain_length=5000,printout_model_info=false,estimate_variance=false,user_σ2_yobs=user_σ2_yobs,user_σ2_weightsNN=user_σ2_weightsNN);
else  #NN-BayesC
	model = build_model(model_equation,R;
                        num_hidden_nodes=1200,
                        latent_traits=omics_name,
                        nonlinear_function="linear");
	out   = runMCMC(model,omics,chain_length=5000,printout_model_info=false,estimate_variance=false,user_σ2_yobs=user_σ2_yobs,user_σ2_weightsNN=user_σ2_weightsNN);
end


# Step 5: Check Accuruacy
#read true bv
bv=vec(readdlm("/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$rep/tbv_m.txt"));
accuruacy_all  = cor(out["EBV_NonLinear"][!,:EBV],bv)
accuruacy_test  = cor(out["EBV_NonLinear"][!,:EBV][956:1055],bv[956:1055])  #testing individuals

@show method,rep
@show accuruacy_all
@show accuruacy_test

res=[rep,accuruacy_test, accuruacy_all]
writedlm("res$rep.txt",res)





