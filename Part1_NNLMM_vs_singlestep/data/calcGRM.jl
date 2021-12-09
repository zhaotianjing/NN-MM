# create genomic relationship matrix based on the genotype
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,Random, LinearAlgebra, StatsBase
Random.seed!(1)

for rep in 1:20
    genofile="/Users/tianjing/Box/Omicsdata/andrew1055/data$rep/geno.csv"
    genotypes=get_genotypes(genofile;separator=',',method="RR-BLUP");

    sum2pq=genotypes.sum2pq

    geno_centered=genotypes.genotypes
    round.(mean(geno_centered,dims=1),digits=2)

    GRM=(geno_centered*geno_centered')/sum2pq
    writedlm("/Users/tianjing/Box/Omicsdata/andrew1055/data$rep/GRM_geno.txt",GRM)
end
