#this file is to create omcis data for missing pattern (ii)
# for each omics feature, some random individuals have missing omics data
using Distributions, DelimitedFiles,Random,CSV,DataFrames
Random.seed!(123)

cd("/Users/tianjing/Box/Omicsdata/nonlinear/varmi_2/")

pct_all=[0.1,0.3,0.5,0.7,0.9]
repAll=20

omicsfile="data1/omics.csv"
omics_all= CSV.read(omicsfile,DataFrame,delim = ',',header=true,missingstrings=["NA"])
for i in 2:size(omics_all,2) #skip ID
    omics_all[!,i]= convert(Vector{Union{Missing,Float64}}, omics_all[!,i]);
end

@show size(omics_all)
@show omics_all[1:5,1:5];
@show omics_all[end-4:end,1:5];

Random.seed!(123)
for i in 1:length(pct_all)
    pct=pct_all[i]
    num_missing_per_gene=Int(floor(955*pct))  #1-955 are training individuals
    for rep in 1:repAll
        omics=copy(omics_all)
        for i in 2:1201
            omics[sample(1:955,num_missing_per_gene,replace=false),i].=missing;
        end
        @show [count(ismissing,col) for col in eachcol(omics)]
        CSV.write("data1/partial_missing_data/omics_partialmissing.pct$pct.rep$rep.csv",omics,missingstring="NA")
    end
end
