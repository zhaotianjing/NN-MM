#this file is to randomly selected some individuals to have no omics data.
using Distributions, DelimitedFiles,Random
Random.seed!(123)

cd("/Users/tianjing/Box/omics/omics_missing_pattern_legarra")

pct_all=[0.1,0.3,0.5,0.7,0.8,0.9,0.95,0.99]
repAll=20  #rep1,...,rep20; in each rep, randomly selected some individuals to have no omics data.

for i in 1:length(pct_all)
    pct=pct_all[i]
    nMissingInd=Int(floor(955*pct))
    println(pct*100," % training individuals (i.e., $nMissingInd training individuals) has no missing data.")

    for rep in 1:repAll
        missingIndIndex=sort(sample(1:955,nMissingInd,replace=false))
        println(length(missingIndIndex), " individuals are sampled in rep$rep.")
        writedlm("missing_index/missingIndex.pct$pct.rep$rep.txt",missingIndIndex)
    end
end
