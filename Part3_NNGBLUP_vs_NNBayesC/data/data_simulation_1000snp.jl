#this file is to simulate data with 1000 SNPs
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,Random, LinearAlgebra, StatsBase
Random.seed!(1)
cd("/Users/tianjing/Box/Omicsdata/andrew_1000snp")

#read all 15000 SNP
snp_all=readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/geno.txt")  #genotype data in linear system
id=vec(Int.(readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/select_1055_index_data1.txt")))
selectall_id=Int.(collect(15000/1000:15000/1000:15000)) #evenly selected 1000 SNPs from 15000 SNPs
snp_1k=snp_all[:,selectall_id]
#save SNP1k
writedlm("data1/geno1k.txt",snp_1k)
snp_names=["snp$i" for i in 1:1000]
snp_1k_df=DataFrame(snp_1k,snp_names)
insertcols!(snp_1k_df, 1, :ID => id)
CSV.write("data1/geno1k.csv",snp_1k_df)



#select 250 QTL data set
nQTL_all=250
Random.seed!(2)
QTL_all_index=sort(sample(1:1000,nQTL_all,replace=false))
writedlm("data1/QTL250_index.txt",QTL_all_index)
QTL=snp_1k[:,QTL_all_index]
writedlm("data1/QTL250.txt",QTL)
QTL_df=snp_1k_df[:,[1;QTL_all_index.+1]]
CSV.write("data1/QTL.csv",QTL_df)


noQTL_index=setdiff(collect(1:1000),QTL_all_index)
geno_noqtl=snp_1k[:,noQTL_index]
writedlm("data1/geno750_noqtl.txt",geno_noqtl)
geno_noqtl_df=snp_1k_df[:,[1;noQTL_index.+1]]
CSV.write("data1/geno750_noqtl.csv",geno_noqtl_df)


nQTL_per_omics=Int(0.05*1000)
nOmics=20

nInd=size(QTL,1)
σ2_mi=2
h2_m=0.61
σ2_gi=σ2_mi*h2_m
σ2_ei=σ2_mi*(1-h2_m)

Random.seed!(3)
σ2_w1=0.01
w1=randn(nOmics)*sqrt(σ2_w1)
var(w1)
writedlm("data1/w1.txt",w1)


omics_names=["gene$i" for i in 1:nOmics]
h2=0.337
#h2=var_omics_g/(var_omics+σ2_e)
var_omics_g=nOmics*σ2_gi*σ2_w1
var_omics=nOmics*σ2_mi*σ2_w1
σ2_e=var_omics_g/h2-var_omics
@show var_omics_g/(σ2_e+var_omics)



#index of QTL for each omics
Random.seed!(4)
omics_qtl_index=Int.(zeros(nQTL_per_omics,nOmics))
for j in 1:nOmics
    omics_qtl_index[:,j]=sort(sample(1:250,nQTL_per_omics,replace=false))
end
writedlm("data1/omics_qtl_index.txt",omics_qtl_index)

#marker effects
Random.seed!(5)
w0=randn(nQTL_per_omics,nOmics)
writedlm("data1/w0.txt",w0)


#create omics
Random.seed!(6)
nInd=1055
Omics_G=zeros(nInd,nOmics)
Omics_E=zeros(nInd,nOmics)
for j in 1:nOmics
    QTL_j = QTL[:,omics_qtl_index[:,j]]
    w0_j  = w0[:,j]
    Omics_gi = QTL_j*w0_j
    Omics_gi = Omics_gi/std(Omics_gi)*sqrt(σ2_gi)
    Omics_G[:,j]=Omics_gi

    Omics_ei=randn(nInd)*sqrt(σ2_ei)
    Omics_E[:,j]=Omics_ei
end

var(Omics_G,dims=1)
var(Omics_E,dims=1)
writedlm("data1/omics_g.txt",Omics_G)
writedlm("data1/omics_e.txt",Omics_E)

Omics=Omics_G+Omics_E
var(Omics,dims=1)
writedlm("data1/omics.txt",Omics)
Omics_df=DataFrame(Omics,omics_names)
insertcols!(Omics_df,1,:ID => id)
CSV.write("data1/omics.csv",Omics_df)


y=Omics*w1+randn(nInd)*sqrt(σ2_e)
writedlm("data1/y.txt",y)
y_df=DataFrame(ID=id,y=y)
CSV.write("data1/y.csv",y_df)

bv=Omics_G*w1
writedlm("data1/bv.txt",bv)

@show var(bv)/var(y)


#make GRM
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,Random, LinearAlgebra, StatsBase
Random.seed!(1)

genofile="data1/geno1k.csv"
genotypes=get_genotypes(genofile;separator=',',method="RR-BLUP")

sum2pq=genotypes.sum2pq

geno_centered=genotypes.genotypes
round.(mean(geno_centered,dims=1),digits=2)

GRM=(geno_centered*geno_centered')/sum2pq
writedlm("data1/GRM.txt",GRM)
