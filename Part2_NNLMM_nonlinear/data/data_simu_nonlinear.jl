#this file simulate omics data for nonlinear system
# the QTLs ("kk"), QTL effects ("u.txt"), QTL position ("qtl.txt"), omics effects ("alpha.txt") in Christensenet al.(2021) are required
# above data can be found in Christensenet al.(2021) or below link:
# http://genoweb.toulouse.inra.fr/~alegarra/GOBLUP/

using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,Random, LinearAlgebra, StatsBase

cd("/Users/tianjing/Box/Omicsdata")

Random.seed!(123)
### step1. get genetic value of omics (since it needs to be transformed to get true bv)
#read QTL
nQTL=20_000 #there are 20000 QLTs in Christensenet al.(2021), but only 5000 will be used
inter=4     #select 1 every 4 QTLs (5000/20000 QTLs)
nOmics=1200
nInd=21100
σ2_mi=2
h2_m=0.61
σ2_gi=h2_m*σ2_mi      #genetic variance of each omics
σ2_ei=(1-h2_m)*σ2_mi  #residual variance of each omics

QTL=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/andrew_code/kk") #QTLs in Christensenet al.(2021)
QTL=QTL[:,2:end] #remove ID
pos=collect(inter:inter:nQTL)
QTL=QTL[:,pos]
writedlm("nonlinear/QTL5000.txt",QTL) #evenly select 5000 QTLs as in Christensenet al.(2021)

#read 500 selected QTL position for all 1200 omics features
QTL_select_pos=Int.(readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/andrew_code/qtl.txt"))
maximum(QTL_select_pos)
minimum(QTL_select_pos)
#read 500 selected QTL effects for all 1200 omics features
QTL_select_effect=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/andrew_code/u.txt")
#calculate omics
# 1st layer
G=zeros(nInd,nOmics) #genetic value of omics
E=zeros(nInd,nOmics) #residuals of omics
Random.seed!(1)
for i in 1:nOmics
	W=QTL[:,QTL_select_pos[:,i]]
	@show size(W),W[1:2,1:2]
    Gi=W*QTL_select_effect[:,i]
	Gi=Gi/std(Gi)*sqrt(σ2_gi)
	G[:,i]=Gi
	E[:,i]=randn(nInd)*sqrt(σ2_ei)
end

var(G,dims=1)
var(E,dims=1)

writedlm("nonlinear/varmi_2/omics_bv.txt",G)
writedlm("nonlinear/varmi_2/omics_residual.txt",E)

omics=G+E
writedlm("nonlinear/varmi_2/omics.txt",omics)
var(omics,dims=1)


mysigmoid(x) = 1/(1+exp(-x))
#save g(omics_g)
G_nonlinear=mysigmoid.(G)
writedlm("nonlinear/varmi_2/omics_bv_nonlinear.txt",G_nonlinear)
#save omics nonlinear
omics_nonlinear=mysigmoid.(omics)
writedlm("nonlinear/varmi_2/omics_nonlinear.txt",omics_nonlinear)


#read neural network weight between omics and phenotype (omics effects on phenotype)
w1=vec(readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/andrew_code/alpha.txt"))

#save true bv
tbv=G_nonlinear*w1
writedlm("nonlinear/varmi_2/tbv.txt",tbv)


#create y_nonlinear
omics_contribution=omics_nonlinear*w1
σ2_z=var(omics_contribution)
σ2_g_nonlinear=var(G_nonlinear*w1)
#h2=σ2_g_nonlinear/(σ2_z+σ2_e)
h2=0.337
σ2_e=σ2_g_nonlinear/h2 - σ2_z #residual variance of y

Random.seed!(1)
y_residual=randn(nInd)*sqrt(σ2_e)
var(y_residual)
y_nonlinear=omics_contribution+y_residual
var(y_nonlinear)
σ2_g_nonlinear/var(y_nonlinear)
writedlm("nonlinear/varmi_2/y.txt",y_nonlinear)

#move genotype from old to new folders
for i in 1:20
	mkdir("nonlinear/varmi_2/data$i");
end
#create 1055 dataset
y_nonlinear=vec(readdlm("nonlinear/varmi_2/y.txt"))
omics=readdlm("nonlinear/varmi_2/omics.txt")
tbv=vec(readdlm("nonlinear/varmi_2/tbv.txt"))
for i in 1:20
	id=Int.(vec(readdlm("andrew1055/data$i/select_1055_index_data$i.txt")))
	# save y_nonlinear
	yi=y_nonlinear[id]
	writedlm("nonlinear/varmi_2/data$i/y.txt",yi)
	yi_df=DataFrame(ID=id,y=yi)
	CSV.write("nonlinear/varmi_2/data$i/y.csv", yi_df)
	# save omics (should use omics, not g(omics) in the model)
	omicsi=omics[id,:]
	writedlm("nonlinear/varmi_2/data$i/omics.txt",omicsi)
	names=["gene$i" for i in 1:1200]
    omicsi_df=DataFrame(omicsi,names)
	insertcols!(omicsi_df, 1, :ID => id)
	CSV.write("nonlinear/varmi_2/data$i/omics.csv", omicsi_df)
	#save true bv
	bvi=tbv[id]
	writedlm("nonlinear/varmi_2/data$i/tbv.txt",bvi)
	bvi_df=DataFrame(ID=id,tbv=bvi)
	CSV.write("nonlinear/varmi_2/data$i/tbv.csv", bvi_df)
end

