# this file is the Julia code for the system of two approaches in Christensenet al.(2021) when all individuals have all omics data
# the R code can be found in Christensenet al.(2021)
using JWAS,DataFrames,CSV,Statistics,JWAS.Datasets,DelimitedFiles,Random, LinearAlgebra, StatsBase
Random.seed!(1)

res_all=zeros(20)

for data in 1:20  #data1,...,data20
    n      = 1055
    nTrain = 955
    nTest  = n-nTrain
    nfeat  = 1200  #nOmics features

    σ2_α  = 0.01   #variance for omics effect α
    σ2_mi = 2      #variance for omics each feature mi
    h2    = 0.47   #total heritability in equation (1)
    h2_m  = 0.61   #omics heritability in equation (2)
    h2_ar = 0.2    #residual bv heritability, part of h2.

    nOmics = 1200

    # calculate variance component for data simulation
    σ2_omics_gi  = h2_m * σ2_mi
    σ2_omics_ei  = (1-h2_m)*σ2_mi
    σ2_am        = nOmics * σ2_omics_gi * σ2_α #total omics genetic variance (=var(Gα))
    σ2_y         = σ2_am/(h2-h2_ar)
    σ2_ar        = h2_ar * σ2_y
    σ2_ϵ         = σ2_y - nOmics*σ2_mi*σ2_α - σ2_ar

    zeta   = σ2_omics_ei/σ2_omics_gi
    eta1   = σ2_ϵ/σ2_α


    #read GRM
    GRM=readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data$data/GRM_geno.txt",'\t')
    size(GRM)
    GRM=GRM+I*0.00001
    Ginv=inv(GRM)

    #read omics data
    M = readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data$data/omics.txt",'\t')

    #read phenotype
    yy = vec(readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data$data/y_no_ar.txt",','))

    ### computing G.hat:  (5) in paper Christensenet al.(2021)
    LHS11 = ones(n)'ones(n)
    LHS12 = ones(n)'I
    LHS22 = I'I + zeta*Ginv

    LHS = [LHS11  LHS12
           LHS12' LHS22]
    RHS=[ones(n)'M
         I'M]
    sol_all=inv(LHS)*RHS
    G_hat = sol_all[2:end,:]
    size(G_hat)

    ### (4) in paper Christensenet al.(2021)
    LHS11 = ones(nTrain)'ones(nTrain)
    LHS12 = ones(nTrain)'M[1:nTrain,:]
    LHS22 = M[1:nTrain,:]' * M[1:nTrain,:] + eta1*I

    LHS = [LHS11  LHS12
           LHS12' LHS22]
    RHS=[ones(nTrain)'yy[1:nTrain]
         M[1:nTrain,:]'yy[1:nTrain]]

    SOL = inv(LHS)*RHS
    alpha=SOL[2:end]

    hat_a_m = G_hat*alpha



    ### prediction accuracy
    tbv=vec(readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data$data/tbv_m.txt"))

    @show cor(hat_a_m,tbv)
    @show cor(hat_a_m[956:1055],tbv[956:1055])

    println("----data$data")
    @show round(cor(hat_a_m,tbv),digits=3)
    @show round(cor(hat_a_m[956:1055],tbv[956:1055]),digits=3)
    res_all[data]=cor(hat_a_m[956:1055],tbv[956:1055])
end


res_all

writedlm("/Users/tianjing/Box/omics/legarra_fullomics_res_20data.txt",res_all)
