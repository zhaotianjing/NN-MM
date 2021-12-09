# this file is the Julia code for the single-step approach in Christensenet al.(2021) when some individuals have no omics data
# the R code can be found in Christensenet al.(2021)
using DelimitedFiles,LinearAlgebra

pct=0.99  #also: 0.1,0.3,0.5,0.7,0.8,0.9,0.95,0.99
repAll=20 #rep1,...,rep20
accu_all=zeros(repAll)
for rep in 1:repAll
    #read missing index
    missing_index=Int.(vec(readdlm("/Users/tianjing/Box/omics/omics_missing_pattern_legarra/missing_index/missingIndex.pct$pct.rep$rep.txt")))
    nInd_noomics = length(missing_index)   #ind with no omics
    all_index=collect(1:1055)
    no_missing_index=setdiff(all_index,missing_index)
    new_order=[missing_index;no_missing_index]  
    println("missing pct $pct, rep $rep, #missing ind $nInd_noomics")

    ###read omics
    M_all=readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/omics.txt")
    #reorder
    M_all=M_all[new_order,:]

    ###read y
    yy = vec(readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/y_no_ar.txt"))  #phenotype
    #reorder
    yy=yy[new_order]

    ### read GRM
    GRM=readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/GRM_geno.txt")
    #reorder
    GRM=GRM[new_order,new_order]

    nfeat = 1200          # omics feature
    nInd=1055            # total ind
    #############
    #set missing patter (the last 10 testing inds must have omics data)
    #first nInd_noomics ind without omics data

    #############
    nInd_omics = nInd-nInd_noomics  # ind with omics
    nInd_train = 955       # first 90 ind for training
    nInd_test  =  nInd-nInd_train # last 10 ind for testing
    println("#train ind $nInd_train, #test ind $nInd_test")


    M = M_all[(nInd_noomics+1):end,:]  #remove omics data for first nInd_noomics indviduals
    size(M)

    #standardization
    sdM=vec(std(M,dims=1))
    Mmean = M-ones(size(M,1))*mean(M,dims=1)
    for i in 1:nfeat
        Mmean[:,i]=Mmean[:,i]/sdM[i]/sqrt(nfeat)
    end
    M=Mmean
    size(M)
    sum(round.(mean(M,dims=1),digits=2))
    sum(round.(var(M,dims=1),digits=10))
    M[1:5,1:5]

    MMpsinv = inv(M*M'+I*0.0001)
    MMpsinv[1:5,1:5]

    σ2_am = 2
    h2_m = 0.61
    h2 = 0.47
    h2_ar = 0.2

    c2m = (h2-h2_ar)/h2_m
    s2y = h2_m*σ2_am/(h2-h2_ar)
    σ2_ar = h2_ar*s2y
    σ2_ϵ =(1-c2m-h2_ar)*s2y

    eta1 = σ2_ϵ/(c2m*s2y)
    ## eta1 also equals σ2_ϵ/σ2_am


    zeta = (1-h2_m)/h2_m


    GGzeta = h2_m*GRM

    GGzeta = GGzeta + I*(1-h2_m)
    ## note: M is standardised such that sum over omics features of mean.diag equals 1.
    ## Therfore a priori expectation of MM' is h2_m*GG+(1-h2)*Identity .

    Omega_inv = inv(GGzeta)
    GGzeta22_inv = inv(GGzeta[(nInd_noomics+1):end,(nInd_noomics+1):end])
    Omega_inv[(nInd_noomics+1):end,(nInd_noomics+1):end] = Omega_inv[(nInd_noomics+1):end,(nInd_noomics+1):end] + MMpsinv - GGzeta22_inv

    Ginv=inv(GRM+I*0.0001)
    ## X contains separate intercept for those with/without omics


    RHS = [sum(yy[1:nInd_noomics]) ; sum(yy[(nInd_noomics+1):nInd_train]) ; yy[1:nInd_train] ; zeros(nInd_test)]


    LHS11 = Diagonal([nInd_noomics,nInd_train-nInd_noomics])


    LHS12 = [ [ones(nInd_noomics);zeros(nInd_omics)]'
              [zeros(nInd_noomics) ; ones(nInd_train-nInd_noomics) ; zeros(nInd_test)]' ]


    LHS22 = [     [I
                   zeros(nInd_test,nInd_train)]  zeros(nInd,nInd_test)] + eta1*Omega_inv


    LHS = Matrix([LHS11  LHS12
                  LHS12' LHS22])

    SOL1 = inv(LHS)*RHS   #u=M*alpha, pluggin u into next MME
    # write(SOL1, file="SOL1_mm")

    ## SOL1 is phenotype prediction from omics M%*%alpha

    ## now we aim for the breeding values


    RHS = [sum(SOL1[3:nInd+2]) ; SOL1[3:nInd+2]]

    LHS11 = nInd
    LHS12 = ones(nInd)'
    LHS22 = I + zeta*Ginv

    LHS = [ LHS11  LHS12
            LHS12' LHS22]

    SOL2 = inv(LHS)*RHS
    SOL2[1:5]

    hat_a_m = SOL2[2:(nInd+1)]


    # write.table(hat.a.m, file="EBVss2", quote=FALSE, row.names=FALSE)


    ##prediction accuracy
    tbv=vec(readdlm("/Users/tianjing/Box/Omicsdata/andrew1055/data1/tbv_m.txt"))
    tbv=tbv[new_order]
    accuruacy_all  = cor(hat_a_m,tbv)
    accuruacy_test  = cor(hat_a_m[956:1055],tbv[956:1055])

    println("----------")
    @show pct,rep
    @show round(accuruacy_all,digits=3)
    @show round(accuruacy_test,digits=3)
    accu_all[rep]=accuruacy_test
end

accu_all
round.(accu_all,digits=3)

writedlm("/Users/tianjing/Box/omics/omics_missing_pattern_legarra/legarra_res/res_legarra_ptc$pct.txt",accu_all)
