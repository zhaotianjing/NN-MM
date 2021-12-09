#This file is to randomly select 5% individuals (i.e., 1055) from 21100 individuals in Christensenet al.(2021).
# genotypes("genotypes.txt"), omics("M.txt"), phenotypes("y.txt"), breeding values("ybv.txt") data in Christensenet al.(2021) were used
# above data can be found in Christensenet al.(2021) or below link:
# http://genoweb.toulouse.inra.fr/~alegarra/GOBLUP/

using DataFrames,CSV,Random,DelimitedFiles,Distributions
cd("/Users/tianjing/Box/Omicsdata/andrew1055")

Random.seed!(123)

nData=20  #create 20 replicates
#selet 5% animals per generation
for rep in 1:nData
    mkdir("data$rep")
    println("rep$rep------------")
    pct=0.05
    selectIndiv=sample(1:1100,Int(1100*pct),replace=false, ordered=true)  #1100 individuals in generation1 in Christensenet al.(2021)
    println("generation:",1," start:",1," ","end:",1100," ","all:",1100, " selected:",length(selectIndiv))
    for i in 2:11 #11 generation in total
        start_i = 1100+(i-2)*2000+1  #2000 individuals in generation1 in Christensenet al.(2021)
        end_i = start_i+2000-1
        selectIndiv_i = sample(start_i:end_i,Int(2000*pct),replace=false, ordered=true)
        selectIndiv = [selectIndiv ; selectIndiv_i]
        println("generation:",i," start:",start_i," ","end:",end_i," ","all:",end_i-start_i+1, " selected:",length(selectIndiv_i))

    end
    @show length(selectIndiv)
    println(selectIndiv[1:5])
    writedlm("data$rep/select_1055_index_data$rep.txt",selectIndiv)
end


#save y (phenotype), omics
y=vec(readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/y.txt"))  #phenotype data in Christensenet al.(2021)
M=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/M.txt")       #omics data in Christensenet al.(2021)
for rep in 1:nData #20 reps
    selectIndiv=Int.(vec(readdlm("data$rep/select_1055_index_data$rep.txt")))
    # save y
    println("number of individuals is ", length(y))
    y_select=y[selectIndiv]
    println("number of individuals after selection is ", length(y_select))
    writedlm("data$rep/y.txt",y_select)
    y_df=DataFrame(ID=selectIndiv,y=y_select)
    CSV.write("data$rep/y.csv", y_df)

    # save omics
    println("number of individuals is ", size(M,1), ", number of genes is ", size(M,2))
    M_select=M[selectIndiv,:]
    println("after selection, number of individuals is ", size(M_select,1), ", number of genes is ", size(M_select,2))
    writedlm("data$rep/omics.txt",M_select)
    names=["gene$i" for i in 1:size(M_select,2)]
    M_df=DataFrame(M_select,names)
    insertcols!(M_df, 1, :ID => selectIndiv)
    CSV.write("data$rep/omics.csv", M_df)
end



#remove ar from y 
# ar is polygenetic effect whose covariance matrix is defined by the pedigree or/and genotypes in Christensenet al.(2021)
# this part is ignored here for simplicity
using DataFrames,CSV,Random,DelimitedFiles,Distributions
cd("/Users/tianjing/Box/Omicsdata/andrew1055")
tbv=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/tbv.txt")  #breeding values in Christensenet al.(2021)
ar=tbv[:,2]
y=vec(readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/y.txt")) #phenotype data in Christensenet al.(2021)
y_no_ar= y-ar

for rep in 1:20
    println("rep$rep------------")
    selectIndiv=Int.(vec(readdlm("data$rep/select_1055_index_data$rep.txt")))
    # save y_no_ar
    y_no_ar_select = y_no_ar[selectIndiv]
    writedlm("data$rep/y_no_ar.txt",y_no_ar_select)
    y_no_ar_df=DataFrame(ID=selectIndiv, y=y_no_ar_select)
    CSV.write("data$rep/y_no_ar.csv", y_no_ar_df)
end


# save breeding values
tbv=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/tbv.txt") #breeding values in Christensenet al.(2021)
tbv_m=tbv[:,1]  #breeding values from omics (i.e., no ar)
for rep in 1:20
    println("rep$rep------------")
    selectIndiv=Int.(vec(readdlm("data$rep/select_1055_index_data$rep.txt")))
    bv_select = tbv_m[selectIndiv]
    writedlm("data$rep/tbv_m.txt",bv_select)
    bv_select_df=DataFrame(ID=selectIndiv, tbv_m=bv_select)
    CSV.write("data$rep/tbv_m.csv", bv_select_df)
end


#save genotype
using DataFrames,CSV,Random,DelimitedFiles,Distributions
cd("/Users/tianjing/Box/Omicsdata/andrew1055")

Random.seed!(123)

geno=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/genotypes.txt")[:,2] #genotypes data in Christensenet al.(2021)
genptype = zeros(length(geno),15000)
for i in 1:21100
    indi=parse.(Int32,split(geno[i],""))
    genptype[i,:] = indi
end


for rep in 1:20
    selectIndiv=Int.(vec(readdlm("data$rep/select_1055_index_data$rep.txt")))

    # save genotype
    println("number of individuals is ", size(genptype,1), ", number of SNP is ", size(genptype,2))
    genptype_select=genptype[selectIndiv,:]
    println("after selection, number of individuals is ", size(genptype_select,1), ", number of SNP is ", size(genptype_select,2))
    writedlm("data$rep/geno.txt",genptype_select)
    names=["snp$i" for i in 1:size(genptype_select,2)]
    genptype_df=DataFrame(genptype_select,names)
    insertcols!(genptype_df, 1, :ID => selectIndiv)
    CSV.write("data$rep/geno.csv", genptype_df)
end



#save tbv (=tbv_m+ar)
using DataFrames,CSV,Random,DelimitedFiles,Distributions
cd("/Users/tianjing/Box/Omicsdata/andrew1055")
tbv=readdlm("/Users/tianjing/Box/Omicsdata/andrewlegarra_simulated_omics/raw/tbv.txt")
tbv=tbv[:,3]
#save tbv_m
for rep in 1:20
    println("rep$rep------------")
    selectIndiv=Int.(vec(readdlm("data$rep/select_1055_index_data$rep.txt")))
    tbv_select = tbv[selectIndiv]
    writedlm("data$rep/tbv.txt",tbv_select)
    tbv_select_df=DataFrame(ID=selectIndiv, tbv=tbv_select)
    CSV.write("data$rep/tbv.csv", tbv_select_df)
end
