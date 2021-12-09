ENV["GKSwstype"]="100" #to enable saving figure in server

#this file is to calculate prediction accuracy for NN-LMM
using DataFrames,CSV,Random,DelimitedFiles,Statistics
using Plots

Random.seed!(123)
method="omics_gblup_1200"  #NN-GBLUP
main_path="/group/qtlchenggrp/tianjing/omics_fix_all/full_omics/"
cd(main_path*method)

nRep=20  #data1,...,data20
chain_length=5000
final_accuracy_all=ones(nRep,2)*999.0

for i in 1:nRep
    #read bv
    bv=vec(readdlm("/group/qtlchenggrp/tianjing/omics_data/andrew1055/data$i/tbv_m.txt"));

    #read results for repi
    if ispath("data$i/results/MCMC_samples_EBV_NonLinear.txt")
        mcmc=readdlm("data$i/results/MCMC_samples_EBV_NonLinear.txt",',')[2:end,:]
        if size(mcmc,1)==1000 #1000 iterations were saved
            mcmc_mean=zeros(size(mcmc))
            for i in 1:1000
                mcmc_mean[i,:]=mean(mcmc[1:i,:],dims=1)
            end
            accuracy_datai=zeros(1000)
            for i in 1:1000
                accuracy_datai[i]=cor(mcmc_mean[i,:][956:1055],bv[956:1055]) #testing individual
            end

            #save results
            open(main_path*"accuracy_$method/accuracy.$method.data$i.txt", "w") do io
                writedlm(io, accuracy_datai)
            end
            final_accuracy=round(accuracy_datai[end],digits=3)
            final_accuracy_all[i,:]=[i,accuracy_datai[end]]

            #plot results
            myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="full_omics.data$i.accu=$final_accuracy")
            savefig(myfig,main_path*"accuracy_$method/accuracy.$method.data$i.png")
        end
    end
end
open(main_path*"accuracy_$method/accuracy_final.$method.txt", "w") do io
    writedlm(io, final_accuracy_all)
end
