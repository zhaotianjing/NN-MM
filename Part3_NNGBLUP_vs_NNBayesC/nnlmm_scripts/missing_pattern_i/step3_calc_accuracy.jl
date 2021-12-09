ENV["GKSwstype"]="100" #to enable saving figure in server

using DataFrames,CSV,Random,DelimitedFiles,Statistics
using Plots

Random.seed!(123)
methodALL=["omics_gblup","omics_bayescpi"]
main_path="/group/qtlchenggrp/tianjing/res_andrew_1000snp/missing_pattern_legarra/"
pctAll=["0.1","0.3","0.5","0.7","0.9"]
nRep=20
chain_length=5000
#read bv
bv=vec(readdlm("/group/qtlchenggrp/tianjing/omics_data/andrew_1000snp/data1/bv.txt"));

for method in methodALL
    for pct in pctAll
        cd(main_path*"$pct/$method")
        final_accuracy_all=ones(nRep,2)*999.0

        for i in 1:nRep
            #read results for repi
            if ispath("rep$i/results/MCMC_samples_EBV_NonLinear.txt")
                mcmc=readdlm("rep$i/results/MCMC_samples_EBV_NonLinear.txt",',')[2:end,:]
                if size(mcmc,1)==1000
                    mcmc_mean=zeros(size(mcmc))
                    for i in 1:1000
                        mcmc_mean[i,:]=mean(mcmc[1:i,:],dims=1)
                    end
                    accuracy_datai=zeros(1000)
                    for i in 1:1000
                        accuracy_datai[i]=cor(mcmc_mean[i,:][956:1055],bv[956:1055])
                    end

                    #save results
                    open(main_path*"$pct/accuracy_$method/accuracy.$method.rep$i.txt", "w") do io
                        writedlm(io, accuracy_datai)
                    end
                    final_accuracy=round(accuracy_datai[end],digits=3)
                    final_accuracy_all[i,:]=[i,accuracy_datai[end]]

                    #plot results
                    myfig=plot(collect(chain_length/1000:chain_length/1000:chain_length),accuracy_datai,xlabel="iteration",ylabel="accuracy",title="ind-level missing.pct$pct.rep$i.accu=$final_accuracy")
                    savefig(myfig,main_path*"$pct/accuracy_$method/accuracy.$method.pct$pct.rep$i.png")
                end
            end
        end
        open(main_path*"$pct/accuracy_$method/accuracy_final.$method.pct$pct.txt", "w") do io
            writedlm(io, final_accuracy_all)
        end
    end
end