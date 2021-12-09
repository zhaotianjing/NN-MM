#this file add the ID column into the genomic relationship matrix
setwd("/Users/tianjing/Box/Omicsdata/andrew1055")

for (rep in 1:20){
  GRM=as.matrix(read.table(paste0("data",rep,"/GRM_geno.txt"),header=F))
  cat(dim(GRM),"\n")
  ID=read.table(paste0("data",rep,"/select_1055_index_data",rep,".txt"))$V1
  GRM_id= cbind(ID,GRM)
  cat(dim(GRM_id),"\n")
  write.table(GRM_id,paste0("data",rep,"/GRM_geno.csv"),sep = ",",col.names=F,row.names=F)
}
