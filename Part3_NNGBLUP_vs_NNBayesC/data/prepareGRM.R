#this file is to add the ID column to GRM
setwd("/Users/tianjing/Box/Omicsdata/andrew_1000snp/data1")


GRM=as.matrix(read.table("GRM.txt"),header=F)
ID=read.table("/Users/tianjing/Box/Omicsdata/andrew1055/data1/select_1055_index_data1.txt")$V1
GRM_id= cbind(ID,GRM)
write.table(GRM_id,"GRM.csv",sep = ",",col.names=F,row.names=F)

