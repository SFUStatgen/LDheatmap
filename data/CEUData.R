CEUData<-read.table("hapmapCEU.txt", header=TRUE)
CEUSNP <- CEUData[-1,]
for (LDindex in 1:length(CEUSNP))
     CEUSNP[,LDindex] <- genetics::as.genotype(CEUSNP[,LDindex])
CEUDist<- as.vector(as.matrix(CEUData[1,]), mode="numeric")
rm(LDindex, CEUData)
