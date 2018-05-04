#Downloaded SNP data for 45 unrelated Chinese and 45 (actually really 44)
#unrelated Japanese from HapMap, Data Release #21/phaseII Jul06, base-pair 
#range same as in original example data for LDheatmap from Utah CEPH families.

ch.dat<-utils::read.table("hapmapCHB.txt",header=TRUE) #rows are SNPs!
names(ch.dat) #to see which columns give the people for a SNP (row)
snpnames.ch<-ch.dat[,"rsnum"] #names of snps
ch.snps<-t(ch.dat[,12:56]) #snp data - rows index people, cols index SNPs

jp.dat<-utils::read.table("hapmapJPT.txt", header=TRUE)
names(jp.dat) #to see which columns give the people for a SNP (row)
snpnames.jp<-jp.dat[,"rsnum"] #names of snps
jp.snps<-t(jp.dat[,12:56]) #snp data - rows index people, cols index SNPs

#Are the SNPs genotyped in the Chinese and Japanese the same?
all(snpnames.ch==snpnames.jp) #Yes, all are the same

#Make a variable that indicates chinese vs japanese
#Glue together the data for chinese and japanese
CHBJPTSNP<-data.frame(rbind(ch.snps,jp.snps))
names(CHBJPTSNP)<-c(as.character(snpnames.jp))

#Loop over the columns of CHBJPTSNP to replace the missing 
#genotypes  (coded as NN) and convert them to genotype objects
library(genetics)
for(i in 1:41) {
 tem<-as.character(CHBJPTSNP[,i]) #Currently columns are factors
 tem[tem=="NN"]<-""
 CHBJPTSNP[,i]<-genetics::genotype(tem,sep="")
}

#Remove any SNPs that are not sufficiently polymorphic in either 
#Chinese or Japanese
ind<-rep(NA,41)
for(i in 1:41) 
ind[i]<-summary(CHBJPTSNP[1:45,i])$Hu>0&summary(CHBJPTSNP[46:90,i])$Hu>0
CHBJPTSNP<-CHBJPTSNP[,ind] #only 13!

#Create a vector of base-pair positions for the 13 SNPs. Need to go 
#back to the ch.dat (or jp.dat) data frames to find the positions in
#a column named "pos".

CHBJPTDist<-ch.dat[,"pos"] 
CHBJPTDist<-CHBJPTDist[ind]

rm(ch.dat,ch.snps,i,ind,jp.dat,jp.snps,snpnames.ch,snpnames.jp,tem)
