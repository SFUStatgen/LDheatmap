# Potential code to modularize

#################################################################################
# The code below is one of the first pieces of code in LDHeatmap
# Its main function is a distance calculation
# Could be easily rewritten as getDist(genetic.distances, LDmeasure){}

# ------------------------------
## Calculate or extract LDmatrix

if(inherits(gdat,"SnpMatrix")){
  ## Exclude SNPs with less than 2 alleles:
  # NOT YET IMPLEMENTED for snp.matrix
  #gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 2))
  #genetic.distances <- genetic.distances[gvars]
  #gdat <- gdat[gvars]
  
  ## Sort data in ascending order of SNPs map position:
  if(!is.vector(genetic.distances))
  {stop("Distance should be in the form of a vector")}
  o<-order(genetic.distances)
  genetic.distances<-genetic.distances[o]
  gdat<-gdat[,o]
  #myLD <- snpStats::ld(gdat,depth=ncol(gdat))
  if(LDmeasure=="r")
    LDmatrix <- snpStats::ld(gdat,depth=ncol(gdat)-1,stats="R.squared")
  else if (LDmeasure=="D")
    LDmatrix <- snpStats::ld(gdat,depth=ncol(gdat)-1,stats="D.prime")
  else 
    stop("Invalid LD measurement, choose r or D'.")      
  LDmatrix <- as.matrix(LDmatrix)
  LDmatrix[lower.tri(LDmatrix,diag=TRUE)] <- NA
  ## LDmatrix is upper-left-triangular, rather than the usual upper-right.
  #nsnp<-length(genetic.distances)
  #tem<-matrix(NA,nrow=nsnp,ncol=nsnp)
  #for(i in 1:(nsnp-1)) { tem[i,(i+1):nsnp]<-LDmatrix[i,1:(nsnp-i)] }
  #LDmatrix<-tem # need something faster than the for loop
  #row.names(LDmatrix)<-attr(myLD,"snp.names")
}
######################################################################


######################################################################
# The following is a for loop, meaning we can probably squeeze a better run time out of our package by replacing with an apply or some vector check
#------------------------
for(i in 1:ncol(gdat)) {
  if(!genetics::is.genotype(gdat[,i]))
    stop("column ",i," is not a genotype object\n")
}
######################################################################


######################################################################
# The following code is repeated nearly identically in the LDheatmap function, decent opportunity to modularize
# Repeated due to an if/else block
#--------------------------
o<-order(genetic.distances)
genetic.distances<-genetic.distances[o]
gdat<-gdat[,o]

######################################################################