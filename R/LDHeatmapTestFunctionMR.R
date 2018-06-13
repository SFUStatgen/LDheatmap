# File will be mostly copied LDheatmap code with my changes
# Looking to build the function again to improve functionality and my understanding of what does what

LDTest <- function(gdat, genetic.distances=NULL, 
                         distances="physical", LDmeasure="r", title="Pairwise LD",
                         add.map=TRUE, add.key=TRUE, geneMapLocation=0.15, 
                         geneMapLabelX=NULL, geneMapLabelY=NULL, 
                         SNP.name=NULL, color=NULL,
                         newpage=TRUE, name="ldheatmap", vp.name=NULL,
                         pop=FALSE, flip=NULL, text=FALSE){
  
  #_______________________Color Key__________________________________________##
  # Draw the Color Key
  if (is.null(color)) {
    if (inherits(gdat, "LDheatmap")) color <- gdat$color
    else  color <- grey.colors(20)
  }
  
  #_______________________Genetic Map________________________________________##
  # adds a genetic map to the heatmap plot along the diagonal
  # This part only identifies if the genemap will need to be flipped, does not do anything else yet
  if (is.null(flip)) {
    if (inherits(gdat, "LDheatmap") && !is.null(gdat$flipVP)) flip <- TRUE 
    else flip <- FALSE
  }
  
  #____________________________________________________________________________#
  ## If genetic.distances is missing, calculate an equispaced default:
  if(is.null(genetic.distances)) {
    if (inherits(gdat,"data.frame"))
      genetic.distances=1000*(1:ncol(gdat))
    else if(inherits(gdat,"matrix"))
      genetic.distances=1000*(1:length(gdat[1,]))
    else    # gdat is of class LDheatmap
      genetic.distances = gdat$genetic.distances
  }
  
  #____________________________________________________________________________#
  ## Calculate or extract LDmatrix, then stored in LDMatrix as an upper triangular matrix
  
  if(inherits(gdat,"SnpMatrix")){
    ## Exclude SNPs with less than 2 alleles:
    # NOT YET IMPLEMENTED for SnpMatrix, is implemented for data.frame, done in else structure below
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
  
  # Similar to the above but using a data.frame instead of a SnpMatrix
  else if(inherits(gdat,"data.frame")){
    for(i in 1:ncol(gdat)) {
      if(!genetics::is.genotype(gdat[,i]))
        stop("column ",i," is not a genotype object\n")
    }
    
    ## Exclude SNPs with less than 2 alleles:
    gvars <- unlist(sapply(gdat, function(x) genetics::nallele(x) == 2))
    genetic.distances <- genetic.distances[gvars]
    gdat <- gdat[gvars]
    
    ## Sort data in ascending order of SNPs map position:
    if(!is.vector(genetic.distances))
    {stop("Distance should be in the form of a vector")}
    o<-order(genetic.distances)
    genetic.distances<-genetic.distances[o]
    gdat<-gdat[,o]
    myLD <- genetics::LD(gdat)
    if(LDmeasure=="r")
      LDmatrix <- myLD[[LDmeasure]]^2   
    else if (LDmeasure=="D'")
      LDmatrix <- abs(myLD[[LDmeasure]])  
    else 
      stop("Invalid LD measurement, choose r or D'.")      
  }
  # Same as above but with LDheatmap structure
  else if(inherits(gdat,"LDheatmap")){
    LDmatrix <- gdat$LDmatrix
    distances <- gdat$distances
  }
  # Same as above but with matrix structure
  else if(inherits(gdat,"matrix")){
    if(nrow(gdat) != ncol(gdat))
      stop("The matrix of linkage disequilibrium measurements must be a square matrix")
    LDmatrix <- gdat
    LDmatrix[lower.tri(LDmatrix, diag=TRUE)] <- NA
  }
  else if(!missing(gdat))  
    stop(paste("No method for an object of class",class(gdat)))
  else
    stop("Need to supply LD matrix or genotypes")
  
  
  #____________________________________________________________________________#
  ## Draw the heat map
  heatmapVP <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
                        name=vp.name)
  flipVP <- viewport(width = unit(.8, "snpc"), height= unit(.8, "snpc"), y=0.6, angle=-45, name="flipVP")
  
}


testNorm <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
                     name="Normal")
testFlip <- viewport(width = unit(.8, "snpc"), height= unit(.8, "snpc"), y=0.57, angle=-45, name="flipVP") # y = 0.57 is measurement of triangle side with hypoteneuse 0.8
testFlipExtra <- viewport(width = unit(1.13, "snpc"), height = unit(.43, "snpc"), y = 0.91) # This extra VP might be the key to the added sections

testNorm
grid.newpage()
pushViewport(testNorm)
#grid.rect(gp = gpar())
grid.draw(ll$LDheatmapGrob$children$heatMap)
pushViewport(testFlip)
#grid.rect(gp = gpar())
upViewport()
grid.draw(globalGenemap)
pushViewport(testFlipExtra)
#grid.rect(gp = gpar())

tg <- textGrob("sample text")
rg <- rectGrob(width=1.1*grobWidth(tg),
                 height=1.3*grobHeight(tg))
boxedText <- gTree(children=gList(tg, rg))

label <- textGrob("A\nPlot\nLabel ",
                  x=0, just="left")
x <- seq(0.1, 0.9, length=50)
y <- runif(50, 0.1, 0.9)
gplot <-
  gTree(
    children=gList(rectGrob(gp=gpar(col="grey60",
                                    fill="white")),
                   linesGrob(x, y),
                   pointsGrob(x, y, pch=16,
                              size=unit(1.5, "mm"))),
    vp=viewport(width=unit(1, "npc") - unit(5, "mm"),
                height=unit(1, "npc") - unit(5, "mm")))
grid.draw(gplot)
