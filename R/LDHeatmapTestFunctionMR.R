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
  
  # Colour selection scaling
  if(color[1]=="blueToRed") color = rainbow(20, start=4/6, end=0, s=.7)[20:1]
  if(newpage)grid.newpage()
  mybreak <- 0:length(color)/length(color)
  
  imgLDmatrix <- LDmatrix
  
  # Bug testing
  globalImgLDmatrix <<- 1 - imgLDmatrix
  #
  
  # Flip or not, determines way data is read into the display
  byrow<-ifelse(flip,FALSE,TRUE) #FALSE if flip=TRUE
  
  colcut <- as.character(cut(1-imgLDmatrix,mybreak,labels=as.character(color), include.lowest=TRUE))
  
  # Bug testing
  globalColcut <<- colcut
  #
  
  # Determines if colour is done as an integer or as a colour code, updates accordingly
  if(is.numeric(color)) colcut <- as.integer(colcut)
  ImageRect<-makeImageRect(dim(LDmatrix)[1],dim(LDmatrix)[2],colcut, name="heatmap",byrow)
  
  # Bug testing
  globalRect <<- ImageRect
  #
  
  # Controls text placement
  ImageText <- NULL
  if (text) ImageText<-makeImageText(dim(LDmatrix)[1],dim(LDmatrix)[2], round(imgLDmatrix, digits = 2), name="heatmaptext")
  title <- textGrob(title, 0.5, 1.05, gp=gpar(cex=1.0), name="title")

  if (flip) {
    ImageRect <- editGrob(ImageRect, vp=flipVP)
    ## Minor modification ##
    globalRotatedRect <<- ImageRect
    ## Minor modification ##
    if (text)
      # Added flip = TRUE parameter to better utilize makeImageText() in the flipped case
      ImageText <- makeImageText(dim(LDmatrix)[1],dim(LDmatrix)[2], round(imgLDmatrix, digits = 2), name="heatmaptext", flip = TRUE)
      textVal <<- ImageText
      ImageText <- editGrob(ImageText, vp=flipVP, rot=0, just=c("right", "top"))
  }
  
  # Updates heatmap in the gTree
  heatMap <- gTree(children=gList(ImageRect, ImageText, title), name="heatMap")
  
  # Draw a diagonal line indicating the physical or genetic map positions of the SNPs
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  ind <- match(SNP.name, row.names(LDmatrix), nomatch=0)
  geneMapVP <- NULL
  if (flip) geneMapVP <- flipVP
  geneMap <- LDheatmapMap.add (nsnps, genetic.distances=genetic.distances,
                               geneMapLocation=geneMapLocation,add.map,  
                               geneMapLabelX=geneMapLabelX,
                               geneMapLabelY=geneMapLabelY,
                               distances=distances, vp=geneMapVP, 
                               SNP.name=SNP.name, ind=ind, flip=flip)
  
  # Bug testing
  globalGenemap <<- geneMap
  #
  
  # Draw the Color Key
  if(add.key) Key <- LDheatmapLegend.add(color, LDmeasure, heatmapVP)
  else Key <- NULL
  
  # Assemble the heatmap, genetic map and color key into a grob and draw it
  LDheatmapGrob<-gTree(children=gList(heatMap, geneMap, Key),
                       vp=heatmapVP, name=name, cl="ldheatmap")
  grid.draw(LDheatmapGrob)
  if(pop){
    downViewport(heatmapVP$name)
    popViewport()} #pop the heat map viewport
  
  ldheatmap <- list(LDmatrix=LDmatrix, LDheatmapGrob=LDheatmapGrob, heatmapVP=heatmapVP, flipVP=geneMapVP,                
                    genetic.distances=genetic.distances, distances=distances, color=color)
  class(ldheatmap) <- "LDheatmap"
  invisible(ldheatmap)
}


######### MAIN ###########
# Library and base object
library(LDheatmap)
library(grid)
source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
library(rtracklayer)
library(snpStats)

data(GIMAP5.CEU)
# If error "could not find makeImageRect", run the LDheatmapHelpers.R script
ll <- LDTest(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)

# Viewports to explore, note the extra viewport that is not normally used in LDheatmap
testNorm <- viewport(width = unit(.8, "snpc"), height = unit(.8, "snpc"),
                     name="Normal")
testFlip <- viewport(width = unit(.8, "snpc"), height= unit(.8, "snpc"), y=0.57, angle=-45, name="flipVP") # y = 0.57 is measurement of triangle side with hypoteneuse 0.8
testFlipExtra <- viewport(width = unit(1.13, "snpc"), height = unit(.43, "snpc"), y = 0.91) # This extra VP might be the key to the added sections
legendVP <- viewport(x = 0.9, y = 0.1, height=.60, width=.6, just=c("right","bottom"), name="keyVP")
titleVP <- viewport(x = 0, y = 0, height =.10, width=.1, just=c("left", "bottom"), name="titleVP")

########### FOR THE FLIPPED IMAGE ############### 
# Need to properly scale the Key in the bottom right

# Organization of viewports, with layers
grid.newpage()
pushViewport(testNorm) # Layer 1, default viewport, same as used in LDheatmap function
  #grid.rect(gp = gpar())
  grid.draw(ll$LDheatmapGrob$children$heatMap) # Added to layer 1, is the heatmap as created by LDheatmap. Means it will be rotated as desired (or not)

pushViewport(testFlip) # Layer 2, viewport fits the rotated heatmap closely and can be used for desired modifications
  #grid.rect(gp = gpar())

upViewport() # Return to Layer 1 to reorient and draw genemap. Could also be done during the testNorm viewport step earlier
  # If we are adding stuff above the flipped heatmap, we will likely need to relocate the title for the genemap
  ll$LDheatmapGrob$children$geneMap$children$title$x <- unit(0.35, "snpc")
  ll$LDheatmapGrob$children$geneMap$children$title$y <- unit(0.71, "snpc")
  grid.draw(ll$LDheatmapGrob$children$geneMap) # Adds Genemap
  
pushViewport(legendVP)
  grid.draw(ll$LDheatmapGrob$children$Key) # Draws the key just off the bottom right of the viewport window created
  #grid.rect(gp = gpar())

# Probably not necessary or could be done better: relocates the title in the situation we have a flipped heatmap
upViewport()
pushViewport(titleVP)
  grid.draw(ll$LDheatmapGrob$children$heatMap$children$title)
##

upViewport()
pushViewport(testFlipExtra)
  #grid.rect(gp = gpar())

# Test plot just to verify what can be added above, need to work on ways to split this window based on number of calls
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

# Working with the rTrackLayer component
llgenesTest <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")

grid.newpage()
pushViewport(testNorm)
grid.rect(gp = gpar())
grid.draw(llgenesTest$LDheatmapGrob$children$transcripts)

grid.newpage()
pushViewport(testFlipExtra)
grid.rect(gp = gpar())
grid.draw(llgenesTest$LDheatmapGrob$children$transcripts)

# Lower y testflipextra
lowerTestFlipExtra <- viewport(width = unit(1.13, "snpc"), height = unit(0.43, "snpc"), y = 0.75)

grid.newpage()
pushViewport(lowerTestFlipExtra)
grid.rect(gp = gpar())
grid.draw(llgenesTest$LDheatmapGrob$children$transcripts)

# Reduced size testFlipExtra
#testFlipExtra <- viewport(width = unit(1.13, "snpc"), height = unit(.43, "snpc"), y = 0.91) # This extra VP might be the key to the added sections
perc  <- 0.25 # perc stands for percentage of window
smallTestFlipExtra <- viewport(width = unit(1.13, "snpc"), height = unit(.43*perc, "snpc"), y = 0.91 - 0.43*(1 - perc)/2)

grid.newpage()
pushViewport(smallTestFlipExtra)
grid.rect(gp = gpar())
grid.draw(llgenesTest$LDheatmapGrob$children$transcripts)

# Other testing with viewports that dont use unit: Makes no difference compared to previous unit(, "snpc") work
otherTestFlipExtra <- viewport(width = 1.13, height = 0.43, y = 0.5)

# Trying to modify the vp's internally: looking not very likely to work
testObj <- llgenesTest$LDheatmapGrob$children$transcripts
testObj$children$gene_plot_title$vp <- testFlipExtra
testObj$children$uc010lpu.1$vp <- testFlipExtra
testObj$children$uc003whr.1$vp <- testFlipExtra

grid.newpage()
pushViewport(testNorm)
grid.rect(gp = gpar())
pushViewport(testFlipExtra)
grid.rect(gp = gpar())
grid.draw(testObj)
upViewport()
grid.draw(testObj)
pushViewport(otherTestFlipExtra)
grid.rect(gp = gpar())
grid.draw(testObj)


# Attempting to modify the grob, not working too well
test2 <- editGrob(testObj, vp = otherTestFlipExtra)
test3 <- editGrob(testObj$children$uc010lpu.1, vp = otherTestFlipExtra)
test4 <- editGrob(testObj$children$uc003whr.1, vp = otherTestFlipExtra)
grid.newpage()
pushViewport(otherTestFlipExtra)
grid.rect(gp = gpar())
grid.draw(test2)
grid.draw(test3)
grid.draw(test4)

######### Expand this approach ###############
# It appears the images work best when displayed in unit(0.8, "snpc") square windows, experiment with moving the windows based on the number
# of tracks added. Something along the lines of a global counter that updates everytime a track is added (unique function call)
globalTrack <<- 0 # Might need to be added to the LDheatmap function instead as the extra functions all go through helpers
# In function calls, globalTrack <<- globalTrack + 1 
globalTrackTest <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")
if(exists("globalTrack")){
  globalTrack <<- globalTrack + 1
} 
# else create globalTrack
######## FOR THE NON-FLIPPED IMAGE ############
# Works by default, can duplicate using the same approach outlined above
llNoFlip <- LDTest(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=FALSE)

grid.newpage()
pushViewport(testNorm) # Layer 1, default viewport, same as used in LDheatmap function
  #grid.rect(gp = gpar())
  grid.draw(llNoFlip$LDheatmapGrob$children$heatMap)
  grid.draw(llNoFlip$LDheatmapGrob$children$geneMap)
  grid.draw(llNoFlip$LDheatmapGrob$children$Key)


  