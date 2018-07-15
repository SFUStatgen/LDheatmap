# Functions below are used in conjunction with LDHeatmap.R
# Provide supporting functionality for major function

makeImageRect <- function(nrow, ncol, cols, name, byrow=TRUE) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  # Creates coordinate pairs, repeating either column numbers (if byrow = TRUE) or row numbers (if byrow = FALSE) to force that type of fill
  if(byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }
  rectGrob(x=right, y=top, 
           width=1/ncol, height=1/nrow, 
           just=c("right", "top"), 
           gp=gpar(col=NA, fill=cols),
           name=name)
}

makeImageText <- function(nrow, ncol, cols, name, flip = FALSE) {
  cols <- as.character(cols)
  cols[is.na(cols)] <- ""
  cols <- paste("   ", cols)
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  
  # Need to fill cells in different order, as was done to generate image
  if(flip){
    right <- rep(xx, each = nrow)
    top <- rep(yy, ncol)
  }
  else{
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  }
  textGrob(cols, x=right, y=top, 
           gp=gpar(cex=0.3),
           just=c("right", "top"), 
           name=name)
}

LDheatmapLegend.add <- function(color, LDmeasure, vp){
  ImageRect<- makeImageRect(2,length(color), cols=c(rep(NA,length(color)),color[length(color):1]),
                            "colorKey")
  keyVP <- viewport(x=1.1, y=-.10, height=.10, width=.5, just=c("right","bottom"), name="keyVP")
  #Adding the label 'Color key'
  if(LDmeasure=="r") {
    ttt<-expression(paste(R^2," Color Key"))
  } else {
    ttt<-"D' Color Key"
  }
  title<-textGrob(ttt, x=0.5, y=1.25, name="title", gp=gpar(cex=0.8))
  
  #Adding labels to the color key
  labels<-textGrob(paste(0.2*0:5), x=0.2*0:5,y=0.25, gp=gpar(cex=0.6), name="labels")
  
  #Drawing ticks at the bottom axis of the color key
  ticks<-segmentsGrob(x0=c(0:5)*0.2 , y0=rep(0.4,6), x1=c(0:5)*0.2 , y1=rep(0.5,6),name="ticks")
  
  #Drawing a box around the color key
  box <- linesGrob(x=c(0,0,1,1,0), y=c(0.5,1,1,0.5,0.5), name="box")
  
  key <- gTree(children=gList(ImageRect, title, labels, ticks, box), name = "Key", vp=keyVP)
  key
}

LDheatmapMapNew.add <- function(nsnps, add.map, genetic.distances, 
                             geneMapLocation=0.15,
                             geneMapLabelX=NULL, geneMapLabelY=NULL,
                             distances="physical", vp=NULL, 
                             SNP.name=NULL, ind=0, flip=FALSE){
  snp <- ((1:nsnps-1) + 0.5) / nsnps
  #####################
  if(add.map){
    min.dist <- min(genetic.distances) 
    max.dist <- max(genetic.distances)
    total.dist <- max.dist - min.dist
    
    if(flip) geneMapLocation<- (-geneMapLocation) # geneMapLocation gets flipped, reflects the gene bar
    
    # Drawing the diagonal line 
    seq.x <- c(0.5*geneMapLocation + 1/(nsnps*2),
               1+0.5*geneMapLocation - 1/(nsnps*2))
    seq.y <- c(-0.5*geneMapLocation + 1/(nsnps*2),
               1-0.5*geneMapLocation - 1/(nsnps*2))
    diagonal<-linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal", vp=vp) # we draw the line with linesGrob, based on geneMapLocation seq
    
    ## Adding line segments to the plot: (point1 <-> point2) 
    ## point1: relative position of a SNP on the scaled line
    ## point2: position of that SNP on the LD image  
    regionx <- seq.x[1] +
      ((genetic.distances-min.dist)/total.dist)*(seq.x[2]-seq.x[1])
    regiony <- seq.y[1] +
      ((genetic.distances-min.dist)/total.dist)*(seq.y[2]-seq.y[1]) 
    segments <- segmentsGrob(snp, snp, regionx, regiony, name="segments", vp=vp)
    
    ## Adding the text indicating Physical length of the region under study
    if (distances=="physical")
      mapLabel <- paste("Physical Length:", round((total.dist/1000),1),
                        "kb", sep="")
    else 
      mapLabel <- paste("Genetic Map Length:", round(total.dist,1),"cM",sep="")
    
    if (!flip) {
      if(is.null(geneMapLabelY)) geneMapLabelY <- 0.3
      if(is.null(geneMapLabelX)) geneMapLabelX <- 0.5
    }
    else {
      if(is.null(geneMapLabelY)) geneMapLabelY <- 0.8
      if(is.null(geneMapLabelX)) geneMapLabelX <- 0.4
    }
    title <- textGrob(mapLabel, geneMapLabelX, geneMapLabelY,
                      gp=gpar(cex=0.9), just="left", name="title")
    
    geneMap <- gTree(children=gList(diagonal, segments, title), name="geneMap")
    
    ## Labelling some SNPs 
    if (!is.null(SNP.name) && (any(ind!=0))){
      if (flip) {
        length_SNP_name <- max(nchar(SNP.name)) 
        long_SNP_name <- paste(rep(8,length_SNP_name), collapse="")
        name_gap <- convertWidth(grobWidth(textGrob(long_SNP_name)), "npc",valueOnly=TRUE)/sqrt(2)
        diagonal<-linesGrob(seq.x, seq.y, gp=gpar(lty=1), name="diagonal", vp=vp)
        #diagonal<-linesGrob(seq.x+name_gap, seq.y-name_gap, gp=gpar(lty=1), name="diagonal", vp=vp)
        segments <- segmentsGrob(snp, snp, regionx, regiony, name="segments", vp=vp)
        #segments <- segmentsGrob(snp+name_gap, snp-name_gap, regionx+name_gap, regiony-name_gap, name="segments", vp=vp)
        
        ############################################
        # Bug: symbols was set to NULL here for some reason
        symbols <- pointsGrob(snp[ind], snp[ind], pch="*",
                              gp=gpar(cex=1.25, bg="blue", col="blue"), name="symbols", vp=vp)
        ############################################
        # Figure out exact necessary coefficient for regionx and regiony with name_gap
        SNPnames <- textGrob(SNP.name, just="left", rot=-45,
                             regionx[ind]-sqrt(2+0.5)*name_gap, regiony[ind]+sqrt(2+0.5)*name_gap, gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
        # Think of better reason to use the +0.5
        # snp[ind], snp[ind], gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
        title <- editGrob(title, y=unit(geneMapLabelY+name_gap, "npc"))
        }
      else{
        symbols <- pointsGrob(snp[ind], snp[ind], pch="*",
                              gp=gpar(cex=1.25, bg="blue", col="blue"), name="symbols", vp=vp)
        SNPnames <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
                             regionx[ind], regiony[ind], gp=gpar(cex=0.6, col="blue"), name="SNPnames", vp=vp)
      }
        geneMap <- gTree(children=gList(diagonal, segments, title, symbols, SNPnames),name="geneMap")
    }} # if(add.map) end
  
  else if (!add.map && !is.null(SNP.name) && (any(ind!=0))){
    geneMap <- textGrob(paste(" ", SNP.name), just="left", rot=-45,
                        snp[ind], snp[ind], gp=gpar(cex=0.6, col="blue"),
                        name="SNPnames")
    if (flip) geneMap <- editGrob(geneMap, vp=vp)
  }
  else geneMap <- NULL
  
  geneMap
}

