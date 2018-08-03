# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


###########################################################################
#' @name LDheatmap
#' @aliases LDheatmap
#' @title This function produces a pairwise LD plot.
#' @description \code{LDheatmap()} is used to produce a graphical display, as a heat map,
#'of pairwise linkage disequilibrium (LD) measurements for SNPs.
#'The heat map is a false color image in the upper-left diagonal of a square plot.
#'Optionally, a line parallel to the diagonal of the image indicating
#'the physical or genetic map positions of the SNPs may be added, along
#'with text reporting the total length of the genomic region considered.
#'
#' @usage LDheatmap(gdat, genetic.distances=NULL, distances="physical",
#'LDmeasure="r", title="Pairwise LD", add.map=TRUE, add.key=TRUE,
#'geneMapLocation=0.15, geneMapLabelX=NULL, geneMapLabelY=NULL,
#'SNP.name=NULL, color=NULL, newpage=TRUE,
#'name="ldheatmap", vp.name=NULL, pop=FALSE, flip=NULL, text=FALSE)
#'
#' @param gdat SNP data: a data frame of genotype objects, a \code{SnpMatrix} object, a square matrix of
#'pairwise linkage disequilibrium measurements or an object of
#'class \code{"LDheatmap"} (the returned object of this function).
#'
#' @param genetic.distances A numeric vector of map locations of the SNPs,
#'in the same order as SNPs listed in \code{gdat},
#'in terms of genetic or physical distances.
#'Physical distances should be in bases, genetic distances should be in
#'centiMorgans (cM).
#' When \code{gdat} is \emph{not} an object of class \code{LDheatmap}, the
#'default is a vector that represents equi-spaced markers, 1kb (1000 bases)
#'apart. When \code{gdat} \emph{is} an object of class \code{LDheatmap},
#'the \code{genetic.distances} argument is taken to be the
#'\code{genetic.distances} list item of \code{gdat}.
#'
#'@param distances A character string to specify whether the provided map locations
#'are in physical or genetic distances.
#'If \code{distances="physical"} (default), the text describing the total
#'length of the region will be \dQuote{Physical Length:XXkb} where XX is the
#'length of the region in kilobases. If \code{distances="genetic"}, the
#'text will be \dQuote{Genetic Map Length:YYcM} where YY is
#'the length of the region in centiMorgans.
#'If \code{gdat} is an object of class LDheatmap,
#'\code{distances} is taken from \code{gdat}.
#'
#'@param LDmeasure A character string specifying the measure of LD
#'- either allelic correlation \eqn{r^2} or Lewontin's
#'|D\eqn{'}|; default = \code{"r"} for \eqn{r^2};
#'type \code{"D'"} for |D\eqn{'}|.  This argument is ignored when the user has already
#'supplied calculated LD measurements through \code{gdat} (i.e., when \code{gdat}
#'is a matrix of pairwise LD measurements or an object of class \code{"LDheatmap"}).
#'
#'@param title A character string for the main title of the plot.
#'Default is \dQuote{Pairwise LD}.
#'
#'@param add.map If \code{TRUE} (default) a diagonal line indicating
#'the physical or genetic map positions of the SNPs will be added to 
#'the plot, along with text indicating the total length of the 
#'genetic region.
#'@param add.key If \code{TRUE} (default) the color legend is drawn.
#'@param geneMapLocation A numeric value specifying the position of the line
#'parallel to the diagonal of the matrix; the larger the value, the
#'farther it lies from the matrix diagonal. Ignored when \code{add.map=FALSE}.
#'
#'@param geneMapLabelX A numeric value specifying the x-coordinate
#'of the text indicating the total length of the genomic region
#'being considered. Ignored when \code{add.map=FALSE}.
#'
#'@param geneMapLabelY A numeric value specifying the y-coordinate
#'of the text indicating the total length of the genomic region
#'being considered. Ignored when \code{add.map=FALSE}.
#'
#'@param SNP.name A vector of character string(s) of SNP name(s) to
#'be labelled. Should match the names of SNPs in the provided object
#'\code{gdat}, otherwise nothing is done.
#'@param color A range of colors to be used for drawing the heat map. Default
#is \code{grey.colors(20)}.
#'@param newpage If \code{TRUE} (default), the heat map will be drawn on a new page.
#'@param name A character string specifying the name of the LDheatmap
#'graphical object (\code{grob}) to be produced.
#'@param vp.name A character string specifying the name of the viewport
#'where the heat map is going to be drawn.
#'@param pop If \code{TRUE}, the viewport where the heat map is drawn is
#'\code{pop}ped (i.e. removed) from the viewport tree after drawing. Default=\code{FALSE}.
#'
#'@param flip If \code{TRUE}, the LDheatmap plot is flipped below a horizontal line, in the style of Haploview. Default is \code{FALSE}.
#'@param text If \code{TRUE}, the LD measurements are printed on each cell.
#'
#'@details
#'The input object \code{gdat} can be a data frame of \code{genotype} objects
#'(a data structure from the \pkg{genetics} package), a \code{SnpMatrix} object (a
#'                                                                               data structure from the \pkg{snpStats} package), or
#'any square matrix with values between 0 and 1
#'inclusive.
#'LD computation is much faster for \code{SnpMatrix} objects than for
#'\code{genotype} objects.
#'In the case of a matrix of LD values between 0 and 1,
#'the values above the diagonal will be plotted.
#'In the display of LD, SNPs appear in the order supplied by the
#'user as the horizontal and vertical coordinates are increased
#'and one moves along the off-diagonal line, from the bottom-left
#'to the top-right corner.  To achieve this, the conventions of
#'the \code{image()} function have been adopted, in which horizontal
#'coordinates correspond to the rows of the matrix and vertical coordinates
#'correspond to columns, and vertical coordinates are indexed in increasing
#'order from bottom to top.
#'For the argument \code{color}, an appropriate
#'color palette for quantitative data is recommended,
#'as outlined in the help page of
#'the \code{\link[RColorBrewer:ColorBrewer]{brewer.pal}()} function of
#'the
#'\pkg{RColorBrewer} package.
#'See the package vignette \code{LDheatmap} for more examples and details
#'of the implementation. Examples of adding ``tracks'' of genomic
#'annotation above a flipped heatmap are in the package vignette
#'\code{addTracks}.
#'
#'
#'@return An object of class \code{"LDheatmap"} which contains the following components:
#' \item{LDmatrix}{ The matrix of pairwise LD measurements plotted in the heat map. }
#' \item{LDheatmapGrob}{ A grid graphical object (grob) representing the produced heat map. }
#' \item{heatmapVP}{ The viewport in which the heat map is drawn. See \link[grid:viewport]{viewport}.}
#' \item{genetic.distances}{The vector of the supplied physical or
#'genetic map locations, or the vector of equispaced marker distances
#'when no distance vector is supplied.}
#' \item{distances}{ A character string specifying whether the provided map
#'distances are physical or genetic. }
#' \item{color}{ The range of colors used for drawing the heat map. }
#' The \code{grob} \code{LDheatmapGrob} has three \code{grob}s as its children (components).
#'They are listed below along with their own children and respectively represent
#'the color image with main title, genetic map and color key of the heat map:
#'  \code{"heatMap"} - \code{"heatmap"}, \code{"title"};
#'\code{"geneMap"} - \code{"diagonal"}, \code{"segments"},
#'\code{"title"}, \code{"symbols"}, \code{"SNPnames"}; and
#'\code{"Key"} - \code{"colorKey"}, \code{"title"}, \code{"labels"},
#'\code{"ticks"}, \code{"box"}.
#'
#'@references Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap:
#'An R Function for Graphical Display of Pairwise Linkage
#'Disequilibria Between Single Nucleotide Polymorphisms.
#'Journal of Statistical Software, \bold{16} Code Snippet 3
#'
#'@note The produced heat map can be modified in two ways.
#'First, it is possible to edit \emph{interactively} the grob components of the heat map,
#'by using the function \code{\link[grid:grid.edit]{grid.edit}};
#'the function will not work if there is no
#'open graphical device showing the heat map.
#'Alternatively, the user can use the function
#'\code{\link[grid:grid.edit]{editGrob}} and work with
#'the grob \code{LDheatmapGrob} returned by \code{LDheatmap}.
#'See Examples for usage.

#'\code{LDheatmap()} uses \code{\link[grid:Grid]{Grid}}, which
#'does not respond to \code{par()} settings.
#'Hence modifying \code{par()} settings of \code{mfrow} and \code{mfcol}
#'will not work with \code{LDheatmap()}. The Examples section shows how to
#'display multiple heat maps on one plot without the use
#'of \code{par()}.
#'
#'
#'@author Ji-hyung Shin <shin@sfu.ca>, Sigal Blay <sblay@sfu.ca>,  Nicholas
#'Lewin-Koh <nikko@hailmail.net>, Brad McNeney <mcneney@stat.sfu.ca>, Jinko
#'Graham <jgraham@cs.sfu.ca>
#'
#'@seealso  \code{\link[genetics]{LD}}, \code{\link[genetics]{genotype}},
#'\code{\link[grid:Grid]{Grid}}, \code{\link{LDheatmap.highlight}},
#'\code{\link{LDheatmap.marks}}
#'
#'@examples # Pass LDheatmap a SnpMatrix object
#'set.seed(1)
#'#make an example matrix of genotypes, coded as 0, 1 2 copies of an index allele
#'gdat<-matrix(rbinom(n=500,size=2,prob=.5),ncol=5)
#'require(snpStats)
#'gdat<-as(gdat,"SnpMatrix")
#'LDheatmap(gdat,genetic.distances=c(0,1000,3000,4000,10000))
#'#Load the package's data set
#'data(CEUData)
#'#Creates a data frame "CEUSNP" of genotype data and a vector "CEUDist"
#'#of physical locations of the SNPs
#'# Produce a heat map in a grey color scheme

#'MyHeatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist,
#'                       color = grey.colors(20))

#'# Same heatmap, flipped below a horizontal gene map -- for examples of
#'# adding genomic annotation tracks to a flipped heatmap see
#'# vignette("addTracks")

#'# flippedHeatmap<-LDheatmap(MyHeatmap,flip=TRUE)

#'# Prompt the user before starting a new page of graphics output
#'# and save the original prompt settings in old.prompt.
#'old.prompt <- devAskNewPage(ask = TRUE)

#'# Highlight a certain LD block of interest:
#'LDheatmap.highlight(MyHeatmap, i = 3, j = 8, col = "black", 
#'fill = "grey",flipOutline=FALSE, crissCross=FALSE)
#'# Plot a symbol in the center of the pixel which represents LD between
#'# the fourth and seventh SNPs:
#'LDheatmap.marks(MyHeatmap,  4,  7,  gp=grid::gpar(cex=2),  pch = "*")


#'#### Use an RGB pallete for the color scheme ####
#'rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")
#'LDheatmap(MyHeatmap, color=rgb.palette(18))


#'#### Modify the plot by using 'grid.edit' function ####
#'#Draw a heat map where the SNPs "rs2283092" and "rs6979287" are labelled.
#'require(grid)
#'LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"))

#'#Find the names of the top-level graphical objects (grobs) on the current display
#'getNames()
#'#[1] "ldheatmap"

#'# Find the names of the component grobs of "ldheatmap"
#'childNames(grid.get("ldheatmap"))
#'#[1] "heatMap" "geneMap" "Key"

#'#Find the names of the component grobs of heatMap
#'childNames(grid.get("heatMap"))
#'#[1] "heatmap" "title"

#'#Find the names of the component grobs of geneMap
#'childNames(grid.get("geneMap"))
#'#[1] "diagonal" "segments" "title"    "symbols"  "SNPnames"

#'#Find the names of the component grobs of Key
#'childNames(grid.get("Key"))
#'#[1] "colorKey" "title"    "labels"   "ticks"    "box"

#'#Change the plotting symbols that identify SNPs rs2283092 and rs6979287
#'#on the plot to bullets
#'grid.edit("symbols", pch = 20, gp = gpar(cex = 1))

#'#Change the color of the main title
#'grid.edit(gPath("ldheatmap", "heatMap", "title"), gp = gpar(col = "red"))

#'#Change size of SNP labels
#'grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=1.5))

#'#Add a grid of white lines to the plot to separate pairwise LD measures
#'grid.edit(gPath("ldheatmap", "heatMap", "heatmap"), gp = gpar(col = "white",
#'                                                              lwd = 2))


#'#### Modify a heat map using 'editGrob' function ####
#'MyHeatmap <- LDheatmap(MyHeatmap, color = grey.colors(20))

#'new.grob <- editGrob(MyHeatmap$LDheatmapGrob, gPath("geneMap", "segments"),
#'                     gp=gpar(col="orange"))

#'##Clear the old graphics object from the display before drawing the modified heat map:
#'grid.newpage()

#'grid.draw(new.grob)
#'# now the colour of line segments connecting the SNP
#'# positions to the LD heat map has been changed from black to orange.


#'#### Draw a resized heat map (in a 'blue-to-red' color scale ####
#'grid.newpage()

#'pushViewport(viewport(width=0.5, height=0.5))
#'LDheatmap(MyHeatmap, SNP.name = c("rs2283092", "rs6979287"), newpage=FALSE,
#'          color="blueToRed")
#'popViewport()


#'#### Draw and modify two heat maps on one plot ####
#'grid.newpage()

#'##Draw and the first heat map on the left half of the graphics device
#'pushViewport(viewport(x=0, width=0.5, just="left"))
#'LD1<-LDheatmap(MyHeatmap, color=grey.colors(20), newpage=FALSE,
#'               title="Pairwise LD in grey.colors(20)",
#'               SNP.name="rs6979572", geneMapLabelX=0.6,
#'               geneMapLabelY=0.4, name="ld1")
#'upViewport()

#'##Draw the second heat map on the right half of the graphics device
#'pushViewport(viewport(x=1,width=0.5,just="right"))
#'LD2<-LDheatmap(MyHeatmap, newpage=FALSE, title="Pairwise LD in heat.colors(20)",
#'               SNP.name="rs6979572", geneMapLabelX=0.6, geneMapLabelY=0.4, name="ld2")
#'upViewport()

#'##Modify the text size of main title of the first heat map.
#'grid.edit(gPath("ld1", "heatMap","title"), gp=gpar(cex=1.5))

#'##Modify the text size and color of the SNP label of the second heat map.
#'grid.edit(gPath("ld2", "geneMap","SNPnames"), gp=gpar(cex=1.5, col="DarkRed"))

#'#### Draw a lattice-like plot with heat maps in panels ####
#'# Load CHBJPTSNP and CHBJPTDist
#'data(CHBJPTData)
#'# Make a variable which indicates Chinese vs. Japanese
#'pop <- factor(c(rep("chinese",45), rep("japanese",45)))
#'require(lattice)

#'xyplot(1:nrow(CHBJPTSNP) ~ 1:nrow(CHBJPTSNP) | pop,
#'       type="n", scales=list(draw=FALSE), xlab="", ylab="",
#'       panel=function(x, y, subscripts,...) {
#'         LDheatmap(CHBJPTSNP[subscripts,], CHBJPTDist, newpage=FALSE) })

#'data(GIMAP5)
#'require(lattice)
#'n<-nrow(GIMAP5$snp.data)
#'xyplot(1:n ~ 1:n | GIMAP5$subject.support$pop,
#'       type="n", scales=list(draw=FALSE), xlab="", ylab="",
#'       panel=function(x, y, subscripts,...) {
#'         LDheatmap(GIMAP5$snp.data[subscripts,],
#'                   GIMAP5$snp.support$Position, SNP.name="rs6598", newpage=FALSE) })



#'#Reset the user's setting for prompting on the graphics output
#'#to the original value before running these example commands.
#'devAskNewPage(old.prompt)
#'
#'@keywords hplot
#' @export

"LDheatmap"<- function (gdat, genetic.distances=NULL, 
             distances="physical", LDmeasure="r", title="Pairwise LD",
             add.map=TRUE, add.key=TRUE, geneMapLocation=0.15, 
             geneMapLabelX=NULL, geneMapLabelY=NULL, 
             SNP.name=NULL, color=NULL,
             newpage=TRUE, name="ldheatmap", vp.name=NULL,
             pop=FALSE, flip=NULL, text=FALSE)
{
  
  requireNamespace("grid")
  
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

  # Flip or not, determines way data is read into the display
  byrow<-ifelse(flip,FALSE,TRUE) #FALSE if flip=TRUE
  
  colcut <- as.character(cut(1-imgLDmatrix,mybreak,labels=as.character(color), include.lowest=TRUE))
  
  
  # Determines if colour is done as an integer or as a colour code, updates accordingly
  if(is.numeric(color)) colcut <- as.integer(colcut)
  ImageRect<-makeImageRect(dim(LDmatrix)[1],dim(LDmatrix)[2],colcut, name="heatmap",byrow)
  
  
  # Controls text placement
  ImageText <- NULL
  if (text) ImageText<-makeImageText(dim(LDmatrix)[1],dim(LDmatrix)[2], round(imgLDmatrix, digits = 2), name="heatmaptext")
  title <- textGrob(title, 0.5, 1.05, gp=gpar(cex=1.0), name="title")
  
  if (flip) {
    ImageRect <- editGrob(ImageRect, vp=flipVP)
    if (text){
      # Added flip = TRUE parameter to better utilize makeImageText() in the flipped case
      ImageText <- makeImageText(dim(LDmatrix)[1],dim(LDmatrix)[2], round(imgLDmatrix, digits = 2), name="heatmaptext", flip = TRUE)
      textVal <- ImageText
      ImageText <- editGrob(ImageText, vp=flipVP, rot=0, just=c("right", "top"))
    }
  }
  
  # Updates heatmap in the gTree
  heatMap <- gTree(children=gList(ImageRect, ImageText, title), name="heatMap")
  
  # Draw a diagonal line indicating the physical or genetic map positions of the SNPs
  nsnps <- ncol(LDmatrix)
  step <- 1/(nsnps-1)
  ind <- match(SNP.name, row.names(LDmatrix), nomatch=0)
  geneMapVP <- NULL
  if (flip) geneMapVP <- flipVP
  geneMap <- LDheatmapMapNew.add (nsnps, genetic.distances=genetic.distances,
                                  geneMapLocation=geneMapLocation,add.map,  
                                  geneMapLabelX=geneMapLabelX,
                                  geneMapLabelY=geneMapLabelY,
                                  distances=distances, vp=geneMapVP, 
                                  SNP.name=SNP.name, ind=ind, flip=flip)
  
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
} # function LDheatmap ends



preDrawDetails.ldheatmap <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}


postDrawDetails.ldheatmap <- function(x) {
  popViewport()
}


preDrawDetails.symbols <- function(x) {
  fontsize <- convertX(unit(1/20,"grobwidth", rectGrob()), "points")
  pushViewport(viewport(gp=gpar(fontsize=fontsize)))
}

postDrawDetails.symbols <- function(x) {
  popViewport()
}


