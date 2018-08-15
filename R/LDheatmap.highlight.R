#' @name LDheatmap.highlight
#' @aliases LDheatmap.highlight
#' @title Highlight a genetic region in the linkage disequilibrium heat map
#' @description The function \code{LDheatmap.highlight()} is used to highlight a
#'specified genetic region in the linkage disequilibrium (LD)
#'heat map drawn with the \code{\link{LDheatmap}()} function.
#' @usage LDheatmap.highlight(LDheatmap, i, j, fill = "NA", col = "black", lwd = 1, 
#'lty = 1,flipOutline=FALSE, crissCross = FALSE)
#' @param LDheatmap An object of class \code{"LDheatmap"} returned
#'by the function \code{LDheatmap()}.
#' @param i A numeric value specifying the index of the first
#'SNP to be in the highlighted region.
#' @param j A numeric value specifying the index of the last SNP,
#'which must be different from \code{i}, to be in the highlighted region.
#' @param fill Color to fill the highlighted area with.
#' @param col A character string specifying the color of the line
#'segments defining the boundary of highlighted region; see
#'\code{\link[graphics]{par}()} for possible values.
#' @param lwd A \emph{positive} number specifying the width of the
#'boundary segments.
#' @param lty Either an integer or a character string specifying the
#'line type of the boundary segments; see \code{\link[graphics]{par}()}
#'for possible values.
#' @param flipOutline A Boolean variable that flips the outlined section 
#'over the diagonal of the heatmap.
#' @param crissCross A Boolean variable that controls whether a contiguous 
#'selection of SNPs are outlined only on their polygonal boundary or at 
#'individual SNP levels.
#' @return A data frame of the x and y coordinates of points defining the
#'border of the highlighted area.
#' @note The function \code{LDheatmap.highlight()} highlights the cells representing
#'the pairwise LD for the SNPs located between \code{i}-th and \code{j}-th (inclusive)
#'SNPs in the genomic region of interest.
#'The order of indices has no effect on the plot.  For example,
#'\code{LDheatmap.highlight(LDheatmap, i=2, j=4)} is the same as
#'\code{LDheatmap.highlight(LDheatmap, i=4, j=2)}, which highlights
#'the cells representing the pairwise LD for the second,
#'third and fourth SNPs.
#' @section Warning: By default, \code{LDheatmap.highlight()} finds the viewport to draw on from
#'the \code{LDheatmap} object passed to it as an argument.
#'However, if \code{LDheatmap()} was called with the option \code{pop=TRUE},
#'the resulting \code{LDheatmap} object is not assigned a
#'viewport. In this case, \code{LDheatmap.highlight()} assumes
#'the user wishes to highlight in the current viewport.
#'Therefore, if \code{LDheatmap()}
#'has been called with the option \code{pop=TRUE},
#'the user must navigate to the correct viewport
#'before calling \code{LDheatmap.highlight()}.
#' @author Nicholas Lewin-Koh
#'<nikko@hailmail.net>,
#'Ji-Hyung Shin <shin@sfu.ca>, Sigal Blay <sblay@sfu.ca>
#' @examples data(CEUData)
#'tt <- LDheatmap(CEUSNP, genetic.distances=CEUDist)
#'LDheatmap.highlight(tt, 3, 8, col="blue", fill="green", lwd=3, flipOutline=FALSE, crissCross=FALSE)
#' @keywords aplot
#' @export


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


#_______________________Highlight a region in the heatmap____________________________##
LDheatmap.highlight <- function(LDheatmap, i, j, fill="NA", col="black", lwd=1, 
                                lty=1, flipOutline=FALSE, crissCross = FALSE){

  
  requireNamespace("grid")
  # Highlights the perimeter of selected cells in the heatmap as a block
  backbone <- function(i,j,nSNP){
    x <- c(i-1,i-1,j-1)/nSNP
    y <- c(i,j,j)/nSNP
    cbind(x,y)
  }
  backboneFlip <- function(i,j,nSNP){
    x <- c(i,j,j)/nSNP
    y <- c(i-1,i-1,j-1)/nSNP
    cbind(x,y)
  }
  
  # # # #
  rectVert <- function(i, j , nSNP){
    rectangles <- data.frame()
    for( k in i:(j-1)){
      x <- c(k, k, k + 1, k + 1) / nSNP
      y <- c(k, j, k, j) / nSNP
      coords <- cbind(x, y)
      rectangles <- rbind(rectangles, coords)
    }
    return(rectangles)
  }
  
  rectHorizontal <- function(i, j ,nSNP){
    rectangles <- data.frame()
    for(m in i:(j-1)){
      x <- c(i-1, m , i-1, m) / nSNP
      y <- c(m+1, m+1, m+2, m+2) / nSNP
      coords <- cbind(x, y)
      rectangles <- rbind(rectangles, coords)
    }
    return(rectangles)
  }
  # # # #
  
  zigzag <- function(i,j,nSNP){
    c1 <- j-i
    nvert <- (2*c1)-1
    x <-c(j-1,rep((j-2):(j-c1),each=2))
    y <- c(rep((j-1):(j-(c1-1)),each=2),j-c1)
    cbind(x,y)/nSNP 
  }
  zigzagFlip <- function(i,j,nSNP){
    c1 <- j-i
    nvert <- (2*c1)-1
    y <-c(j-1,rep((j-2):(j-c1),each=2))
    x <- c(rep((j-1):(j-(c1-1)),each=2),j-c1)
    cbind(x,y)/nSNP 
  }
  
  nSNP <- dim(LDheatmap$LDmatrix)[1]
  if(length(i)>1 | length (j) > 1) stop("i and j must be scalar indices")
  if((i<1 | i>nSNP) |(j<1 | j>nSNP) )
    stop(paste("index out of bounds, i and j must be in (1,",nSNP,")",sep=""))
  if(i==j) stop("i cannot be equal to j")
  if(i>j){
    h<-i
    i <- j
    j <- h
  }
  pgon <- data.frame(rbind(backbone(i,j,nSNP), zigzag(i,j,nSNP)))
  if(!is.null(LDheatmap$flipVP)) pgon <- data.frame(rbind(backboneFlip(i,j,nSNP), zigzagFlip(i,j,nSNP)))
  ## Square or almost square interior Blocks
  names(pgon) <- c("x","y")
  
  
  # # # # For the grid highlight case
  vertRectangles <- rectVert(i, j, nSNP = dim(LDheatmap$LDmatrix)[1]) 
  horizonRectangles <- rectHorizontal(i, j, nSNP = dim(LDheatmap$LDmatrix)[1])
  names(vertRectangles) <- c("x", "y")
  names(horizonRectangles) <- c("x", "y")
  # # # #
  
  
  heatmap.vp <- LDheatmap$heatmapVP$name
  #If heatmap.vp is on the grid display list, i.e., it is included in the 
  #returned value of current.vpTree(), a[1]=1 else a[1]=NA:
  a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), fixed=TRUE)
  if(!is.na(a[1]))   seekViewport(heatmap.vp)
  else               pushViewport(LDheatmap$heatmapVP)
  if (!is.null(LDheatmap$flipVP)) pushViewport(LDheatmap$flipVP)
  
  # Added section #
  if(flipOutline == T){
    tempy <- pgon$y
    tempx <- pgon$x
    pgon$y <- tempx
    pgon$x <- tempy
  }
  highlight <- polygonGrob(x=pgon$x, y=pgon$y, 
                           gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="highlight")
  
  # # # #
  if(crissCross == TRUE){
    
    for(i in 1:(dim(vertRectangles)[1]/4)){
      width <- vertRectangles$x[(i-1)*4 + 3] - vertRectangles$x[(i-1)*4 + 1]
      height <- vertRectangles$y[(i-1)*4 + 1] - vertRectangles$y[(i-1)*4 + 2]
      if(is.null(LDheatmap$flipVP)){
        oneRect <- rectGrob(x = vertRectangles$x[(i-1)*4+1] - width/2, y = vertRectangles$y[(i-1)*4+1] - height/2, 
                            width = width,
                            height= height,
                            gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="rect")
      }
      else{
        # Flip is swapping of x-y coordinates, therefore reverse assignment of width and height
        width <- vertRectangles$y[(i-1)*4 + 1] - vertRectangles$y[(i-1)*4 + 2]
        height <- vertRectangles$x[(i-1)*4 + 3] - vertRectangles$x[(i-1)*4 + 1]
        oneRect <- rectGrob(x = vertRectangles$y[(i-1)*4+1] - width/2, y = vertRectangles$x[(i-1)*4+1] - height/2,
                            width = width,
                            height= height,
                            gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="rect")
      }
      
      grid.draw(oneRect)
    }
    for(j in 1:(dim(horizonRectangles)[1]/4)){
      width <- horizonRectangles$x[(j-1)*4 + 2] - horizonRectangles$x[(j-1)*4 + 1]
      height<- horizonRectangles$y[(j-1)*4 + 3] - horizonRectangles$y[(j-1)*4 + 1]
      if(is.null(LDheatmap$flipVP)){
        oneRect <- rectGrob(x = horizonRectangles$x[(j-1)*4 + 2] - width/2, y = horizonRectangles$y[(j-1)*4 + 1] - height/2,
                            width = width,
                            height = height,
                            gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="rect")
      }
      else{
        # Flip is swapping of x-y coordinates, therefore reverse assignment of width and height
        width <- horizonRectangles$y[(j-1)*4 + 3] - horizonRectangles$y[(j-1)*4 + 1]
        height<- horizonRectangles$x[(j-1)*4 + 2] - horizonRectangles$x[(j-1)*4 + 1]
        oneRect <- rectGrob(x = horizonRectangles$y[(j-1)*4 + 1] - width/2, y = horizonRectangles$x[(j-1)*4 + 2] - height/2,
                            width = width,
                            height = height,
                            gp=gpar(col=col, fill=fill, lwd=lwd, lty=lty), name="rect")
      }
      grid.draw(oneRect)
    }
  }
  # # # #
  
  
  grid.draw(highlight)
  if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
  else              popViewport() 
  invisible(pgon)
}


