#' @name LDheatmap.marks
#' @aliases LDheatmap.marks
#' @title Plots a symbol in the centers of cells of the heat map image
#' @description The function \code{LDheatmap.marks()} is used to plot
#'a symbol in the centers of cells representing the
#'pairwise linkage disequilibria of specified pairs of SNPs.
#' @usage LDheatmap.marks(LDheatmap, i, j = NULL, pch = 20, gp=gpar(...), ...)
#' @param LDheatmap An object of class \code{"LDheatmap"}
#'returned by the function \code{\link{LDheatmap}()}.
#' @param i A vector of indices of the first set of SNPs.
#' @param j A vector of indices of the second set of SNPs.
#' @param pch Either an integer value or a single character specifying
#'the symbol to be plotted. See \code{\link[graphics]{points}()}
#'for possible values and their corresponding symbols.
#' @param gp Graphical parameters; See \code{\link[grid]{gpar}()}.
#' @param ... Graphical parameter settings to be passed on to the \code{\link[grid]{gpar}()}
#'function.
#' @return \item{ x }{ The vector of x coordinate(s) of the plotted symbol(s).}
#' \item{ y }{ The vector of y coordinate(s) of the plotted symbol(s). }
#' @details The lengths of the vectors \code{i} and \code{j} must be the same and
#'greater than or equal to 1.
#'If the lengths are greater than 1, the function plots the specified
#'symbol in the centers of the (i\eqn{\mbox{\textasciicircum}}{^}k,
#'                              j\eqn{\mbox{\textasciicircum}}{^}k)-th cells (for k=1,...K; K =
#'                                                                              length of the vectors \code{i} and \code{j}), where
#'i\eqn{\mbox{\textasciicircum}}{^}k and
#'j\eqn{\mbox{\textasciicircum}}{^}k are
#'the k-th elements of vectors \code{i} and \code{j}, respectively.
#'For example, if \code{i=c(1,2)} and \code{j=c(3,5)}, \code{LDheatmap()}
#'plots a symbol in the centers of the cells representing pairwise
#'linkage disequilibria between the first and third SNPs and between the
#'second and fifth SNPs in the genome of interest.  Note that the order
#'of the sets of indices does not matter; for example,
#'\code{LDheatmap.marks(LDheatmap, i=c(1,2), j=c(3,5))} is equivalent
#'to \code{LDheatmap.marks(LDheatmap, i=c(3,5), j=c(1,2))}.
#' @section Warning: By default, \code{LDheatmap.marks()} finds the viewport to draw on from
#'the \code{LDheatmap} object passed to it as an argument.
#'However, if \code{LDheatmap()} was called with the option \code{pop=TRUE},
#'the resulting \code{LDheatmap} object is not assigned a
#'viewport. In this case, \code{LDheatmap.marks()} assumes
#'the user wishes to highlight in the current viewport.
#'Therefore, if \code{LDheatmap()}
#'has been called with the option \code{pop=TRUE},
#'the user must navigate to the correct viewport
#'before calling \code{LDheatmap.marks()}.
#' @author Nicholas Lewin-Koh <nikko@hailmail.net>, Ji-Hyung Shin <shin@sfu.ca>,
#'Sigal Blay <sblay@sfu.ca>
#' @examples data(CEUData)
#'tt <- LDheatmap(CEUSNP, genetic.distances=CEUDist)
#'LDheatmap.marks(tt, 15, 3, cex=1.6, col="blue")
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



#_______________________Add marks to heatmap__________________________________##
# Adds a symbol to the i,jth cell of the heatmap. 
# Does this work when flip=TRUE?
# The default is to add a symbol to the diagonal (j=i). 
# i and j can be vectors.
LDheatmap.marks <- function(LDheatmap, i, j=NULL, pch=20, gp=gpar(...), ...){
    requireNamespace("grid")
    nSNP <- dim(LDheatmap$LDmatrix)[1]
    if(is.null(j)) j<-i
    ind <- i>j
    if(any(ind)){
      ind <- as.numeric(ind)
      ind[ind>0] <- i[ind>0]
      i[ind>0] <- j[ind>0]
      j[ind>0] <- ind[ind>0]
    }
    pts<-list(x=(i-0.5)*1/nSNP,y=(j-0.5)*1/nSNP)
    if(!is.null(LDheatmap$flipVP)) pts<-list(x=(j-0.5)*1/nSNP,y=(i-0.5)*1/nSNP)
    heatmap.vp <- LDheatmap$heatmapVP$name
    #If heatmap.vp is on the grid display list, i.e., it is included in the 
    #returned value of current.vpTree(), a[1] <- 1 else a[1] <- NA
    a <- grep(paste("[", heatmap.vp, "]", sep=""), as.character(current.vpTree()), fixed=TRUE)
    if(!is.na(a[1]))   seekViewport(heatmap.vp)
    else               pushViewport(LDheatmap$heatmapVP)
    if (!is.null(LDheatmap$flipVP)) pushViewport(LDheatmap$flipVP)
    Symbols <- pointsGrob(pts$x, pts$y, pch=pch, gp=gp, name="symbols")
    symbols <- gTree(children=gList(Symbols), name="Symbols", cl="symbols")
    grid.draw(symbols)
    if(!is.na(a[1]))  upViewport(0)  #back to the root viewport
    else              popViewport() 
    invisible(pts)
}


