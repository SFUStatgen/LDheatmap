#' @name LDheatmap.addScatterplot
#' @aliases LDheatmap.addScatterplot
#' @title Add a scatter plot to an LDheatmap object
#' @description Add a scatter plot to an LDheatmap object.
#'The x axis is the map of genetic distances of the SNPs.
#' @usage LDheatmap.addScatterplot(LDheatmap, P, height = 0.2, ylab = NULL, ylim=NULL,
#'type = "points")
#' @param LDheatmap An object of class LDheatmap.
#' @param P A vector with the values to be plotted as the y axis.
#' @param height The height of the plot.
#' @param ylab The y axis label.
#' @param ylim The y axis limits.
#' @param type Plot type. Possible values are \code{"points"} (the default), \code{"lines"} or \code{"both"}.
#' @details The function creates an \code{"association"} \code{grob} and adds it to the \code{LDheatmap} object.
#'Then it pushes a \code{viewport} and draws the \code{LDheatmapGrob} onto it.
#' @return An object of class LDheatmap given as an argument, with the \code{grob}
#'\code{LDheatmapGrob} modified to inclue the \code{"association"} child grob.
#' @note The function can display at most two scatter plots in the default setting. To add three or more scatter plots 
#'in the same viewport, the user can change the input value \code{"location"} from function \code{constructVP} which is 
#'the function inside the \code{LDheatmap.addScatterplot}. The default \code{"location"} value is 0.03, for adding the 
#'third scatter plot, user needs to set the \code{"location"} to 0.23, where 0.2 units is the height of the scatter plot.
#'For the fourth scatter plot, set the \code{"location"} to 0.43 etc. 
#'See Examples for usage.
#' @author Sigal Blay <sblay@sfu.ca> and more
#' @seealso \code{\link{LDheatmap}}
#' @examples # Load the package's data set
#'data("CEUData")
#'# Produce an LDheatmap object
#'MyLDheatmap <- LDheatmap(CEUSNP, genetic.distances = CEUDist, flip = TRUE)
#'# Generate an arbitrary vector of values to plot
#'Yvalues <- seq(length = length(MyLDheatmap$genetic.distances), from = 0.01, to = 0.5)
#'# Add scatter plot
#'assoc <- LDheatmap.addScatterplot(MyLDheatmap, Yvalues)
#'######## Adding three or more scatter plots ########
#'# Redefine LDheatmap.addScatterplot() to display the third scatter plot
#'LDheatmap.addScatterplot_test3 <- function(LDheatmap, P, height=0.2, ylab=NULL, 
#'ylim=NULL, type="points",color,pch) {
#'if (dim(LDheatmap$LDmatrix)[1] != length(P)) {
#'print("Length of vector not equal number of SNPs in LDheatmap")
#'return()
#'
#'flip <- !is.null(LDheatmap$flipVP)
#'vp <- constructVP(LDheatmap$LDheatmapGrob, 0.23, flip)
#'......
#'return(LDheatmap)
#'}}
#'environment(LDheatmap.addScatterplot_test3) <- asNamespace('LDheatmap')
#' 
#' @keywords hplot
#' @export


# LDheatmap - Plots measures of pairwise linkage disequilibria for SNPs
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


#______________________________________Add Scatter Plot________________________________##
LDheatmap.addScatterplot <- function(LDheatmap, P, height=0.2, ylab=NULL, ylim=NULL, type="points") {
  requireNamespace("grid")
  if (dim(LDheatmap$LDmatrix)[1] != length(P) ) {
     print ("Length of vector not equal number of SNPs in LDheatmap")
     return()
  }

  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, 0.03, flip)
  vp$height <- unit(height, "npc")
  vp$name <- "associationVP"
  if(is.null(ylim)) ylim <- c(floor(min(P)), ceiling(max(P)))
  vp$yscale <- ylim
  vp$xscale <- c(min(LDheatmap$genetic.distances), max(LDheatmap$genetic.distances))
  xaxis <- linesGrob(x=vp$xscale, y=0, default.units="native", name="xaxis")
  yaxis <- linesGrob(x=min(LDheatmap$genetic.distances), y=vp$yscale, default.units="native",
		name="yaxis") 
  yaxisT <- yaxisGrob(name="yaxis_ticks", gp=gpar(fontsize=7))
  ylab <- textGrob(ylab, rot=90, gp=gpar(fontsize=9), name="yaxis_title",
	x=unit(min(LDheatmap$genetic.distances), "native")- unit(10, "millimeters"))
  vpstack <- vp; if(flip) vpstack <- vpStack(LDheatmap$flipVP,vp)
  association <- gTree(children=gList(xaxis, yaxis, yaxisT, ylab), name="association", vp=vpstack)
  if (type=="points" || type == "both") {
     graph_points <- pointsGrob(LDheatmap$genetic.distances, P, size=unit(2, "millimeters"),
 		name="points")
     association <- addGrob(association, graph_points)
  }
  if (type=="lines" || type == "both") {
     graph_lines <- linesGrob(LDheatmap$genetic.distances, P, default.units="native", name="lines")
     association <- addGrob(association, graph_lines)
  }
  LDheatmap$LDheatmapGrob <- addGrob(LDheatmap$LDheatmapGrob, association)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, vp)
  return(LDheatmap)
}


