#' @name LDheatmap.addGenes
#' @aliases LDheatmap.addGenes
#' @title Add gene plot to an LDheatmap object.
#' @description Retrieve genes from the UCSC Genome Browser, generate the genes plot and add it to
#'an LDheatmap object.
#' @usage LDheatmap.addGenes(LDheatmap, chromosome, genome = NULL, genesLocation = 0.02,
#'splice_variants = TRUE, non_coding = TRUE)
#' @param LDheatmap An object of class LDheatmap.
#' @param chromosome A character string that identifies the chromosome.
#' @param genome The genome assembly to use. The default is the most recent human genome assembly on the UCSC genome browser.
#' @param genesLocation The gene plot distance from the LD heat map gene map.
#' @param splice_variants If \code{FALSE}, exclude gene splice variants.
#' @param non_coding If \code{FALSE}, exclude non-coding genes.
#' @details Note: The \code{LDheatmap} object should have a non-NULL \code{genetic.distances}
#'component. Otherwise the gene map will not be placed correctly.
#'The genes are color coded as follows:
#'  black -- feature has a corresponding entry in the Protein Data Bank (PDB);
#'dark blue -- transcript has been reviewed or validated by either the RefSeq, SwissProt or CCDS staff;
#'medium blue -- other RefSeq transcripts; and
#'light blue -- non-RefSeq transcripts.
#'
#'For assemblies older than hg18, all genes are plotted in grey.
#' @return An object of class LDheatmap given as an argument, with the \code{grob}
#'\code{LDheatmapGrob} modified to inclue the \code{"transcripts"} child grob.
#' @references \url{http://genome.ucsc.edu/cgi-bin/hgTrackUi?g=knownGene}
#' @author Sigal Blay <sblay@sfu.ca>
#' @seealso \code{\link{LDheatmap}}, \code{\link{plotGenes}}
#' @examples \dontrun{
#'data(GIMAP5.CEU)
#'ll<-LDheatmap(GIMAP5.CEU$snp.data,GIMAP5.CEU$snp.support$Position,flip=TRUE)
#'# Add gene plot
#'llplusgenes <- LDheatmap.addGenes(ll, chr="chr7", genome="hg18")
#'}
#' @keywords hplot
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

############################################################################



#_______________________Add genes from UCSC genome Browser to an LDheatmap_____________##
LDheatmap.addGenes <- function(LDheatmap, chromosome,  genome=NULL, genesLocation=0.02, splice_variants = TRUE, non_coding = TRUE) {
  if(is.null(LDheatmap$genetic.distances)) stop("LDheatmap must have genetic distances")
  requireNamespace("grid")
  minRange <- min(LDheatmap$genetic.distances)
  maxRange <- max(LDheatmap$genetic.distances)
#  minRange <- 150124000 #150434000
#  maxRange <- 150154477 #150220000

# gimap5
# minRange <- 149656848
# maxRange <- 150154477

# a range with a single gene
# minRange <- 149656848
# maxRange <- 149750000

  flip <- !is.null(LDheatmap$flipVP)
  vp <- constructVP(LDheatmap$LDheatmapGrob, genesLocation, flip)
  pushViewport(LDheatmap$LDheatmapGrob$vp) # the vps are necessary for the calculation of grobWidth
  if (flip) pushViewport(LDheatmap$flipVP) # in the plotGenes function
  Transcripts <- plotGenes(minRange, maxRange, chromosome, genome, plot_lines_distance=0.04, vp=vp, splice_variants)
  vp <- Transcripts$vp
  vpstack <- vp; if(flip) vpstack <- vpStack(LDheatmap$flipVP,vp)
  grobT <- editGrob(Transcripts, vp=vpstack)
  LDheatmap$LDheatmapGrob <-addGrob(LDheatmap$LDheatmapGrob, grobT)
  LDheatmap$LDheatmapGrob <- moveTitles(LDheatmap$LDheatmapGrob, vp)
  return(LDheatmap)
}


