#' @name convertVCF
#' @aliases convertVCF
#' @title Convert a VCF file into genotype data
#' @description Convert a VCF file into a list composed of genetic distances, subject IDs, and a SnpMatrix/XSnpMatrix object.
#' 
#' @usage convertVCF(file, phased = NULL, subjects = NULL, ...)
#' 
#' @param file A filename for a variant call format (vcf) file.
#' 
#' @param phased If \code{TRUE} the output genotype data are in the form of a XSnpMatrix object. Otherwise, they are in the form of a SnpMatrix object. 
#'If it is unspecified, the phasing status will be determined by checking the first entry in the GT section of the VCF file.
#'
#' @param subjects A character or factor containing subject IDs. If supplied, genotype info of only those subjects will be returned. 
#'This should be a subset of the sample IDs in the VCF file.
#'
#' @param ... Additional arguments to be passed to vcfR::read.vcfR().
#' 
#' @return A list which contains the following components:
#' \item{genetic.distances}{ A numeric vector of the reference positions of SNPs. }
#' \item{subjectID}{ A character vector of IDs of the subjects which the returned genotype data belong to. }
#' \item{data}{An object of SnpMatrix/XSnpMatrix class from the 'snpStats' package containing genotype data. }
#' 
#'@seealso  \code{\link[vcfR]{read.vcfR}}, \code{\link[snpStats]{SnpMatrix-class}}, \code{\link[snpStats]{XSnpMatrix-class}}
#'
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

###########################################################################

convertVCF <- function(file, phased = NULL, subjects = NULL, ...) {
  
  # read vcf file
  snp <- vcfR::read.vcfR(file, ...)
  
  # check validity of phased
  if (!is.null(phased) & !is.logical(phased)) stop("Invalid input for parameter phased, must be a logical constant or NULL.")
  
  # extract genetic distances
  if ("POS"%in%colnames(snp@fix)) genetic.distance = as.numeric(snp@fix[,"POS"])
  
  # extract genotype info
  GT <- snp@gt[,!colnames(snp@gt)%in%"FORMAT"]
  
  # extract subject IDs
  subjectID = colnames(GT)
  
  # subset GT if sample names are provided
  if (!is.null(subjects)) {
    GT <- GT[,colnames(GT)%in%subjects]
    if (ncol(GT) == 0) stop ("could not find entries with the provided sample names") 
  }
  
  # extract phasing info if it is not provided
  if (is.null(phased)) {
    sep <- unlist(strsplit(GT[1,1], ""))[[2]]
    if (sep == "|") {
      phased <- TRUE
    } else if (sep == "/") {
      phased <- FALSE
    }
  }
  
  # extract and add snp identifiers if any
  if ("ID"%in%colnames(snp@fix)) rownames(GT) <- snp@fix[,"ID"]
  
  # convert GT to SnpMatrix if phased if FALSE, else convert to XSnpMatrix
  mat <- GT_to_SnpMatrix(GT, phased)
  
  
  return(list(genetic.distances = genetic.distance, subjectID = subjectID, data = mat))
  
}