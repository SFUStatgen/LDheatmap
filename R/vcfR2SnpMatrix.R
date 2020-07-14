#' @name  vcfR2SnpMatrix
#' @aliases  vcfR2SnpMatrix
#' @title Extract genotype information from a vcfR object
#' @description Convert a \code{\link[vcfR]{vcfR-class}} object into a list composed of genetic distances, subject IDs, and a \code{\link[snpStats]{SnpMatrix-class}}/\code{\link[snpStats]{XSnpMatrix-class}} object.
#' 
#' @usage  vcfR2SnpMatrix(obj, phased = NULL, subjects = NULL)
#' 
#' @param obj An instance or path of the \code{\link[vcfR]{vcfR-class}} object to be processed.
#' 
#' @param phased If \code{TRUE} the output genotype data are in the form of a \code{\link[snpStats]{XSnpMatrix-class}} object. Otherwise, they are in the form of a \code{\link[snpStats]{SnpMatrix-class}} object. 
#'If it is unspecified, the phasing status will be determined by checking the first entry in the gt slot of the vcfR object. If the first entry is also missing, the value will be set to \code{FALSE}.
#'
#' @param subjects A character or factor containing subject IDs. If supplied, genotype info of only those subjects will be returned. 
#'This should be a subset of the sample IDs in the vcfR object.
#' 
#' @details
#' In order to let \code{ vcfR2SnpMatrix} function properly, the input \code{\link[vcfR]{vcfR-class}} object is expected to be generated from a valid VCF file which contains only biallelic SNPs and includes a GT section.
#' 
#' @return A list which contains the following components:
#' \item{genetic.distances}{ A numeric vector of the reference positions of SNPs. }
#' \item{subjectID}{ A character vector of IDs of the subjects which the returned genotype data belong to. }
#' \item{data}{An object of \code{\link[snpStats]{SnpMatrix-class}}/\code{\link[snpStats]{XSnpMatrix-class}} containing genotype data. }
#' 
#'@seealso  \code{\link[vcfR]{io_vcfR}}, \code{\link[snpStats]{SnpMatrix-class}}, \code{\link[snpStats]{XSnpMatrix-class}}
#'
#' @examples # Load the vcfR object -- requires vcfR
#'  if (requireNamespace("vcfR", quietly = TRUE)) {
#'    require(vcfR)
#'    data(vcfR_example)
#'    vcf <- vcf[8:12,]
#'    # Extract needed genotype information
#'    alist <- vcfR2SnpMatrix(vcf)
#' 
#'    # Draw a pairwise LD plot using the extracted data
#'    LDheatmap(alist$data, alist$genetic.distance)
#'  }
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

 vcfR2SnpMatrix <- function(obj, phased = NULL, subjects = NULL) {
   
  requireNamespace("snpStats")
   
  # check validity of phased
  if (!is.null(phased) & !is.logical(phased)) stop("Invalid input for parameter phased, must be a logical constant or NULL.")
  
  # extract genetic distances
  if ("POS"%in%colnames(obj@fix)) genetic.distance = as.numeric(obj@fix[,"POS"])
  
  # extract genotype info
  GT <- obj@gt[,!colnames(obj@gt)%in%"FORMAT"]
  
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
    } else {
      phased <- FALSE
    }
  }
  
  # extract and add snp identifiers
  if ("ID"%in%colnames(obj@fix)) rownames(GT) <- obj@fix[,"ID"]
  
  # convert GT to SnpMatrix if phased is FALSE, else convert to XSnpMatrix
  mat <- GT_to_SnpMatrix(GT, phased)
  
  
  return(list(genetic.distances = genetic.distance, subjectID = subjectID, data = mat))
  
}