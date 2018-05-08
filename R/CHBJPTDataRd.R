#' @name CHBJPTData
#' @aliases CHBJPTData
#' CHBJPTSNP
#' CHBJPTDist
#' hapmapJPT
#' hapmapCHB
#' @docType data
#' @title Example of data set for LDHeatmap
#' @description CHBJPTSNP:  Genotypes on 13 SNPs for 45 Chinese and 45 Japanese people
#' CHBJPTDist: Physical map positions of the 13 SNPs
#' @usage data(CHBJPTData)
#' @format CHBJPTSNP: A dataframe of SNP genotypes.
#'  Each row represents an individual.
#'  Each column represents a SNP.
#'
#'  CHBJPTDist: a vector of integers, representing SNP physical map
#'  locations on the chromosome.
#' @details The data frame \code{CHBJPTSNP} contains genotypes for 13 SNPs on chromosome 7,
#'  from 45 Chinese and 45 Japanese individuals.
#'  The Chinese individuals were unrelated residents of the community at Beijing Normal
#'  University with at least 3 Han Chinese grandparents.
#'  The Japanese individuals were unrelated residents of the Tokyo metropolitan
#'  area with all grandparents from Japan. The data are from release 21 of the
#'  International HapMap project (The International HapMap Consortium 2005).
#' @source International HapMap Project \url{www.hapmap.org}
#' @references The International HapMap Consortium. A haplotype map of
#' the human genome. Nature 437, 1299-1320. 2005.
#' @examples  data(CHBJPTData)
#'#Now do our panel plot with LDheatmaps in the panels
#'library(lattice)
#'pop<-factor(c(rep("chinese", 45), rep("japanese", 45)))
#'xyplot(1:nrow(CHBJPTSNP) ~ 1:nrow(CHBJPTSNP) | pop, type="n",
#'       scales=list(draw=FALSE), xlab="", ylab="",
#'       panel=function(x, y, subscripts,...) {
#'         LDheatmap(CHBJPTSNP[subscripts,], CHBJPTDist, newpage=FALSE)})
#' @keywords datasets
NULL
