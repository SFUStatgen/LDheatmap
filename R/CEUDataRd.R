#' @name CEUData
#' @aliases CEUData
#' CEUSNP
#' CEUDist
#' hapmapCEU
#' @docType data
#' @title Example data set for LDheatmap
#'@description  CEUSNP:  Genotypes on 15 SNPs for 60 people
#'
#'CEUDist: Physical map positions of the 15 SNPs in CEUSNP
#'
#' @usage data(CEUData)
#' @format CEUSNP: A dataframe of SNP genotypes.
#' Each row represents an individual.
#' Each column represents a SNP.
#'
#'CEUDist: A vector of integers, representing SNP physical map locations on the chromosome.
#' @details Data on SNPs with minor allele frequency greater
#'than 5\% from a 9kb region of chromosome 7 (base positions 126273659
#'                                            through 126282556 from release 7 of the International HapMap Project).
#'Genotypes from 30 parent-offspring trios (both
#'                                          parents, one offspring) were obtained.
#'The 30 trios are taken from the so-called CEPH families, a set of
#'multi-generational families from Utah with ancestry from northern and
#'western Europe. From this set of 90 people, 60 parents were extracted.
#' @source International HapMap Project \url{www.hapmap.org}
#' @references The International HapMap Consortium. A haplotype map of
#'the human genome. Nature 437, 1299-1320. 2005.
#'@examples data(CEUData)
#'@keywords datasets
NULL
