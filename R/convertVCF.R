
convertVCF <- function(vcf, phased = TRUE, samples = NULL, verbose = FALSE) {
  
  snp <- vcfR::read.vcfR(vcf, verbose = verbose)
  GT <- snp@gt[,!colnames(snp@gt)%in%"FORMAT"]
  
  if (!is.null(samples)) {
    GT <- GT[,colnames(GT)%in%samples]
    if (ncol(GT) == 0) stop ("could not find entries with the provided sample names") 
  }
  
  GT <- t(GT)
  if ("ID"%in%colnames(snp@fix)) colnames(GT) <- snp@fix[,"ID"]
  
  mat <- GT_to_numeric(GT, phased)
  
  if (phased == FALSE) {
    requireNamespace("snpStats")
    mat <- as(mat, "SnpMatrix")
  }
  
  return(mat)
  
}