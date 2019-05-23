
convertVCF <- function(vcf, phased = TRUE, samples = NULL, verbose = FALSE) {
  
  # read vcf file and extract genotype info
  snp <- vcfR::read.vcfR(vcf, verbose = verbose)
  GT <- snp@gt[,!colnames(snp@gt)%in%"FORMAT"]
  
  # subset GT if sample names are provided
  if (!is.null(samples)) {
    GT <- GT[,colnames(GT)%in%samples]
    if (ncol(GT) == 0) stop ("could not find entries with the provided sample names") 
  }
  
  GT <- t(GT)
  # extract and add snp identifiers if any
  if ("ID"%in%colnames(snp@fix)) colnames(GT) <- snp@fix[,"ID"]
  
  # convert GT to numeric
  mat <- GT_to_numeric(GT, phased)
  
  # convert GT to snpMatrix if genotypes are unphased
  # else convert to raw type
  mat <- numeric_to_raw_prep(mat)
  mode(mat) <- "raw"
  
  if (phased != TRUE) mat <- new("SnpMatrix", mat)
  
  return(mat)
  
}