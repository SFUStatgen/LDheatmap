
convertVCF <- function(vcf, phased = NULL, samples = NULL, ...) {
  
  # read vcf file
  snp <- vcfR::read.vcfR(vcf, ...)
  
  # extract genetic distances
  if ("POS"%in%colnames(snp@fix)) genetic.distance = as.numeric(snp@fix[,"POS"])
  
  # extract genotype info
  GT <- snp@gt[,!colnames(snp@gt)%in%"FORMAT"]
  
  # extract subject IDs
  subjectID = colnames(GT)
  
  # subset GT if sample names are provided
  if (!is.null(samples)) {
    GT <- GT[,colnames(GT)%in%samples]
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
  
  
  return(list(genetic.distance = genetic.distance, subjectID = subjectID, data = mat))
  
}