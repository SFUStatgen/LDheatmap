
convertVCF <- function(vcf, phased = TRUE, samples = NULL, verbose = TRUE) {
  
  # read vcf file
  snp <- vcfR::read.vcfR(vcf, verbose = verbose)
  
  # extract genetic distances
  if ("POS"%in%colnames(snp@fix)) genetic.distance = snp@fix[,"POS"]
  
  # extract genotype info
  GT <- snp@gt[,!colnames(snp@gt)%in%"FORMAT"]
  
  # extract subject IDs
  subjectID = colnames(GT)
  
  # subset GT if sample names are provided
  if (!is.null(samples)) {
    GT <- GT[,colnames(GT)%in%samples]
    if (ncol(GT) == 0) stop ("could not find entries with the provided sample names") 
  }
  
  # extract and add snp identifiers if any
  if ("ID"%in%colnames(snp@fix)) rownames(GT) <- snp@fix[,"ID"]
  
  # convert GT to numeric
  # mat now has subjects in rows and mutations in columns
  mat <- GT_to_numeric(GT, phased)
  
  # convert GT to snpMatrix if genotypes are unphased
  # else convert to raw type
  mat <- numeric_to_raw_prep(mat, phased)
  mode(mat) <- "raw"
  
  if (phased != TRUE) mat <- new("SnpMatrix", mat)
  
  return(list(genetic.distance = genetic.distance, subjectID = subjectID, data = mat))
  
}