#' Retrieve base substitutions from vcf
#' 
#' A function to extract base substitutions of each position in vcf
#' @param vcf A CollapsedVCF object
#' @return Character vector with base substitutions
#' @import GenomicRanges
#' @export

get_muts = function(vcf) 
{
  ref = as.character(ref(vcf))
  alt = as.character(unlist(alt(vcf)))
  muts = paste(ref, alt, sep=">")
  return(muts)
}

