#' Retrieve context of base substitutions
#' 
#' A function to extract the bases 3' upstream and 5' downstream of the base substitutions
#' @param vcf A CollapsedVCF object
#' @param ref_genome Reference genome
#' @return Character vector the context of the base substitutions
#' @export

get_mut_context = function(vcf, ref_genome) 
{
  vcf_context = as.character(getSeq(get(ref_genome), seqnames(vcf), start(vcf) - 1, end(vcf) + 1))
  return(vcf_context)
}