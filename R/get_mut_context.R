#' Retrieve context of base substitutions
#' 
#' A function to extract the bases 3' upstream and 5' downstream of each position in vcf 
#' @param vcf A CollapsedVCF object
#' @param ref_genome 
#' @return trinucleotides
#' @export

# ref_genome eigenlijk meegeven als variable in functie? En dan ook in functies die deze gebruiken ????????????

get_mut_context = function(vcf) 
{
  vcf_context = as.character(getSeq(get(ref_genome), seqnames(vcf), start(vcf) - 1, end(vcf) + 1))
  return(vcf_context)
}