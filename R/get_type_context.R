#' Retrieve context of base substitution types
#' 
#' A function to extract the bases 3' upstream and 5' downstream of a base substitution type
#' @param vcf A CollapsedVCF object
#' @return trinucleotides
#' @export

get_type_context = function(vcf)
{
  mut_context = get_mut_context(vcf)
  muts = get_muts(vcf)
  types = get_types(vcf)  
  # find the mutations for which the context needs to be adjusted
  x = which(muts != types)
  # subset mut_context
  y = mut_context[x]
  # change the context of these mutations to reverse complement of the context
  y = chartr('ATGC', 'TACG', y)
  y = reverse(y)
  # replace subset with reverse complement
  mut_context[x] = y
  return(list(types, mut_context))
}