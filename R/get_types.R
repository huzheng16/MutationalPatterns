#' Retrieve base substitution types from vcf file
#' 
#' A function to extract the base substitutions and translate to the 6 base substitution types
#' @param vcf A CollapsedVCF object
#' @return types Base substitution types
#' @export

get_types = function(vcf) 
{
  muts = get_muts(vcf)
  types = unlist(muts)
  types = gsub('G>T', 'C>A', types)
  types = gsub('G>C', 'C>G', types)
  types = gsub('G>A', 'C>T', types)
  types = gsub('A>T', 'T>A', types)
  types = gsub('A>G', 'T>C', types)
  types = gsub('A>C', 'T>G', types)
  return(types)
}