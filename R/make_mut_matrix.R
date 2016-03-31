#' Make 96 trinucleotide mutation frequency matrix
#'  
#' @description
#' @param vcf_list List of collapsed vcf objects from which one would like to contstruct a count matrix
#' @return 96 mutation count matrix
#' @import GenomicRanges
#' @export

make_mut_matrix = function(vcf_list, ref_genome)
{
  df = data.frame()
  for(vcf in vcf_list)
  {
    type_context = get_type_context(vcf, ref_genome)
    row = mut_96_occurences(type_context)
    df = rbind(df, row)
  }
  names(df) = names(row)
  row.names(df) = names(vcf_list)
  return(df)
}