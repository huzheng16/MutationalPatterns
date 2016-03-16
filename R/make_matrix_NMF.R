#'Construct trinucleotide changes frequency matrix
#'  
#' @param vcf_list List of collapsed vcf objects from which one would like to contstruct a count matrix
#' @return df Count data.frame 
#' @import GenomicRanges
#' @export

make_matrix_NMF = function(vcf_list, ref_genome, type = "normal")
{
  df = data.frame()
  if(type == "normal")
  {
    for(vcf in vcf_list)
    {
      type_context = get_type_context(vcf, ref_genome)
      row = make_sample_row_96(type_context)
      df = rbind(df, row)
    }
  }
  if(type == "strand bias")
  {
    for(vcf in vcf_list)
    {
      row = make_sample_row_192(vcf)
      df = rbind(df, row)
    }
  }
  else if(type == "replication bias")
  {
    for(vcf in vcf_list)
    {
      row = make_sample_row_288(vcf)
      df = rbind(df, row)
    }
  }
  names(df) = names(row)
  row.names(df) = names(vcf_list)
  return(df)
}

