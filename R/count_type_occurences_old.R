#' Count the occurences of each substitution type
#' 
#' @param vcf_list A list of CollapsedVCF object
#' @return dataframe
#' @export

count_type_occurences = function(vcf_list)
{  
  n_samples = length(vcf_list)
  df = data.frame()
  print("Counting base substitution type occurrences for sample:")
  for(i in 1:n_samples)
  {
    print(i)
    vcf = vcf_list[[i]]
    types = get_types(vcf)
    df = rbind(df,as.vector(table(types)))
  }
  row.names(df) = names(vcf_list)
  colnames(df) = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  return(df)
}
