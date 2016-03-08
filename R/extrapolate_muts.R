#' Extract sample name from path to file
#' 
#' Function to extract the sample name from path to file
#' @param vcf CollapsedVCF object
#' @param GR_surveyed GRanges object of surveyed bed file
#' @return df Dataframe with number of mutations extrapolated to whole genome and other stats
#' @export


extrapolate_muts = function(vcf, GR_surveyed, genome_length)
{
  n_muts = dim(vcf)[1]
  surveyed_length = sum(as.numeric(width(GR_surveyed)))
  percentage_surveyed = round( (surveyed_length / genome_length) * 100 ,1)
  # extrapolate n mutations to the whole nonN autosomal genome
  extrapolated_n_muts = round( (n_muts * genome_length) / surveyed_length ,1)
  df = data.frame(n_muts, surveyed_length, percentage_surveyed, extrapolated_n_muts)
  return(df)
}





