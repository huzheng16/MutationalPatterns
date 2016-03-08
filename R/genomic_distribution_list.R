#' Find overlaps between mutations and a genomic region
#' 
#' Function finds the number of mutations that reside in genomic region and takes surveyed area of genome into account
#' @param MutPat_object A list...
#' @param region_list List with GRanges objects containing locations of genomic regions
#' @export


genomic_distribution_list = function(MutPat_object, region_list)
{
  df = data.frame()
  # for each region j
  for(j in 1:length(region_list))
  {
    # for each sample i
    for(i in 1:length(MutPat_object$vcf))
    {
      sample = names(MutPat_object$vcf)[i]
      res = genomic_distribution(MutPat_object$vcf[[i]], MutPat_object$surveyed[[i]], region_list[[j]])
      res$region = names(region_list)[j]
      res$sample = sample
      res$type = MutPat_object$type[i]
      res$individual = MutPat_object$individual[i]
      df = rbind(df, res)
    }
  }
  return(df)
}






