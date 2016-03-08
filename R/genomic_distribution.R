#' Find overlaps between mutations and a genomic region
#' 
#' Function finds the number of mutations that reside in genomic region and takes surveyed area of genome into account
#' @param vcf CollapsedVCF object containing mutations
#' @param surveyed GRanges object containing regions of the genome that were surveyed
#' @param region GRanges object containing locations of genomic region(s)
#' @export


genomic_distribution = function(vcf, surveyed, region)
{
  # number of mutations in vcf file
  n_muts = dim(vcf)[1]
  # number of bp that was surveyed
  surveyed_length = sum(as.numeric(width(surveyed)))
  
  # Intersect genomic region and surveyed region
  surveyed_region = intersect(surveyed, region)
  surveyed_region_length = sum(width(surveyed_region))
  
  # Find which mutations lie in surveyed genomic region
  overlap = findOverlaps(vcf, surveyed_region)
  muts_in_region = as.data.frame(as.matrix(overlap))$queryHits
  
  observed = length(muts_in_region)
  prob = n_muts / surveyed_length
  expected = prob * surveyed_region_length
  
  res = data.frame(n_muts, surveyed_length, surveyed_region_length, observed)
  return(res)
}