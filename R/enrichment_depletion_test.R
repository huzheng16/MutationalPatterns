#' Enrichment depletion test
#' 
#' Test whether specific genomic region is depleted or enriched using binomial test
#' @param n_muts Number of mutations in whole surveyed genome
#' @param observed Observed number of mutations in surveyed region
#' @param surveyed_length Number of bp that was surveyed in whole genome
#' @param surveyed_region_length Number of bp that is surveyed in region
#' @export


enrichment_depletion_test = function(prob, observed, surveyed_region_length)
{
  expected = prob * surveyed_region_length
  
  if(observed < expected)
  {
    # For depletion
    # do lower tail test
    pval = pbinom(observed, surveyed_region_length, prob, lower.tail=TRUE)
  }else{
    # For enrichment
    # do upper tail test
    pval = pbinom(observed-1, surveyed_region_length, prob, lower.tail=FALSE)
  }
  if(pval < 0.05){
    significant="*"
  }else{significant=""} 
  res = data.frame(pval, significant)
  return(res)
}







