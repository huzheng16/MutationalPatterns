#' Read a bed file as a Granges object
#' 
#' Read a bed file as a Granges object
#' @param bed_file_list List of bed files
#' @return granges_list List of Granges objects
#' @import GenomicRanges
#' @export


bed_to_granges_list = function(bed_file_list)
{
  granges_list = list()
  for(i in 1:length(bed_file_list))
  {
    # sample name
    sample = extract_sample_name(bed_file_list[i])
    new_bed = bed_to_granges(bed_file_list[i])
    new_bed = list(new_bed)
    names(new_bed) = sample
    granges_list = c(granges_list, new_bed)
  }
  return(granges_list)
}