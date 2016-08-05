#' Read bed files as a list of GRanges objects
#' 
#' Read bed files as a list of GRanges objects
#' @param bed_files Character vector with bed files
#' @param names Character vector with names of regions in bed files
#' @return granges_list List of Granges objects
#' @import GenomicRanges
#' @export


bed_to_granges = function(bed_files, names)
{
  if(!(length(bed_files) == length(names))) stop("Provide the same number of names as bed files")
  granges_list = list()
  for(i in 1:length(bed_files))
  {
    bed_file = bed_files[i]
    bed = read.table(bed_file, header = F, stringsAsFactors = F)
    chr = paste("chr", bed[,1], sep="")
    # Convert BED (0-based) start postion to Granges (1-based)
    start = bed[,2] + 1
    # In BED end position is excluded, in Granges end position is included -> +1 -1 -> no conversion needed
    end = bed[,3]
    new_bed = GRanges(chr, IRanges(start,end))  
    new_bed = list(new_bed)
    names(new_bed) = names[i]
    granges_list = c(granges_list, new_bed)
  }
  return(granges_list)
}