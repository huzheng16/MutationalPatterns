#' Read a bed file as a Granges object
#' 
#' Read a bed file as a Granges object
#' @param bed_file Location of your BED-file
#' @return GR A GRanges object
#' @import GenomicRanges
#' @export
#' @examples
#' bed_to_granges("my_file.bed")

bed_to_granges = function(bed_file)
{
  bed = read.table(bed_file, header = F, stringsAsFactors = F)
  chr = paste("chr", bed[,1], sep="")
  # Convert BED (0-based) start postion to Granges (1-based)
  start = bed[,2] + 1
  # In BED end position is excluded, in Granges end position is included -> +1 -1 -> no conversion needed
  end = bed[,3]
  GR = GRanges(chr, IRanges(start,end))  
  return(GR)
}


