#' Extract sample name from path to file
#' 
#' Function to extract the sample name from path to file
#' @param file_name String of path to file
#' @return name String of sample name

extract_sample_name = function(filename)
{
  sample = tail(strsplit(filename, "/")[[1]],1)
  name =  unlist(strsplit(sample, ".vcf"))
  name = strsplit(name, "_")[[1]][1]
  return(name)
}