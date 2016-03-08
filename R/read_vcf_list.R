#' Read vcf files into list of CollapsedVCF objects
#' 
#' Function reads Variant Call Format (VCF) files into a CollapsedVCF object and combines them in a list object
#' @param vcf_file_list Character vector of vcf file names
#' @param sample_names Character vector of sample names
#' @return List of CollapsedVCF objects
#' @export
#' @examples 
#' vcf_file_list = list.files("/your_vcf_dir/", full.names = T) 
#' vcf_list = read_vcf_list(vcf_files_list, sample_names)

read_vcf_list = function(vcf_file_list, sample_names)
{
  vcf_list = list()
  for(i in 1:length(vcf_file_list))
  {
    new_vcf = read_vcf(vcf_file_list[i])
    new_vcf = list(new_vcf)
    names(new_vcf) = sample_names[i]
    vcf_list = c(vcf_list, new_vcf)
  }
  return(vcf_list)
}

