#' Read vcf files into list of CollapsedVCF objects
#' 
#' Function reads Variant Call Format (VCF) files into a CollapsedVCF object and combines them in a list object
#' @param vcf_files Character vector of vcf file names
#' @param sample_names Character vector of sample names
#' @return List of CollapsedVCF objects
#' @export
#' @examples 
#' vcf_file_list = list.files("/your_vcf_dir/", full.names = T) 
#' vcf_list = read_vcf_list(vcf_files_list, sample_names)

read_vcf = function(vcf_files, sample_names, genome = "-")
{
  if(!(length(vcf_files) == length(sample_names))){stop("Provide the same number of sample names as vcf files")}
  vcf_list = list()
  for(i in 1:length(vcf_files))
  {
    vcf = readVcf(vcf_files[i], genome)
    # remove positions with multiple alternative alleles
    mult_alt_all = which(lapply(alt(vcf), function(x) length(x) > 1) == TRUE)
    if(length(mult_alt_all > 0)){
      warning(paste(sample_names[i],"contains", length(mult_alt_all), "position(s) with multiple alternative alleles. These positions were excluded."))
      vcf = vcf[-mult_alt_all]
    }
    # remove indel positions
    insertions = which(unlist(lapply(alt(vcf), function(x) width(x))) > 1)
    deletions = which(width(ref(vcf)) > 1)
    indels = c(insertions, deletions)
    if(length(indels > 0)){
      warning(paste(sample_names[i],"contains", length(indels), "indel position(s). These positions were excluded."))
      vcf = vcf[-indels]
    }
    
    # add "chr" to chromosomes if not there already
    # if(length(grep("chr", seqlevels(vcf))) == 0){seqlevels(vcf) = paste('chr', seqlevels(vcf), sep = '')}
    vcf = list(vcf)
    names(vcf) = sample_names[i]
    vcf_list = c(vcf_list, vcf)
  }
  return(vcf_list)
}

