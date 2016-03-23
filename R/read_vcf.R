#' Read vcf file
#' 
#' Function reads Variant Call Format (VCF) file into a CollapsedVCF object for hg19 genome
#' @param vcf_file Vcf file to be read
#' @return A CollapsedVCF object
#' @export

read_vcf = function(vcf_file)
{
  vcf = readVcf(vcf_file, "-")
  # add "chr" to chromosomes if not there already
  if(length(grep("chr", seqlevels(vcf))) == 0){seqlevels(vcf) = paste('chr', seqlevels(vcf), sep = '')}
  return(vcf)
}