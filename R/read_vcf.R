#' Read vcf file
#' 
#' Function reads Variant Call Format (VCF) file into a CollapsedVCF object for hg19 genome
#' @param vcf_file Vcf file to be read
#' @return A CollapsedVCF object
#' @export

read_vcf = function(vcf_file)
{
  vcf = readVcf(vcf_file, "-")
  # check for chr in chromosomes
  # probably don't want this for the package!
  if(length(grep("chr", seqlevels(vcf))) == 0){seqlevels(vcf) = paste('chr', seqlevels(vcf), sep = '')}
  return(vcf)
}