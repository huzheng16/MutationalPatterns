#' Extract mutational signatures using NMF
#' 
#' Decomposes trinucleotide count matrix into signatures and contribution of those signatures to the spectra of the individuals samples/vcf files
#' @param vcf_files Character vector of vcf_files
#' @param rank Number of signatures one would like to extract
#' @param nrun Number of iterations, default = 200
#' @param outdir Specify output directory
#' @export

extract_signatures = function(mut_spectrum, rank, nrun = 200)
{
  mut_spectrum = t(as.matrix(mut_spectrum))
  # add small pseudocount
  mut_spectrum = mut_spectrum + 0.0001
  # Calculate nmf
  print("Decomposing matrix using NMF...")
  res = NMF::nmf(mut_spectrum, rank = rank, method = "brunet", nrun=nrun, seed = 123456)
  print(paste("Number of iterations:", nrun))
  # Find signatures and contribution of signatures
  signatures = NMF::basis(res)
  contribution = NMF::coef(res)
  return(list(mut_spectrum = mut_spectrum, signatures = signatures, contribution = contribution))
}

