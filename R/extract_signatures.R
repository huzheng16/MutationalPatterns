#' Extract mutational signatures from 96 mutation matrix using NMF
#' 
#' Decomposes trinucleotide count matrix into signatures and contribution of those signatures to the spectra of the samples/vcf files
#' @param mut_matrix 96 mutation matrix 
#' @param rank Number of signatures to extract
#' @param nrun Number of iterations, default = 200
#' @return Named list of mutation matrix, signatures and signature contribution
#' @export

extract_signatures = function(mut_matrix, rank, nrun = 200)
{
  mut_matrix = t(as.matrix(mut_matrix))
  # add small pseudocount
  mut_matrix = mut_matrix + 0.0001
  # Calculate nmf
  print("Decomposing matrix using NMF...")
  res = NMF::nmf(mut_matrix, rank = rank, method = "brunet", nrun=nrun, seed = 123456)
  print(paste("Number of iterations:", nrun))
  # Find signatures and contribution of signatures
  signatures = NMF::basis(res)
  contribution = NMF::coef(res)
  return(list(mut_matrix = mut_matrix, signatures = signatures, contribution = contribution))
}

