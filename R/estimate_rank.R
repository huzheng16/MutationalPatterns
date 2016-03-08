#' Estimate optimal rank for NMF decomposition
#' 
#' Find optimal number of signatures for NMF decomposition
#' @param vcf_files Character vector of vcf_files
#' @param rank_range Range of ranks one would like to test 
#' @param outdir Specify output directory
#' @export

estimate_rank = function(mut_spectrum, rank_range, nrun=100, outdir)
{
  # Add small pseudocount
  mut_spectrum = mut_spectrum + 0.0001
  # Check if rank_range is appropriate
  if(nrow(mut_spectrum) < max(rank_range))
  {
    stop("Maximum number of ranks should be smaller than the number of samples")
  }
  # Estimate ranks
  print("Estimating ranks...")
  estim.r = NMF::nmf(mut_spectrum, rank = rank_range, method = "brunet", nrun = nrun, seed = 123456)
  # plot result
  output = paste(outdir,"estim_rank.pdf", sep="")
  pdf(output)
  print(plot(estim.r))
  dev.off()
  print(paste("Output is written to ", output))
}