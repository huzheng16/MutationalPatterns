# plot contribution with residual

plot_contribution = function(nmf_res)
{
  mut_spectrum  = nmf_res$mut_spectrum
  signatures = nmf_res$signatures
  contribution = nmf_res$contribution
  names = sapply(colnames(contribution), function(x) extract_sample_name(x))
  colnames(contribution) = names
  # Reconstruct spectrum by multiplying signature and contribution matrices
  reconstructed = signatures %*% contribution
  # Absolute difference between reconstructed and the original mutation spectrum
  unexplained = abs(reconstructed - mut_spectrum)
  # Total number of mutations not reconstructed
  error = colSums(unexplained)
  # Residual is defined as the fraction of unexplained mutations
  residual = error / colSums(mut_spectrum)
  # Fraction explained by signatures
  expl = 1 - residual
  # Find scaling factor
  factor = expl / colSums(contribution)
  # get new contribution matrix
  for(j in 1:ncol(signatures))
  {
    if(j == 1){new_contribution = factor * contribution[j,]}
    else {new_contribution = rbind(new_contribution, factor * contribution[j,])}
  }
  new_contribution = rbind(new_contribution, residual)
  
  # Assign rownames (depends on no. of signatures)
  for(j in 1:ncol(signatures))
  {
    n = as.character(j)
    name = paste("Signature", n)
    if(j == 1){rownames = name}
    else {rownames = c(rownames, name)}
  }
  rownames = c(rownames, "Residual")
  rownames(new_contribution) = rownames
  
  # Plot contribution including unexplained fraction
  m_contribution = melt(new_contribution)
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  newcolors = c("#7368AC", "#64BFB3","#F9B438", "#BB398D", "#BBD61D", "#706F6F")
  mycolors = newcolors[c(1:ncol(signatures),6)]
  
  plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
    geom_bar(position = "fill", stat="identity") + coord_flip() +  
    xlim(rev(levels(factor(m_contribution$Sample)))) +
    # fill with specified colors  
    scale_fill_manual(values=mycolors) +
    # ylabel
    labs(x = "Sample") +  
    # white background
    theme_bw() +
    # no gridlines
    theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
  return(plot)
}










