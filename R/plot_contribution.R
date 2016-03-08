#' Plot signature contribution
#' 
#' Plot contribution of signatures including unexplained fraction
#' 
#' @param nmf_res NMF result
#' @param index optional subset parameter
#' @export
#' 

plot_contribution = function(nmf_res, index=c())
{
  signatures = nmf_res$signatures
  contribution = nmf_res$contribution
  if(length(index > 0)){contribution = contribution[,index]}
  
  # Assign rownames (depends on no. of signatures)
  for(j in 1:ncol(signatures))
  {
    n = as.character(j)
    name = paste("Signature", n)
    if(j == 1){rownames = name}
    else {rownames = c(rownames, name)}
  }
  rownames(contribution) = rownames
  
  # Plot contribution including unexplained fraction
  m_contribution = melt(contribution)
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  newcolors = c("#7368AC", "#64BFB3","#F9B438", "#BB398D", "#BBD61D", "#706F6F")
  
  mycolors = newcolors[c(1:ncol(signatures),6)]
  
  plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
    geom_bar(position = "fill", stat="identity")  +  
    xlim(rev(levels(factor(m_contribution$Sample)))) +
    # fill with specified colors  
    scale_fill_manual(values=mycolors) +
    # ylabel
    labs(x = "", y = "Relative contribution") +  
    # white background
    theme_bw() +
    # no gridlines
    theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
  return(plot)
}

