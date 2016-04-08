#' Plot signature contribution
#' 
#' Plot contribution of signatures including unexplained fraction
#' 
#' @param nmf_res NMF result
#' @param index optional sample subset parameter
#' @param coord_flip Flip X and Y coordinates, default = F
#' @export
#' 

plot_contribution = function(contribution, index=c(), coord_flip = F)
{
  # optional subsetting
  if(length(index > 0)){contribution = contribution[,index]}
  # Plot contribution
  m_contribution = melt(contribution)
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  
  plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
    geom_bar(position = "fill", stat="identity", colour="black")  +  
    # make sure sample ordering is correct
    xlim(rev(levels(factor(m_contribution$Sample)))) +
    # ylabel
    labs(x = "", y = "Relative contribution") +  
    scale_fill_discrete(name="Signature") +
    # white background
    theme_bw() +
    # no gridlines
    theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank())
  
  # optional coordinate flipping
  if(coord_flip == T){plot = plot + coord_flip()}
  return(plot)
}

