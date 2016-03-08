#' Plot 192 trinucleotide signatures 
#' 
#' @param signatures Dataframe containing mutational sigantures
#' @export

plot_signatures_192 = function(signatures)
{
  # Variables for signature plotting
  substitutions = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
  index = c(rep(1,1,32), rep(2,1,32), rep(3,1,32), rep(4,1,32), rep(5,1,32), rep(6,1,32))
  # Relative contribution
  norm_signatures = apply(signatures, 2, function(x) x / sum(x) )
  # Signature colnames
  names = c()
  for(i in 1:ncol(signatures)){names = c(names,(paste("Signature", as.character(i))))}
  colnames(norm_signatures) = names
  # Context
  context = substring(rownames(norm_signatures),1,3)
  # Replace mutated base with dot
  substring(context,2,2) = "."
  # Add strand information
  Untranscribed = rep("U", 96)
  Transcribed = rep("T", 96)
  strand = c(rbind(Untranscribed, Transcribed))
  # Construct dataframe
  df = data.frame(substitution = substitutions[index], context = context, strand = strand)
  rownames(norm_signatures) = NULL
  df2 = cbind(df, as.data.frame(norm_signatures))
  df3 = melt(df2, id.vars = c("substitution", "context", "strand"))
  # plotting signature
  print("Plotting signatures...")
  plot = ggplot(data=df3, aes(x=context, y=value, fill=substitution, alpha=strand, width=0.6)) +  
          geom_bar(stat="identity", colour="black", size=.2) + 
          scale_fill_manual(values=sub_colors) +
          scale_alpha_manual(values=c(0.2, 1, 1)) + 
          facet_grid(variable ~ substitution) + 
          ylab("Relative contribution") + 
          coord_cartesian(ylim=c(0,0.2)) +
          scale_y_continuous(breaks=c(0.05, 0.1, 0.15, 0.2)) +
          # no legend
          guides(fill=FALSE) + 
          # white background
          theme_bw() +
          # format text
          theme(axis.title.y=element_text(size=14,vjust=1),
                axis.text.y=element_text(size=14),
                axis.title.x=element_text(size=14),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=14),
                strip.text.y=element_text(size=14),
                panel.grid.major.x = element_blank())
  return(plot)
}

