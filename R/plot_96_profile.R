#' Plot 96 trinucleotide profile
#'  
#' Plot relative contribution of 96 trinucleotides      
#' @param mut_matrix 96 trinucleotide profile matrix
#' @return 96 trinucleotide profile plot
#' @export


plot_96_profile = function(mut_matrix, colors = c("#DBD7C8", "#B2D39C", "#71C1BA", "#2DAFCE", "#2476B2", "#737E93")){
  # check color vector length
  if(length(colors) != 6){stop("Color vector length not 6")}
  substitutions = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
  index = c(rep(1,1,16), rep(2,1,16), rep(3,1,16), rep(4,1,16), rep(5,1,16), rep(6,1,16))
  # Relative contribution
  norm_mut_matrix = apply(mut_matrix, 2, function(x) x / sum(x) )
  # Context
  context = rownames(norm_mut_matrix)
  # Replace mutated base with dot
  substring(context,2,2) = "."
  # Construct dataframe
  df = data.frame(substitution = substitutions[index], context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  print(ggplot(data=df3, aes(x=context, y=value, fill=substitution, width=0.6)) +  
          geom_bar(stat="identity", colour="black", size=.2) + 
          scale_fill_manual(values=colors) + 
          facet_grid(variable ~ substitution) + 
          ylab("Relative contribution") + 
          coord_cartesian(ylim=c(0,0.15)) +
          scale_y_continuous(breaks=c(0, 0.1)) +
          # no legend
          guides(fill=FALSE) + 
          # white background
          theme_bw() +
          # format text
          theme(axis.title.y=element_text(size=12,vjust=1),
                axis.text.y=element_text(size=8),
                axis.title.x=element_text(size=12),
                axis.text.x=element_text(size=5,angle=90,vjust=0.4),
                strip.text.x=element_text(size=14),
                strip.text.y=element_text(size=14, angle=360),
                panel.grid.major.x = element_blank()))
}