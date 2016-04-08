#' Compare two 96 mutation profiles
#' 
#' @return plot
#' @export

plot_compare_profiles = function(profile1, profile2, profile_names = c("profile 1", "profile 2"), profile_ymax = 0.15, diff_ylim = c(-0.02, 0.02), colors = c("#DBD7C8", "#B2D39C", "#71C1BA", "#2DAFCE", "#2476B2", "#737E93"))
{
  s1_relative = profile1 / sum(profile1)
  s2_relative = profile2 / sum(profile2)
  diff = s2_relative - s1_relative
  # residual sum of squares
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = T, digits = 3)
  
  x = cbind(s1_relative, s2_relative, diff)
  colnames(x) = c(profile_names, "Difference")
  
  substitutions = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G')
  index = c(rep(1,1,16), rep(2,1,16), rep(3,1,16), rep(4,1,16), rep(5,1,16), rep(6,1,16))
  # Context
  context = rownames(x)
  # Replace mutated base with dot
  substring(context,2,2) = "."
  # Construct dataframe for plotting
  df = data.frame(substitution = substitutions[index], context = context)
  rownames(x) = NULL
  df2 = cbind(df, as.data.frame(x))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # Add dummy non_visible data points to force y axis limits per facet
  df4 = data.frame(substitution = rep("C>A", 4), context = rep("A.A",4), variable = c(profile_names, "Difference", "Difference"), value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2]))
  
  plot = ggplot(data=df3, aes(x=context, y=value, fill=substitution, width=0.6)) +  
    geom_bar(stat="identity", position = "identity", colour="black", size=.2) + 
    geom_point(data = df4, aes(x = context, y = value), alpha = 0) +
    scale_fill_manual(values=colors) + 
    facet_grid(variable ~ substitution, scales = "free_y") + 
    ylab("Relative contribution") + 
    # ylim(-yrange, yrange) +
    # no legend
    guides(fill=FALSE) + 
    # white background
    theme_bw() +
    ggtitle(paste("RSS =", RSS)) + 
    # format text
    theme(axis.title.y=element_text(size=12,vjust=1),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=12),
          axis.text.x=element_text(size=5,angle=90,vjust=0.4),
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14, angle=360),
          panel.grid.major.x = element_blank())
  return(plot)
}