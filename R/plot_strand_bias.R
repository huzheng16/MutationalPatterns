#' Plot strand bias per base substitution type
#' 
#' @description For each base substitution type and transcriptional strand the total number of mutations
#' and the relative contribution within a group is returned
#' @param strand_bias_df Data.frame, result from strand_bias function
#' @param mode Either "absolute" for absolute number of mutations, or "relative" for relative contribution, default = "relative"
#' @param colors Optional color vector for plotting with 6 values
#' @return Barplot
#' @export


plot_strand_bias = function(strand_bias, colors)
{
  # if colors parameter not provided, set to default colors
  if(missing(colors)){colors=COLORS6}
  # determine max y value for plotting
  # = log2 ratio with pseudo counts
  max = round( max(abs(log2(strand_bias$ratio))), digits = 1) + 0.1
  # plot strand bias with poisson test results
  plot = ggplot(strand_bias, aes(x = type, y = log2(ratio), fill = type)) +
    scale_fill_manual(values = COLORS6) +
    geom_bar(colour = "black", stat ="identity", position = "identity") +
    scale_y_continuous(limits = c(-max, max), breaks = seq(-1, 1, 0.2)) +
    geom_text(aes(x = type, y = log2(ratio), ymax = log2(ratio), 
                  label = significant, vjust = ifelse(sign(log2(ratio)) > 0, 0.5, 1)), 
              size = 8, position = position_dodge(width = 1)) +
    facet_grid(. ~ group) +
    theme_bw()  +
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), legend.title = element_blank()) +
    xlab("") + 
    ylab("log2(transcribed/untranscribed)") +
    scale_x_discrete(breaks=NULL)
  return(plot)
}