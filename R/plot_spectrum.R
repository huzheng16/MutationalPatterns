#' Plot point mutation spectrum
#'    
#' @param type_occurences Type occurences matrix
#' @param CT Optional distinction between C>T at CpG
#' @param by Optional grouping variable
#' @param colors Optional color vector with 7 values
#' @return Spectrum plot
#' @export

plot_spectrum = function(type_occurences, CT = F, by = "all", colors = c("#DBD7C8", "#B2D39C", "#B3D9CF","#71C1BA", "#2DAFCE", "#2476B2", "#737E93"))
{
  # check color vector length
  if(length(colors) != 7){stop("Color vector length not 7")}
  if(CT == F){type_occurences = type_occurences[,1:6] }
  if(CT == T){type_occurences = type_occurences[,c(1:2,8,7,4:6)]}
  # relative contribution per sample
  df2 = type_occurences / rowSums(type_occurences)
  # add by info to df
  df2$by = by
  # reshape
  df3 = melt(df2, id.vars = "by")
  # count number of mutations per mutation type
  counts = melt(type_occurences, measure.vars = colnames(type_occurences))
  df4 = cbind(df3, counts$value)
  colnames(df4)[4] = "nmuts" 
  # calculate the mean and the stdev on the value for each group broken down by type + variable
  x = ddply(df4, c("by", "variable"), summarise, mean = mean(value), stdev = sd(value))
  info_x = ddply(df4, c("by"), summarise, total_individuals = sum(value), total_mutations = sum(nmuts))
  x = merge(x, info_x)
  info_type = data.frame(sub_type = c("C>A", "C>G", "C>T", "C>T", "C>T", "T>A", "T>C", "T>G"), variable = c("C>A", "C>G", "C>T", "C>T at CpG", "C>T other", "T>A", "T>C", "T>G"))
  x = merge(x,info_type)
  # define colors for plotting
  if(CT == F){colors = colors[c(1,2,4:7)]}
  # define positioning of error bars
  x$error_pos = x$mean
  # if C>T stacked bar (distinction between CpG sites and other)
  if(CT == T){
    # adjust positioning of error bars for stacked bars
    # mean of C>T at CpG should be plus the mean of C>T at other
    x = x[order(x$by),]
    CpG = which(x$variable == "C>T at CpG")
    other = which(x$variable == "C>T other")
    x$error_pos[CpG] = x$error_pos[other] + x$error_pos[CpG]
    # define order of bars
    order = order(factor(x$variable, levels = c("C>A", "C>G",  "C>T other", "C>T at CpG","T>A", "T>C", "T>G")))
    x = x[order,]
  }
  # make barplot
  plot = ggplot(data=x, aes(x=sub_type, y=mean, fill=variable, group=sub_type)) +
          geom_bar(stat="identity") +
          geom_errorbar(aes(ymin=error_pos-stdev, ymax=error_pos+stdev), width=0.2) + 
          facet_wrap(by ~ total_mutations) + 
          scale_fill_manual(values=colors, name="Point mutation type") +
          theme_bw() +
          xlab("") +
          ylab("Relative contribution") + 
          theme(axis.ticks = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank()) 
  return(plot)
} 