#' Genomic rainfall plot
#' 
#' Rainfall plot visualizes the types of mutations and intermutation distance
#' @param vcf CollapsedVCF object
#' @param ref_genome BSgenome reference genome object
#' @param title Character string optional plot title
#' @param color Vector of 6 colors used for plotting
#' @param chrom Vector of chromosome/contig names in the reference genome
#' @param cex Dot size
#' @param cex_text Text size
#' @param ylim Maximum y value (genomic distance)
#' @return plot Rainfall plot
#' @export

rainfall_plot = function(vcf, ref_genome, chromosomes, title = "", color = c('#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#E6AB02', 'black'), cex = 2.5, cex_text = 3, ylim = 1e+08)
{
  # check color vector length
  if(length(color) != 6){stop("Color vector length not 6")}
  # get chromsome lengths of reference genome
  chr_length = seqlengths(get(ref_genome))
  # subset
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  # cumulative sum
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  names(chr_cum) = names(chr_length)
  labels = gsub("chr", "", names(chr_length))
  # position of chromosome labels
  m=c()
  for(i in 2:length(chr_cum)){m = c(m,(chr_cum[i-1] + chr_cum[i]) / 2)}
  # mutations characteristics
  type = loc = dist = chrom = c()
  # for each chromosome
  for(i in 1:length(chromosomes))
  {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    type = c(type, get_types(chr_subset)[-1])
    loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    n = dim(chr_subset)[1]-1
    if(n<1){n=0}
    chrom = c(chrom, rep(chr[i],n))
  }
  data = data.frame(type = type, location = loc, distance = dist, chromosome = chrom)
  
  plot = ggplot(data, aes(x=location, y=distance)) +
    geom_point(aes(colour=factor(type)), cex=cex) + 
    geom_vline(xintercept = as.vector(chr_cum), linetype="dotted") +
    annotate("text", x = m, y = ylim, label = labels, cex=cex_text) +
    xlab("Genomic Location") +
    ylab("Genomic Distance") +
    scale_y_log10() +
    scale_colour_manual(values=color) +
    scale_x_continuous(expand = c(0,0), limits=c(0, max(chr_cum))) +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position="bottom", legend.title=element_blank(), legend.key = element_blank(), panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank(), axis.ticks.x =element_blank(), axis.text.x = element_blank()) + 
    guides(colour = guide_legend(nrow = 1))
  
  return(plot)
}



