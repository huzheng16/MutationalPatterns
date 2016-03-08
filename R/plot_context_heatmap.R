#' Plot mutation context heatmap
#' 
#' Creates a heatmap of base substitution occurences with their context, normalized per sample and for trinucloetide occurences in the genome
#' @param vcf_files Character vector containing vcf files
#' @export
#' 

plot_context_heatmap = function(vcf_files, output)
{
  data("trinucleotide_fraction_mat")
  matrix = combined_matrix(vcf_files)
  plot_heatmap(matrix, vcf_files, output)
}

make_sample_row = function(type_context)
{
  types = c('C>A','C>G','C>T','T>A','T>C','T>G')
  matrix = matrix(data = 0, nrow = 4, ncol = 24)
  # for all mutations in this sample
  for(i in 1:length(type_context[[1]]))
  {
    # Find mutation type number
    type_number = which(types == type_context[[1]][i])
    # Find mutation context number
    context = as.character(type_context[[2]][i])
    context = chartr('A,C,G,T', '1,2,3,4', context)
    fiveprime = as.integer(substring(context,1,1))
    threeprime = as.integer(substring(context,3,3))
    x = fiveprime
    y = threeprime + (type_number-1)*4
    matrix[x,y] = matrix[x,y] + 1
  }
  # normalize per sample
  matrix_sam = matrix / sum(matrix)
  # normalize for trinucleotide occurences
  matrix_sam_norm = matrix_sam / trinucleotide_fraction_mat
  # return log transformed matrix
  final = log10(matrix_sam_norm)
  # find smallest value after -Inf
  min = min(final[final > -Inf])
  # replace -Inf with smallest value / 2
  final[final == -Inf] = min/2
  return(final)
}

combined_matrix = function(vcf_files)
{
  count = 0
  for(vcf_file in vcf_files)
  {
    vcf = read_vcf(vcf_file)
    type_context = get_type_context(vcf)
    matrix = make_sample_row(type_context)
    count = count + 1
    print(count)
    if(count == 1)
    {
      output_matrix = matrix
    }
    else
    {
      output_matrix = rbind(output_matrix, matrix)
    }
  }
  return(output_matrix)
}

plot_heatmap = function(matrix, vcf_files, output)
{
  # number of samples
  n_samples = length(vcf_files)
  # extract sample names from file names
  samples =  sapply(vcf_files, function(x) tail(strsplit(x, "/")[[1]],1))
  # plotting somehow only works if I use global variables
  NAMES <<- sapply(samples, function(x) strsplit(x, ".vcf"))
  # 6 base substitution types  
  TYPES <<- c('C>A','C>G','C>T','T>A','T>C','T>G')
  # create col and row names
  colnames = rep(c('A','C','G','T'), 6)
  rownames = rep(c('A','C','G','T'), n_samples)
  # determine position of horizontal and vertical lines
  HLINES <<- seq(0.5, 4 * n_samples + 0.5, 4) 
  VLINES <<- seq(0.5, 24.5, 4)
  # create colour palette
  jBuPuFun = colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
  paletteSize = 256
  jBuPuPalette = jBuPuFun(paletteSize)
  # determine position for substitution types
  POS <<- seq(2.5, 22.5, 4)
  # determine position for sample names
  POS2 <<- seq(n_samples * 4 - 1.5, 2.5, -4)
  # output pdf
  pdf(output)
  # plot heatmap
  heatmap.2(as.matrix(matrix), key = T, Rowv = NA, Colv = NA, scale = 'none', dendrogram = 'none', 
            density.info = 'none', col = jBuPuPalette, trace = 'none', keysize = 1,
            labCol = colnames, labRow = rownames, cexRow = 0.6, cexCol = 0.6, srtCol = 0, 
            add.expr = {abline(h=HLINES, v=VLINES, lwd = 1); mtext(text = TYPES, side = 3, at = POS, cex = 0.6); mtext(text = NAMES, side = 2, at = POS2, adj = 1, las = 1, line = 1, cex = 0.6)}, 
            ylab = "5' base", xlab = "3' base", cex = 0.6)
  dev.off()
  # remove global variables
  rm(HLINES, VLINES, NAMES, POS, POS2, TYPES, pos = ".GlobalEnv", inherits = T)
}