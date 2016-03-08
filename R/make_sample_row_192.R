#' Make sample row 192 for NMF
#' 
#' @param vcf Collapsed vcf object
#' @return vector 
#' @import IRanges
#' @export

make_sample_row_192 = function(vcf)
{
  # load gene bodies file
  data("gene_bodies") 
  # make Granges list object from gene_bodies data
  GR.gene_bodies = GRanges(seqnames = paste("chr", gene_bodies[,1], sep="" ), ranges = IRanges(start = gene_bodies[,2], end = gene_bodies[,3]), 
                           strand = gene_bodies[,6], gene = gene_bodies[,4])
  # merge overlapping genes on same strand, gene names are lost
  GR.gene_bodies = reduce(GR.gene_bodies)
  # make Granges object from vcf
  GR.vcf = rowRanges(vcf)
  # overlap mutations with gene bodies
  overlap = findOverlaps(GR.vcf, GR.gene_bodies)
  overlap = as.data.frame(as.matrix(overlap))
  colnames(overlap) = c('vcf_id', 'gene_body_id')
  # find reference allele of mutations (+ strand of reference genome is reported in vcf file)
  vcf_overlap = vcf[overlap$vcf_id]
  ref = as.character(ref(vcf_overlap))
  # Find the strand of C or T (since we regard mutations as C>X or T>X)
  # which mutations have ref allele C or T
  i = which(ref == "C" | ref == "T")
  # store mutation strand info in vector
  strand_muts = rep(0,nrow(overlap))
  strand_muts[i] = "+"
  strand_muts[-i] = "-"
  # find strand of gene bodies
  strand_genebodies = as.character(strand(GR.gene_bodies)[overlap$gene_body_id])
  # find if mut and gene_bodie are on the same strand
  same_strand = (strand_muts  == strand_genebodies)
  # subset vcf object for both untranscribed and transcribed
  U_index = which(same_strand == TRUE)
  T_index = which(same_strand == FALSE)
  vcf_overlap_U = vcf_overlap[U_index]
  vcf_overlap_T = vcf_overlap[T_index]
  # get type context for both vcf subsets
  type_context_U = get_type_context(vcf_overlap_U)
  type_context_T = get_type_context(vcf_overlap_T)
  # make 96-trinucleotide count vector per set
  U_vector = make_sample_row_96(type_context_U)
  T_vector = make_sample_row_96(type_context_T)
  # add names
  C_triplets = c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT") 
  T_triplets = c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT") 
  triplets = c(rep(C_triplets,3), rep(T_triplets,3))
  names_U = paste(triplets, "-u", sep = "")
  names_T = paste(triplets, "-t", sep = "")
  # combine vectors in alternating fashion
  vector = c(rbind(U_vector, T_vector))
  names = c(rbind(names_U, names_T))
  names(vector) = names
  return(vector)
}










