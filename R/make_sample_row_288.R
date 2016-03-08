#' Make sample row 288
#' 
#' @param signatures Dataframe containing mutational sigantures
#' @export
#' 

make_sample_row_288 = function(vcf)
{
  # overlap mutations
  early_muts = subsetByOverlaps(vcf, GR.early)
  intermediate_muts = subsetByOverlaps(vcf, GR.intermediate)
  late_muts = subsetByOverlaps(vcf, GR.late)
  
  type_context_early = get_type_context(early_muts)
  type_context_intermediate = get_type_context(intermediate_muts)
  type_context_late = get_type_context(late_muts)
  
  # make 96-trinucleotide count vector per set
  early_vector = make_sample_row_96(type_context_early)
  intermediate_vector = make_sample_row_96(type_context_intermediate)
  late_vector = make_sample_row_96(type_context_late)
  
  # add names
  C_triplets = c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT") 
  T_triplets = c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT") 
  triplets = c(rep(C_triplets,3), rep(T_triplets,3))
  names_early = paste(triplets, "-1", sep = "")
  names_intermediate = paste(triplets, "-2", sep = "")
  names_late = paste(triplets, "-3", sep = "")
  
  # combine vectors in alternating fashion
  vector = c(rbind(early_vector, intermediate_vector, late_vector))
  names = c(rbind(names_early, names_intermediate, names_late))
  names(vector) = names
  return(vector)
}

