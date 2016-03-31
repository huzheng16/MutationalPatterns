#' Count 96 trinucleotide mutation occurences
#'  
#' @param type_context result from get_type_context function
#' @return vector with 96 trinucleotide mutation occurences
#' @export

mut_96_occurences = function(type_context)
{
  types = c('C>A','C>G','C>T','T>A','T>C','T>G')
  C_triplets = c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT") 
  T_triplets = c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT") 
  triplets = c(rep(C_triplets,3), rep(T_triplets,3))
  vector = rep(0,96)
  names(vector) = triplets
  
  # for all mutations in this sample
  for(i in 1:length(type_context[[1]]))
  {
    # Find mutation type
    type = which(types == type_context[[1]][i])
    # Find triplet
    if(type < 4)
    {
      context = which(C_triplets == type_context[[2]][i])
    } else
    {
      context = which(T_triplets == type_context[[2]][i])
    }   
    pos = (type - 1)*16 + context
    vector[pos] = vector[pos] + 1
  }
  return(vector)
}