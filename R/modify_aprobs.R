# Preserve and modify probs object attributes. 
modify_aprobs <- function(probs, new_probs) {
  x_attr <- attributes(probs)
  x_class <- class(probs)
  attrs <- names(x_attr)
  attrs <- attrs[!(attrs %in% c("class", "names", "dim", "dimnames", "alleles", "alleleprobs"))]
  
  new_probs <- list(new_probs)
  names(new_probs) <- names(probs)[1]
  
  for(obj in attrs) {
    attr(new_probs, obj) <- x_attr[[obj]]
  }
  # Override previous values to state these act like allele probs.
  attr(new_probs, "alleles") <- dimnames(new_probs[[1]])[[2]]
  attr(new_probs, "alleleprobs") <- TRUE
  
  class(new_probs) <- c("calc_genoprob", "list")
  
  new_probs
}
