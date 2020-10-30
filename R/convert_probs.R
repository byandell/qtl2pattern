# Convert probs object to new format.
convert_probs <- function(object) {
  # calc_genoprob attributes: is_x_chr, alleles, alleleprobs, and crosstype

  # Could be feather_genoprob.
  if(inherits(object, "feather_genoprob"))
    return(object)

  # If already an array, assume conversion already done.
  if(is.array(object[[1]]))
    return(object)

  out <- object$probs

  for(obj in c("crosstype", "is_x_chr", "alleles", "alleleprobs")) {
    if(!is.null(object[[obj]]))
      attr(out, obj) <- object[[obj]]
  }
  class(out) <- class(object)
  out
}
