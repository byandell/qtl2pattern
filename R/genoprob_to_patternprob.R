#' Collapse genoprob according to pattern
#'
#' @param probs1 object of class \code{\link[qtl2]{calc_genoprob}}
#' @param sdp SNP distribution pattern
#'
#' @return object of class \code{\link[qtl2]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \donttest{genoprob_to_patternprob(probs1, pattern_sets)}
#' @export
genoprob_to_patternprob <- function(probs1, sdp) {
  
  bprobs <- probs1[[1]]
  
  # SDP to pattern_sets
  diplos <- dimnames(bprobs)[[2]]
  haplos <- unique(unlist(stringr::str_split(diplos, "")))
  pattern_sets <- c(pattern_diplos(sdp, haplos, diplos))
    
  # Reduce to array with unique pattern_sets.
  pprobs <- bprobs[, !duplicated(pattern_sets),, drop=FALSE]
  dimnames(pprobs)[[2]] <- upat <- as.character(sort(unique(pattern_sets)))
  ## For duplicated, need to run through set
  for(i in upat) {
    j <- which(duplicated(pattern_sets))
    j <- j[pattern_sets[j] == i]
    if(length(j)) {
      pprobs[,i,] <- pprobs[,i,] + 
        apply(bprobs[, j,, drop=FALSE], c(1,3), sum)
    }
  }
  modify_aprobs(probs1, pprobs)
}
