#' Collapse genoprob according to pattern
#'
#' @param probs1 object of class \code{\link[qtl2scan]{calc_genoprob}}
#' @param pattern_sets vector of pattern sets
#'
#' @return object of class \code{\link[qtl2scan]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{genoprob_to_patternprob(probs1, pattern_sets)}
#' @export
genoprob_to_patternprob <- function(probs1, pattern_sets) {
  pprobs <- probs1[[1]][,!duplicated(pattern_sets),,drop=FALSE]
  dimnames(pprobs)[[2]] <- upat <- as.character(sort(unique(pattern_sets)))
  ## For duplicated, need to run through set
  for(i in upat) {
    j <- which(duplicated(pattern_sets))
    j <- j[pattern_sets[j] == i]
    if(length(j)) {
      pprobs[,i,] <- pprobs[,i,] + 
        apply(probs1[[1]][,j,, drop=FALSE], c(1,3), sum)
    }
  }
  modify_aprobs(probs1, pprobs)
}
