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
  tmp <- probs1[[1]][,!duplicated(pattern_sets),,drop=FALSE]
  dimnames(tmp)[[2]] <- upat <- as.character(sort(unique(pattern_sets)))
  ## For duplicated, need to run through set
  for(i in upat) {
    j <- which(duplicated(pattern_sets))
    j <- j[pattern_sets[j] == i]
    if(length(j)) {
      tmp[,i,] <- tmp[,i,] + 
        apply(probs1[[1]][,j,, drop=FALSE], c(1,3), sum)
    }
  }
  probs1[[1]] <- tmp
  attr(probs1, "alleles") <- dimnames(probs1[[1]])[[2]]
  attr(probs1, "alleleprobs") <- TRUE # <- so that genoprob_to_alleleprob leaves unchanged

  probs1
}
#' Collapse genoprob according to pattern
#'
#' @param probs1 object of class \code{\link[qtl2scan]{calc_genoprob}}
#' @param action SNP gene action type
#'
#' @return object of class \code{\link[qtl2scan]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{snpprob_collapse(snpprobs, snpsets)}
#' @export
snpprob_collapse <- function(snpprobs,
                             action = c("additive","add+dom","non-add",
                                        "recessive","dominant")) {
  if(dim(snpprobs[[1]])[2] == 2)
    return(snpprobs)
  if(is.null(action))
    return(snpprobs)
  if(action == "add+dom")
    return(snpprobs)
  if(action == "additive") {
    tmp <- snpprobs[[1]][,c(1,3),,drop=FALSE]
    dimnames(tmp)[[2]] <- as.character(0:1)
    tmp <- tmp + 0.5 * snpprobs[[1]][,c(2,2),,drop=FALSE]
  } else {
    snpsets <- switch(action,
                      "non-add" = c(0,1,0),
                      "recessive" = c(0,1,1),
                      "dominant" = c(0,0,1))

    ## Assume snpsets only has 0s and 1s
    if(all(snpsets==1) | all(snpsets==0) | length(snpsets) != 3)
      stop("need mix of 0s and 1s")

    snpsets <- as.character(snpsets)
    tmp <- snpprobs[[1]][,!duplicated(snpsets),,drop=FALSE]
    dimnames(tmp)[[2]] <- unique(snpsets)
    tmpa <- snpsets[duplicated(snpsets)]
    tmp[,tmpa,] <- tmp[,tmpa,,drop=FALSE] + snpprobs[[1]][,duplicated(snpsets),,drop=FALSE]
  }
  snpprobs[[1]] <- tmp
  attr(snpprobs, "alleles") <- dimnames(snpprobs[[1]])[[2]]
  attr(snpprobs, "alleleprobs") <- TRUE # <- so that genoprob_to_alleleprob leaves unchanged
  snpprobs
}
