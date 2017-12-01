#' Collapse genoprob according to pattern
#'
#' @param probs1 object of class \code{\link[qtl2]{calc_genoprob}}
#' @param action SNP gene action type
#'
#' @return object of class \code{\link[qtl2]{calc_genoprob}}
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
    patprobs <- snpprobs[[1]][,c(1,3),,drop=FALSE]
    dimnames(patprobs)[[2]] <- as.character(0:1)
    patprobs <- patprobs + 0.5 * snpprobs[[1]][,c(2,2),,drop=FALSE]
  } else {
    snpsets <- switch(action,
                      "non-add" = c(0,1,0),
                      "recessive" = c(0,1,1),
                      "dominant" = c(0,0,1))
    
    ## Assume snpsets only has 0s and 1s
    if(all(snpsets==1) | all(snpsets==0) | length(snpsets) != 3)
      stop("need mix of 0s and 1s")
    
    snpsets <- as.character(snpsets)
    patprobs <- snpprobs[[1]][, !duplicated(snpsets),, drop = FALSE]
    dimnames(patprobs)[[2]] <- unique(snpsets)
    patprobsa <- snpsets[duplicated(snpsets)]
    patprobs[,patprobsa,] <- patprobs[, patprobsa,, drop = FALSE] + 
      snpprobs[[1]][, duplicated(snpsets),, drop = FALSE]
  }
  modify_aprobs(snpprobs, patprobs)
}
