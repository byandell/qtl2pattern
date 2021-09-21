#' Collapse genoprob according to pattern
#'
#' @param probs1 object of class \code{\link[qtl2]{calc_genoprob}}
#' @param sdp SNP distribution pattern
#' @param alleles use allele string if \code{TRUE}
#'
#' @return object of class \code{\link[qtl2]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' \donttest{
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Convert genotype probabilities to pattern probabilities for pattern 1.
#' pattern_pr <- genoprob_to_patternprob(pr, 7, TRUE)
#' 
#' str(pr)
#' str(pattern_pr)
#' }
#' 
#' @export
genoprob_to_patternprob <- function(probs1, sdp, alleles = FALSE) {
  
  bprobs <- probs1[[1]]
  
  # SDP to pattern_sets
  diplos <- dimnames(bprobs)[[2]]
  haplos <- unique(unlist(stringr::str_split(diplos, "")))
  pattern_sets <- c(pattern_diplos(sdp, haplos, diplos))
  upat <- as.character(sort(unique(pattern_sets)))
  if(alleles) {
    pattern_alleles <- 
      c(stringr::str_split(sdp_to_pattern(sdp, haplos), ":")[[1]], "het")
    pattern_alleles <- pattern_alleles[c(1,3,2)]
    pattern_sets <- pattern_alleles[1 + pattern_sets]
    upat <- pattern_alleles
  }
    
  # Reduce to array with unique pattern_sets.
  pprobs <- bprobs[, !duplicated(pattern_sets),, drop=FALSE]
  dimnames(pprobs)[[2]] <- upat
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
