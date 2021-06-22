#' Extract pattern of diplotypes
#'
#' @param sdp vector of straind distribution patterns from \code{\link{top_snps_pattern}}
#' @param haplos vector of haplotype names
#' @param diplos vector of diplotype names
#' @param cont vector of types of contrasts (\code{NULL} or from \code{c("add","dom","b6r","b6d")})
#'
#' @return matrix of diplotype patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @importFrom stringr str_detect
#'
pattern_diplos <- function(sdp, haplos, diplos, cont=NULL) {
  if(!length(sdp))
    return(NULL)
  hap_logic <- sdp_to_logical(sdp, haplos)

  if(is.null(cont))
    cont <- rep("add", length(sdp))
  ## Want to insert logic for contrast type here.
  hap_to_dip_logic <- function(hap_set, cont="add", haplos, diplos) {
    hap <- paste(haplos[hap_set], collapse = "|")
    hap0 <- paste(haplos[!hap_set], collapse = "|")
    str_hap <- stringr::str_detect(diplos, hap)
    str_hap0 <- stringr::str_detect(diplos, hap0)
    switch(cont,
           { #add
             ## Number of copies of 1 haplotype
             (str_hap & str_hap0) + 2 * (str_hap & !str_hap0)
           },
           dom =,
           ovd = str_hap & str_hap0,
           b6r = str_hap,
           b6d = !str_hap0)
  }
  pattern <- sdp_to_pattern(sdp, haplos)
  patternl <- matrix(0, length(sdp), length(diplos),
                     dimnames = list(pattern, diplos))
  for(i in seq_along(pattern)) {
    patternl[i,] <- hap_to_dip_logic(hap_logic[,i], cont[i], haplos, diplos)
  }
  patternl
}
#' Extract pattern of haplotypes
#'
#' @param sdp vector of sdp from \code{\link{top_snps_pattern}}
#' @param haplos vector of haplotype names
#'
#' @return matrix of haplotype patterns
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname pattern_diplos
#'
pattern_haplos <- function(sdp, haplos) {
  if(!length(sdp))
    return(NULL)
  hap_logic <- 1 * t(sdp_to_logical(sdp, haplos))

  pattern <- sdp_to_pattern(sdp, haplos)
  dimnames(hap_logic) <- list(pattern, haplos)
  hap_logic
}
