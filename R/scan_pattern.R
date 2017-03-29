#' Genome scan by pattern set
#'
#' @param probs1 object of class \code{\link[qtl2geno]{calc_genoprob}}
#' @param phe data frame with one phenotype
#' @param K kinship matrix
#' @param covar covariate matrix
#' @param map genome map
#' @param patterns data frame of pattern information
#' @param haplos vector of haplotype names
#' @param diplos vector of diplotype names
#' @param condense_patterns remove snp_action from contrasts if TRUE
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{scan_pattern(probs, phe, K, covar, patterns, haplos, diplos)}
#'
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom stringr str_split
#' @importFrom CCSanger sdp_to_pattern
#'
scan_pattern <- function(probs1, phe, K = NULL, covar = NULL,
                         map, patterns, haplos = NULL, diplos = NULL,
                         condense_patterns = TRUE) {
  if(!nrow(patterns))
    return(NULL)

  if(!("contrast" %in% names(patterns)))
    patterns$contrast <- ""

  if(is.null(diplos)) {
    diplos <- dimnames(probs1[[1]])[[2]]
  }
  if(is.null(haplos)) {
    haplos <- unique(unlist(stringr::str_split(diplos, "")))
  }
  if(ncol(phe) > 1) {
    warning("only using first phenotype")
    phe <- phe[, 1, drop=FALSE]
  }

  ## For now, limit to one phenotype.
  ## But see how to have a list across phenotypes
  ## Also need to take care of covariates properly; see qtl2shiny:::scan1_covar.
  pheno_name <- dimnames(phe)[[2]]

  ## SDP patterns
  patterns <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(
        dplyr::filter(patterns,
                      pheno == pheno_name),
        sdp, max_snp, max_pos),
      founders = CCSanger::sdp_to_pattern(sdp, haplos),
      contrast = paste(contrast, collapse=","),
      max_lod = max(max_lod)))
  if(!condense_patterns) {
    dplyr::mutate(patterns,
                  founders = paste(founders, contrast, sep = "_"))
  }
  pattern_three <- pattern_diplos(patterns$sdp, haplos, diplos)
  npat <- nrow(patterns)

  ## Diplotype sets
  dip_set <- sapply(stringr::str_split(rownames(pattern_three), ":"),
                    function(x) {
                      c(x[1], "het", x[2])
                    })
  dimnames(dip_set) <- list(as.character(seq(0, nrow(dip_set) - 1)),
                            patterns$founders)

  # set up first diplotype set
  probs2 <- genoprob_to_patternprob(probs1, pattern_three[1,])
  coefs <- list()
  coefs[[1]] <- qtl2scan::scan1coef(probs2, phe, K, covar)
  scans <- qtl2scan::scan1(probs2, phe, K, covar)
  lod <- matrix(scans, nrow(scans), ncol(dip_set))
  dimnames(lod) <- list(dimnames(scans)[[1]],
                        patterns$founders)

  # loop through other diplotype sets
  # While scans could be combined with cbind method, this seems more efficient.
  if(npat > 1) {
    dimnames(coefs[[1]])[[2]][1:3] <- c("ref","het","alt")
    for(i in seq(2, npat)) {
      probs2 <- genoprob_to_patternprob(probs1, pattern_three[i,])
      coefs[[i]] <- qtl2scan::scan1coef(probs2, phe, K, covar)
      dimnames(coefs[[i]])[[2]][1:3] <- c("ref","het","alt")
      lod[,i] <- qtl2scan::scan1(probs2, phe, K, covar)
    }
  }
  # rearrange patterns by descending max LOD
  patterns$max_pos <- apply(lod, 2,
                            function(x) map[[1]][which.max(x)])
  patterns <- dplyr::arrange(patterns,
                             dplyr::desc(max_lod))

  ## Make sure we have attributes for scans and coefs
  scans <- modify_object(scans, lod[, patterns$founders, drop=FALSE])

  names(coefs) <- patterns$founders
  coefs <- coefs[patterns$founders]
  class(coefs) <- c("listof_scan1coef", class(coefs))

  # return object.
  out <- list(patterns=patterns,
              dip_set = dip_set[, patterns$founders],
              group = as.numeric(pattern_three[patterns$founders,,
                                               drop=FALSE]),
              scan = scans,
              coef = coefs)

  ## Adjust max position from genome scan to SNP scan.
  ## Used for vertical line at max.
  class(out) <- c("scan_pattern", class(out))
  out
}
#' @param object object of class \code{\link{scan_pattern}}
#'
#' @export
#' @method summary scan_pattern
#' @rdname scan_pattern
summary.scan_pattern <- function(object, map, ...) {
  summary(object$coef, scan1_object = object$scan, map, ...)
}
