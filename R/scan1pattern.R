#' Genome scan by pattern set
#'
#' @param probs1 object of class \code{\link[qtl2]{calc_genoprob}}
#' @param phe data frame with one phenotype
#' @param K kinship matrix
#' @param covar covariate matrix
#' @param map genome map
#' @param patterns data frame of pattern information
#' @param condense_patterns remove snp_action from contrasts if TRUE
#' @param blups Create BLUPs if \code{TRUE}
#' @param do_scans Do scans if \code{TRUE}.
#'
#' @return List containing:
#' \itemize{
#' \item{patterns} Data frame of summary for top patterns (column \code{founders} has pattern)
#' \item{dip_set} Diplotype sets for contrasts
#' \item{group} Group for each founder pattern
#' \item{scan} Object of class \code{\link[qtl2]{scan1}}.
#' \item{coef} Object of class \code{listof_scan1coef}. See package 'qtl2ggplot'.
#' }
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
#' # Download SNP info for DOex from web and read as RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Convert to SNP probabilities
#' snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
#' snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
#' 
#' # Scan SNPs
#' scan_snppr <- qtl2::scan1(snppr, DOex$pheno)
#' top_snps_tbl <- top_snps_pattern(scan_snppr, snpinfo)
#' 
#' # Summarize to find top patterns
#' patterns <- dplyr::arrange(summary(top_snps_tbl), dplyr::desc(max_lod))
#' 
#' # Scan using patterns.
#' scan_pat <- scan1pattern(pr, DOex$pheno, map = DOex$pmap, patterns = patterns)
#' 
#' # Summary of scan1pattern.
#' summary(scan_pat, DOex$pmap)
#' }
#'
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom stringr str_split
#' @importFrom rlang .data
#'
scan1pattern <- function(probs1, phe, K = NULL, covar = NULL,
                         map, patterns,
                         condense_patterns = TRUE,
                         blups = FALSE,
                         do_scans = TRUE) {
  if(!nrow(patterns))
    return(NULL)

  ## For now, limit to one phenotype.
  ## But see how to have a list across phenotypes
  ## Also need to take care of covariates properly; see scan1_covar.
  pheno_names <- colnames(phe)

  diplos <- dimnames(probs1[[1]])[[2]]
  haplos <- unique(unlist(stringr::str_split(diplos, "")))
  
  if(!("contrast" %in% names(patterns))) 
    patterns$contrast <- ""
  
  ## SDP patterns
  patterns <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(
        dplyr::filter(patterns,
                      .data$pheno %in% pheno_names),
        .data$sdp, .data$snp_id, .data$max_pos, .data$pheno),
      founders = sdp_to_pattern(.data$sdp, haplos),
      contrast = paste(.data$contrast, collapse=","),
      max_lod = max(.data$max_lod)))
  
  if(!condense_patterns & !all(patterns$contrast == "")) {
    dplyr::mutate(patterns,
                  founders = paste(.data$founders, .data$contrast, sep = "_"))
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
  probs2 <- genoprob_to_patternprob(probs1, patterns$sdp[1])
  scan1fn <- ifelse(blups, 
                    qtl2::scan1blup, 
                    qtl2::scan1coef)
  coefs <- list()
  coefs[[1]] <- scan1fn(probs2, 
                        phe[, patterns$pheno[1], drop = FALSE],
                        K, covar)
  if(do_scans) {
    scans <- qtl2::scan1(probs2, 
                             phe[, patterns$pheno[1], drop = FALSE],
                             K, covar)
    lod <- matrix(scans, nrow(scans), ncol(dip_set))
    dimnames(lod) <- list(rownames(scans),
                          patterns$founders)
  }

  # loop through other diplotype sets
  # While scans could be combined with cbind method, this seems more efficient.
  dimnames(coefs[[1]])[[2]][1:3] <- c("ref","het","alt")
  if(npat > 1) {
    for(i in seq(2, npat)) {
      probs2 <- genoprob_to_patternprob(probs1, patterns$sdp[i])
      coefs[[i]] <- scan1fn(probs2, 
                            phe[, patterns$pheno[i], drop = FALSE],
                            K, covar)
      dimnames(coefs[[i]])[[2]][1:3] <- c("ref","het","alt")
      if(do_scans)
        lod[,i] <- qtl2::scan1(probs2, 
                                   phe[, patterns$pheno[i], drop = FALSE],
                                   K, covar)
    }
  }
  if(do_scans) {
    # rearrange patterns by descending max LOD
    patterns$max_pos <- apply(lod, 2,
                              function(x) map[[1]][which.max(x)])
    patterns <- dplyr::arrange(patterns,
                               dplyr::desc(.data$max_lod))
  }

  ## Make sure we have attributes for scans and coefs
  if(do_scans)
    scans <- modify_object(scans, lod[, patterns$founders, drop=FALSE])
  else
    scans <- NULL

  names(coefs) <- patterns$founders
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
  attr(out, "blups") <- blups
  class(out) <- c("scan1pattern", class(out))
  out
}
#' @param object object of class \code{\link{scan1pattern}}
#' @param ... additional parameters passed on to other methods
#' @export
#' @method summary scan1pattern
#' @rdname scan1pattern
summary.scan1pattern <- function(object, map, ...) {
  if(exists("summary_listof_scan1coef")) {
    # Only available if qtl2ggplot package is attached
    # Set up unique names as pheno_pattern_contrast
    pheno <- paste(object$patterns$pheno, object$patterns$founders, sep = "_")
    names(object$coef) <- pheno
    colnames(object$scan) <- pheno
    summary(object$coef, scan1_object = object$scan, map, ...)
  } else {
    object$patterns
  }
}
