#' Genome scan by pattern set
#'
#' @param probs1 object of class \code{\link[qtl2scan]{calc_genoprob}}
#' @param phe data frame with one phenotype
#' @param K kinship matrix
#' @param covar covariate matrix
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
scan_pattern <- function(probs1, phe, K, covar,
                         patterns, haplos = NULL, diplos = NULL,
                         condense_patterns = TRUE) {
  if(!nrow(patterns))
    return(NULL)

  if(!("contrast" %in% names(patterns)))
    patterns$contrast <- ""

  if(is.null(diplos)) {
    diplos <- dimnames(probs1$probs[[1]])[[2]]
  }
  if(is.null(haplos)) {
    haplos <- unique(unlist(stringr::str_split(diplos, "")))
  }

  patterns <- dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(patterns,
                      pheno, sdp, max_snp, max_pos),
      contrast = paste(CCSanger::sdp_to_pattern(sdp, haplos),
                       paste(contrast, collapse=","),
                       sep = "_"),
      max_lod = max(max_lod)))
  pattern_three <- pattern_diplos(patterns$sdp, haplos, diplos)
  dip_set <- sapply(stringr::str_split(rownames(pattern_three), ":"),
                    function(x) {
                      c(x[1], "het", x[2])
                    })
  ## Colors from http://colorbrewer2.org qualitative
  out <- list(patterns=patterns,
              dip_set = dip_set,
              scans=list())
  for(i in seq_len(ncol(dip_set))) {
    tmp <- snpscan_pattern(probs1, phe, K, covar,
                           pattern_three[i,])
    dimnames(tmp$coef$coef)[[2]][seq(nrow(dip_set))] <-
      dip_set[,i]
    out$scans[[i]] <- tmp

  }
  dimnames(out$dip_set)[[1]] <- as.character(seq(0, nrow(out$dip_set)-1))
  dimnames(out$dip_set)[[2]] <- names(out$scans) <- patterns$contrast

  out$patterns$max_pos <- sapply(out$scans, function(x)
    max(x$scan)$pos)
  out$patterns <- dplyr::arrange(out$patterns,
                                 dplyr::desc(max_lod))
  contrasts <- out$pattern$contrast

  # For now, refactor this stuff. Later, get rid of snpscan_pattern?
  out$group <- as.matrix(sapply(out$scans, function(x) x$pattern))
  out$group <- out$group[, contrasts, drop=FALSE]

  coefs <- lapply(out$scans, function(x) {
    xcoef <- x$coef
    dimnames(xcoef$coef)[[2]][1:3] <- c("ref", "het", "other")
    xcoef
  })
  coefs <- coefs[contrasts]
  class(coefs) <- c("listof_scan1coef", class(coefs))
  out$coef <- coefs

  tmplod <- out$scans[[1]]$scan
  tmplod$lod <- sapply(out$scans, function(x) x$scan$lod)
  tmplod$lod <- tmplod$lod[, contrasts, drop=FALSE]
  out$scan <- tmplod

  out$scans <- NULL

  if(condense_patterns) {
    short <- stringr::str_replace(out$patterns$contrast, "_.*", "")
    out$patterns$contrast <- names(out$coef) <- dimnames(out$scan$lod)[[2]] <- short
  }

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
summary.scan_pattern <- function(object, ...) {
  object$patterns
}
#' @param x object of class \code{\link{scan_pattern}}
#' @param plot_type type of plot from \code{c("lod","coef")}
#' @param patterns allele patterns to plot (default all)
#' @param columns columns used for coef plot
#' @param min_lod minimum LOD peak for contrast to be retained
#' @param ylim_coef vertical limits for coef plot
#' @param ... additional parameters
#'
#' @export
#' @method plot scan_pattern
#' @rdname scan_pattern
#' @importFrom dplyr bind_cols filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes geom_path geom_vline ggplot ggtitle
plot.scan_pattern <- function(x, plot_type=c("lod","coef","coef_and_lod"),
                              patterns=x$patterns$contrast,
                              columns = 1:3,
                              min_lod = 3,
                              ylim_coef = c(-2,2),
                              ...) {
  plot_type <- match.arg(plot_type)

  x$patterns <- dplyr::filter(x$patterns,
                              max_lod >= min_lod)

  patterns <- patterns[patterns %in% x$patterns$contrast]
  x$scan$lod <- x$scan$lod[,patterns, drop=FALSE]
  tmp <- class(x$coef)
  x$coef <- x$coef[patterns]
  class(x$coef) = tmp

  switch(plot_type,
         lod = plot(x$scan, lodcolumn = seq_along(patterns), ...),
         coef = plot(x$coef, columns, ylim = ylim_coef, ...),
         coef_and_lod = plot(x$coef, columns, scan1_output = x$scan))
}
