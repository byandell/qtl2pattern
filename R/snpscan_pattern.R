#' Plot coefficients by pattern set
#'
#' @param probs1 object of class \code{\link[qtl2scan]{calc_genoprob}}
#' @param phe data frame with one phenotype
#' @param K kinship matrix
#' @param covar covariate matrix
#' @param patternl list of patterns
#' @param patterns data frame of pattern information
#' @param title title of plot (default "")
#' @param col_names names of colors (default from patterns)
#' @param probs2 genotype probabilities for pattern set
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return object of class \code{snpscan_pattern}
#'
#' @importFrom qtl2scan scan1 scan1coef
#'
snpscan_pattern <- function(probs1, phe, K, covar, pattern) {
  if(is.matrix(pattern))
    pattern <- apply(pattern, 1, paste, collapse="")
  pattern <- as.character(pattern)
  tmp <- dim(probs1$probs[[1]])[2]
  if(length(pattern) != tmp)
    stop("pattern must have length ", tmp)

  probs2 = genoprob_to_patternprob(probs1, pattern)
  out <- list(scan = qtl2scan::scan1(probs2, phe, K, covar),
              coef = qtl2scan::scan1coef(probs2, phe, K, covar),
              pattern = pattern)
  class(out) <- c("snpscan_pattern", class(out))
  out
}
#' summary of snpscan_pattern object
#'
#' @param object object of class \code{\link{snpscan_pattern}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return invisible probs collapsed to pattern
#'
#' @method summary snpscan_pattern
#' @rdname snpscan_pattern
#'
summary.snpscan_pattern <- function(object, ...) {
  cat(pattern, "\n")
  print(tmp <- max(out$scan))
  out$coef$coef[rownames(tmp),]
}

#' Plot of snpscan_pattern object
#'
#' @param x object of class \code{\link{snpscan_pattern}}
#' @param title title for plots (default "")
#' @param colors colors for plot (default \code{link[CCSanger]{CCcolors}})
#' @param col_names names of colors (default from patterns)
#' @param color_lod color of LOD plot line
#' @param center center coefficients if \code{TRUE} (default)
#' @param max_coef maximum absolute coefficient for plotting (default 5)
#' @param ... additional parameters
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @return invisible probs collapsed to pattern
#'
#' @method plot snpscan_pattern
#' @rdname snpscan_pattern
#'
plot.snpscan_pattern <- function(x,
                                 title="",
                                 colors = NULL,
                                 col_names = NULL,
                                 color_lod = 1,
                                 coef=TRUE,
                                 lod=TRUE,
                                 center=TRUE,
                                 max_coef=5,
                                 ...) {
  pattern <- as.character(x$pattern)
  upat <- sort(unique(pattern))
  ndim <- length(upat)
  if(is.null(colors))
    colors <- CCSanger::CCcolors[1 + seq(0,ndim-1)%%8]
  if(is.null(col_names))
    col_names <- pattern
  names(colors) <- col_names
  columns <- match(upat, dimnames(x$coef$coef)[[2]])

  # Center coef on mean per locus if TRUE
  if(center) {
    tmp <- x$coef$coef[,upat]
    x$coef$coef[,upat] <- tmp - apply(tmp, 1, mean, na.rm=TRUE)
  }

  ylim <- range(c(x$coef$coef[,upat]), na.rm = TRUE)
  if(any(abs(ylim) > max_coef)) {
    if(ylim[1] < -max_coef)
      ylim[1] <- -max_coef
    if(ylim[2] > max_coef)
      ylim[2] <- max_coef
  }
  # qtl2ggplot::plot_coef
  plot(x$coef, columns,
            col = colors,
            lty = 1 + floor(seq(0, ndim - 1) / 8),
            ylim = ylim, main = title,
            legend.position = "right",
            ...)
}
#' Plot coefficients by pattern set
#'
#' @param probs1 object of class \code{\link[qtl2scan]{calc_genoprob}}
#' @param phe data frame with one phenotype
#' @param K kinship matrix
#' @param covar covariate matrix
#' @param patterns data frame of pattern information
#' @param haplos vector of haplotype names
#' @param diplos vector of diplotype names
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname snpscan_pattern
#'
plot_snpscan <- function(probs1, phe, K, covar,
                         patterns,
                         haplos, diplos) {
  x <- scan_pattern(probs1, phe, K, covar,
                    patterns, haplos, diplos)
  if(!is.null(x)) {
    print(plot(x, "lod"))
    print(plot(x, "coef"))
  }
}
