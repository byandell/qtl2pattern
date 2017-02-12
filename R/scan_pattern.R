#' Genome scan by pattern set
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
#' @examples
#' \dontrun{scan_pattern(probs, phe, K, covar, patterns, haplos, diplos)}
#'
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom stringr str_split
#' @importFrom CCSanger sdp_to_pattern
#'
scan_pattern <- function(probs1, phe, K, covar,
                         patterns, haplos = NULL, diplos = NULL) {
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

  ## Adjust max position from genome scan to SNP scan.
  ## Used for vertical line at max.
  out$patterns$max_pos <- sapply(out$scans, function(x)
    max(x$scan)$pos)
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
#' @param title title for plot
#' @param ... additional parameters
#'
#' @export
#' @method plot scan_pattern
#' @rdname scan_pattern
#' @importFrom dplyr bind_cols filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes geom_path geom_vline ggplot ggtitle
plot.scan_pattern <- function(x, plot_type=c("lod","coef"),
                              patterns=x$patterns$contrast,
                              title = NULL,
                              ...) {
  plot_type <- match.arg(plot_type)
  switch(plot_type,
         lod = {
           ## LOD scans

           ## bind lod scans across patterns
           lod_scans <- dplyr::bind_cols(lapply(x$scans[patterns],
                                                function(x)
                                                  as.data.frame(x$scan$lod)))
           if(!nrow(lod_scans))
             return(NULL)
           names(lod_scans) <- patterns
           lod_scans$pos <- x$scans[[1]]$scan$map[[1]]
           ## gather by pattern
           lod_scans <- tidyr::gather(lod_scans, contrast,lod,-pos)

           p <- list()
           for(i in patterns) {
             if(is.null(title))
               mytitle <- paste((
                 dplyr::select(
                   dplyr::filter(x$patterns, contrast == i),
                   pheno, max_snp, contrast))[1,],
                 collapse = " ")
             else
               mytitle <- title
             xint <- dplyr::filter(x$patterns, contrast == i)$max_pos
             p[[i]] <- ggplot2::ggplot(
               dplyr::filter(lod_scans, contrast == i),
               ggplot2::aes(pos,lod)) +
               ggplot2::geom_path(lwd=1) +
               ggplot2::geom_vline(xintercept=xint, lty="dashed") +
               ggplot2::ggtitle(mytitle)
           }
           if(length(patterns) == 1)
             p[[1]]
           else {
             for(i in patterns)
               print(p[[i]])
           }
         },
         coef = {
           ## Coefficient scans
           p <- list()
           for(i in patterns) {
             if(is.null(title))
               mytitle <- paste((
                 dplyr::select(
                   dplyr::filter(x$patterns, contrast == i),
                   pheno, max_snp, contrast))[1,],
                 collapse = " ")
             else
               mytitle <- title
             xint <- (dplyr::filter(x$patterns, contrast == i))$max_pos
             p[[i]] <- plot(x$scans[[i]]$coef, 1:3,
                            col = c("#1b9e77","#d95f02","#7570b3"),
                            title = mytitle,
                            ...) +
               ggplot2::geom_vline(xintercept = xint, linetype=2)
           }
           if(length(patterns) == 1)
             p[[1]]
           else {
             for(i in patterns)
               print(p[[i]])
           }
         })
}
