#' Plot of merge_feature object
#'
#' @param object of class \code{merge_feature}
#' @param pheno name of phenotype to be plotted
#' @param plot_by element to plot by (one of \code{c("pattern","consequence")})
#' @param ... other arguments not used
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 aes facet_wrap geom_jitter ggplot ggtitle xlab ylab
#' @importFrom rlang .data
#' 
#' @rdname merge_feature
#'
ggplot_merge_feature <- function(object, pheno, plot_by=c("pattern","consequence"), ...) {
  haplos <- attr(object, "haplos")
  plot_by <- match.arg(plot_by)
  object$lod <- object[[pheno]]
  object <- dplyr::filter(
    dplyr::mutate(object, pattern = sdp_to_pattern(.data$sdp, haplos)),
    !is.na(.data$lod))
  switch(plot_by,
         pattern = {
           ggplot2::ggplot(object) +
             ggplot2::aes(x = .data$pos, y = .data$lod, col = .data$pattern) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~ .data$snp_type, scale = "free") +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Consequence for", pheno))
         },
         consequence = {
           ggplot2::ggplot(object) +
             ggplot2::aes(x = .data$pos, y = .data$lod, col = .data$snp_type) +
             ggplot2::geom_jitter() +
             ggplot2::facet_wrap(~ .data$pattern, scale = "free") +
             ggplot2::xlab("Position in Mbp") +
             ggplot2::ylab("LOD") +
             ggplot2::ggtitle(paste("Top SNPs by Allele Pattern for", pheno))
         })
}

#' @method autoplot merge_feature
#' @export
#' @export autoplot.merge_feature
#' @rdname merge_feature
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.merge_feature <- function(object, ...)
  ggplot_merge_feature(object, ...)
