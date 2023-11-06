#' Plot scan pattern usign ggplot2
#' 
#' @param object object of class \code{\link{scan1pattern}}
#' 
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2]{insert_pseudomarkers}}.
#'
#' @param plot_type type of plot from \code{c("lod","coef")}
#' @param patterns allele patterns to plot (default all)
#' @param columns columns used for coef plot
#' @param min_lod minimum LOD peak for contrast to be retained
#' @param lodcolumn columns used for scan1 plot (default all \code{patterns})
#' @param facet Plot facets if multiple phenotypes and patterns provided (default = \code{"pheno"}).
#' @param ... additional parameters
#'
#' @return object of class \code{\link[ggplot2]{ggplot}}
#' 
#' @export
#' 
#' @importFrom dplyr bind_cols filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 aes geom_path geom_vline ggplot ggtitle
#' @importFrom rlang .data
#' 
#' @rdname scan1pattern
#' 
ggplot_scan1pattern <- function(object, map, plot_type = c("lod","coef","coef_and_lod"),
                              patterns = object$patterns$founders,
                              columns = 1:3,
                              min_lod = 3,
                              lodcolumn = seq_along(patterns),
                              facet = "pheno", ...) {
  plot_type <- match.arg(plot_type)
  
  ggplot_scan1pattern_internal(object, map, plot_type, 
                               patterns, columns, min_lod,
                               lodcolumn, facet, ...)
}
ggplot_scan1pattern_internal <- function(object, map, plot_type, 
                                         patterns, columns, min_lod,
                                         lodcolumn, facet,
                                         colors = NULL, ...) {
  
  object$patterns <- dplyr::filter(object$patterns,
                              .data$max_lod >= min_lod)
  
  m <- object$patterns$founders %in% patterns
  patterns <- object$patterns$founders[m]
  pheno <- object$patterns$pheno[m]
  tmp <- object$scan
  colnames(tmp) <- pheno
  object$scan <- modify_object(object$scan, tmp[, m, drop = FALSE])
  
  pattern <- matrix(patterns, nrow(object$scan), ncol(object$scan), byrow = TRUE)
  
  tmp <- class(object$coef)
  object$coef <- object$coef[m]
  if(length(object$coef) > 1 & length(unique(pheno)) > 1) {
    names(object$coef) <- paste0(pheno, "_", names(object$coef))
  }
  class(object$coef) = tmp
  
  switch(plot_type,
         lod = autoplot(object$scan, map, lodcolumn = lodcolumn,
                        pattern = pattern, 
                        facet = facet, ...),
         coef = autoplot(object$coef, map, columns, colors = colors, ...),
         coef_and_lod = autoplot(object$coef, map, columns, 
                                 scan1_output = object$scan,
                                 lodcolumn = lodcolumn,
                                 pattern_lod = pattern, 
                                 facet_lod = facet, colors = colors, ...))
}

#' @method autoplot scan1pattern
#' @export
#' @export autoplot.scan1pattern
#' @rdname scan1pattern
#' 
#' @importFrom ggplot2 autoplot
#' 
autoplot.scan1pattern <- function(object, ...)
  ggplot_scan1pattern(object, ...)
