#' Plot allele1 object
#' 
#' Plot alleles for haplotype, diplotype and top patterns and genome position.
#' 
#' @param x Object of class \code{\link{allele1}}.
#' @param scan1_object Optional object of class \code{\link[qtl2]{scan1}} to find peak.
#' @param map Genome map (required if \code{scan1_object} present)
#' @param pos Genome position in Mbp (supercedes \code{scan1_object})
#' @param trim If \code{TRUE}, trim extreme alleles.
#' @param frame If \code{TRUE}, enable frames for \code{\link[plotly]{ggplotly}}
#' @param ... Other parameters ignored.
#' 
#' @export
#' @importFrom ggplot2 aes element_blank 
#' facet_grid geom_text ggplot scale_x_continuous theme
#' @importFrom dplyr filter group_by mutate ungroup
#' 
ggplot_allele1 <- function(x, scan1_object=NULL, map=NULL, pos=NULL, trim = TRUE, 
                         frame = FALSE, ...) {
  
  if(is.null(pos)) {
    if(is.null(scan1_object))
      pos_center <- median(x$pos)
    else
      pos_center <- summary(scan1_object, map)$pos[1]
  } else {
    pos_center <- pos
    if(pos_center < min(x$pos) | pos_center > max(x$pos))
      stop("position must be within range of scans")
  }
  
  if(!frame) {
    tmpfn <- function(pos, pos_center) {
      a <- abs(pos - pos_center)
      a == min(a)
    }
    x <- dplyr::ungroup(
      dplyr::filter(
        dplyr::group_by(x, source),
        tmpfn(pos, pos_center)))
  }
  
  if(trim & !attr(x, "blups"))
    x <- trim_quant(x)
  else
    x <- dplyr::mutate(x, trim = effect)
  
  x$x <- jitter(rep(1, nrow(x)))
  
  if(frame) {
    # Need devtools::install_github(“ropensci/plotly”)
    p <- ggplot2::ggplot(x,
           ggplot2::aes(x=x, y=trim, value=effect, col = probe, 
                        frame = pos, ids = allele, type=source)) + 
      ggplot2::geom_point(size = 3, shape = 1)
  } else {
    p <- ggplot2::ggplot(x,
           ggplot2::aes(x=x, y=trim, value=effect, col = probe, 
                        label=allele)) + 
      ggplot2::geom_text(size=4)
  }
  p + ggplot2::facet_grid(~source, scales = "free") +
    ggplot2::ylab("allele means") +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(expand=c(0,0.01))
} 
#' @export
#' @rdname ggplot_allele1
autoplot.allele1 <- function(x, ...)
  ggplot_allele1(x, ...)

trim_quant <- function(object, beyond = 3) {
  quant <- quantile(object$effect, c(.25,.75))
  range <- quant + c(-1,1) * beyond * diff(quant)
  object$trim <- pmin(pmax(object$effect, range[1]), range[2])
  object
}
