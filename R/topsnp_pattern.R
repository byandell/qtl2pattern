#' Fine Mapping scans by allele pattern
#'
#' Separate fine mapping scans by allele pattern.
#'
#' @param out_lmm_snps output of linear mixed model for \code{phename} (see \code{\link[qtl2scan]{scan1}})
#' @param proband character string being one or more of the phenotypes in \code{out_lmm}
#' @param drop include all SNPs within \code{drop} of max LOD (default 1.5)
#' @param min_lod minimum LOD retained (default 1.5)
#'
#' @return table of top_snps at maximum lod for \code{pattern}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords hplot
#'
#' @examples
#' dontrun(topsnp_pattern(out_lmm_snps, proband))
#'
#' @export
#' @importFrom dplyr distinct
topsnp_pattern <- function(out_lmm_snps, proband,
                           drop = 1.5, min_lod = 3) {

  ## Want to do pull top SNPs (and InDels and SVs) for combo of
  ##   proband  phenotype namess
  ##   sdp  strain distribution pattern
  ##   type     SNP, InDel, or SVS (DEL,INS,INV,...)
  ## Do not need to use top_snps().
  ## Restructure this routine.

  phename <- dimnames(out_lmm_snps$lod)[[2]]
  chr_id <- names(out_lmm_snps$map)
  if(length(chr_id) != 1)
    stop("only 1 chr allowed")

  #############################
  ## Fudge with scan1 output to get max over multiple proband traits.
  lodcol <- match(proband, phename)
  if(!length(lodcol))
    stop(paste(proband, "not part of scan1 object"))

  ## Following is simplified from qtl2scan::top_snps()
  lod <- out_lmm_snps$lod[, lodcol, drop = FALSE]
  n_mar <- nrow(lod)
  n_phe <- ncol(lod)

  ## String lods out in one column; replicate map and pheno names.
  pheno <- matrix(proband, n_mar, n_phe, byrow = TRUE)

  ## Keep only what is large enough.
  max_lod <- apply(lod, 2, max, na.rm = TRUE)
  keep <- apply(lod, 1, function(x, max_lod, min_lod, drop) {
    !is.na(x) & x > max_lod - drop & x >= min_lod
  }, max_lod, min_lod, drop)
  if(!is.matrix(keep)) { ## case of 1 proband
    keep <- t(as.matrix(keep))
    dimnames(keep)[[1]] <- proband
  }

  ## Get snpinfo and add lod and phenos
  snpinfo <- out_lmm_snps$snpinfo[[chr_id]]
  snpindex <- snpinfo$index
  snpi <- matrix(FALSE, nrow(snpinfo), n_phe,
                 dimnames = list(NULL, proband))
  for(i in proband) {
    snpi[,i] <- snpindex %in% which(keep[i,])
  }
  snpinfo <- snpinfo[unlist(apply(snpi, 2, which)),]
  snpinfo$lod <- lod[snpindex,][snpi]
  snpinfo$pheno <- pheno[snpindex,][snpi]

  #############################
  ## Get top SNPS for proband
  snpinfo <- dplyr::distinct(snpinfo, pheno, snp_id, pos_Mbp, .keep_all=TRUE)

  if(!nrow(snpinfo))
    return(NULL)
  class(snpinfo) <- c("topsnp_pattern", class(snpinfo))
  snpinfo
}
#' Summary of topsnp_pattern object
#'
#' @param object of class \code{topsnp_pattern}
#' @param type type of summary
#' @param ... other arguments not used
#'
#' @return table summary
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method summary topsnp_pattern
#' @rdname topsnp_pattern
#' @export
#' @importFrom dplyr arrange desc group_by mutate summarize ungroup
#' @importFrom CCSanger sdp_to_pattern
summary.topsnp_pattern <- function(object,
                                   type = c("best","common"),
                                   ...) {
  if(!nrow(object))
    return(NULL)

  type <- match.arg(type)
  switch(type,
         best = { ## Top SNPs across all phenotypes.
           if(!nrow(object))
             return(NULL)
           dplyr::arrange(
             dplyr::mutate(object,
                           pattern = CCSanger::sdp_to_pattern(sdp)),
             dplyr::desc(lod))
           },
         common = { ## Find most common patterns by pheno.
           dplyr::arrange(
             dplyr::mutate(
               dplyr::ungroup(
                 dplyr::summarize(
                   dplyr::group_by(object, pheno,sdp),
                   count = n(),
                   pct = round(100 * n() / nrow(object), 2),
                   min_lod = min(lod),
                   max_lod = max(lod),
                   max_snp = snp_id[which.max(lod)],
                   max_pos = pos_Mbp[which.max(lod)])),
               pattern = CCSanger::sdp_to_pattern(sdp)),
             dplyr::desc(max_lod))
         })
}
#' Plot of topsnp_pattern object
#'
#' @param x of class \code{topsnp_pattern}
#' @param ... other arguments not used
#' @param title title of plot (default "")
#' @param group group plot by one of c("pheno","pattern")
#'
#' @return ggplot2 object
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @method plot topsnp_pattern
#' @rdname topsnp_pattern
#' @export
#' @importFrom dplyr arrange desc filter group_by mutate summarize ungroup
#' @importFrom CCSanger sdp_to_pattern
#' @importFrom stringr str_replace
#' @importFrom ggplot2 aes facet_wrap geom_path geom_point ggplot ggtitle scale_shape_manual
plot.topsnp_pattern <- function(x, ...,
                                title="",
                                group=c("pheno","pattern")) {
  ## Lots of missing values.
  ## Plot has turned into spikes somehow.

  if(!nrow(x))
    return(plot_null())

  group <- match.arg(group)

  ## Create pattern from sdp
  x <- dplyr::mutate(x,
                     pattern = CCSanger::sdp_to_pattern(sdp))

  ## Collapse patterns with only 1.
  ones <- (dplyr::ungroup(
    dplyr::filter(
      dplyr::summarize(dplyr::group_by(x, sdp),
                       count = n()),
      count == 1)))$sdp
  if(length(ones) > 1) {
    x <- dplyr::mutate(x,
                       pattern = stringr::str_replace(pattern,
                                                      paste(ones, collapse="|"),
                                                      "rest"))
  }

  ## Type of difference among founders.
  if(is.null(x$type))
    x$type <- "SNP"
  else {
    x$type <- substring(x$type,1,3)
    x$type[x$type == "InD"] <- "indel"
  }
  x$type <- factor(x$type)

  ## See http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  shapes <- c(SNP=96,indel=23,INS=25,DEL=24,INV=22)
  ## Add diamond shape to any overlooked above.
  tmp <- levels(x$type) %in% names(shapes)
  if(any(!tmp)) {
    newshapes <- levels(x$type)[!tmp]
    shapes <- c(shapes, rep(21,length(newshapes)))
    names(shapes)[-(1:5)] <- newshapes
  }
  shapes <- shapes[levels(x$type)]

  ## Arrange patterns from largest to smallest LOD.
  x$pattern <- factor(x$pattern,
                       levels=(dplyr::ungroup(
                         dplyr::select(
                           dplyr::arrange(
                             dplyr::summarize(
                               dplyr::group_by(x, pattern),
                               max = max(lod)),
                             dplyr::desc(max)),
                           pattern)))[[1]])

  ## Plot LOD by position for common top SNP patterns.
  if(group == "pheno") {
    p <- ggplot2::ggplot(dplyr::arrange(x, pos_Mbp),
                         ggplot2::aes(x=pos_Mbp, y=lod, col=pattern)) +
      ggplot2::facet_wrap(~pheno)
  } else { # group == "pattern"
    p <- ggplot2::ggplot(dplyr::arrange(x, pos_Mbp),
                         ggplot2::aes(x=pos_Mbp, y=lod, col=pheno)) +
      ggplot2::facet_wrap(~ pattern)
  }
  p + ggplot2::geom_point(ggplot2::aes(shape = type),
                          size = 3, fill = "grey40") +
    ggplot2::scale_shape_manual(name = "SV Type",
                                labels = names(shapes),
                                values = shapes) +
    ggplot2::geom_path() +
    ggplot2::ggtitle(title)
}
