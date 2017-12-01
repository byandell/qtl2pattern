#' List of scan1coef objects
#'
#' Create a list of scan1coef objects using \code{\link[qtl2]{scan1coef}}.
#'
#' @param phe data frame of phenotypes
#' @param probs genotype probabilities object for one chromosome from \code{\link[qtl2]{calc_genoprob}}
#' @param K list of length 1 with kinship matrix
#' @param covar matrix of covariates
#' @param blups Create BLUPs if \code{TRUE}
#'
#' @return object of class \code{listof_scan1coeff}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{listof_scan1coef(probs, phe)}
#'
#' @export
#' @importFrom qtl2 scan1coef scan1blup
#' 
listof_scan1coef <- function(probs, phe, K=NULL, covar=NULL, blups = FALSE) {
  eff <- list()
  phename <- dimnames(phe)[[2]]
  scan1fn <- ifelse(blups, 
                    qtl2::scan1blup, 
                    qtl2::scan1coef)
  for(pheno in phename)
    eff[[pheno]] <- scan1fn(probs, phe[, pheno, drop=FALSE], K, covar)
  class(eff) <- c("listof_scan1coef", class(eff))
  eff
}

#' Summary of object of class listof_scan1coeff
#'
#' Summary of object of class \code{\link{listof_scan1coeff}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coeff}
#' @param scan1_object object from \code{scan1}
#'
#' @param map A list of vectors of marker positions, as produced by
#' \code{\link[qtl2]{insert_pseudomarkers}}.
#'
#' @param coef_names names of effect coefficients (default is all coefficient names)
#' @param center center coefficients if \code{TRUE}
#' @param ... arguments for \code{\link[qtl2ggplot]{ggplot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
#' @importFrom dplyr bind_cols
summary_listof_scan1coef <-
  function(object, scan1_object, map,
           coef_names = dimnames(object[[1]])[[2]],
           center = TRUE,
           ...) {
  phename <- names(object)
  chr_id <- names(map)
  
  # Summary of phenos on chr; reorder and rename by phename
  sum_chr <- summary(scan1_object, map, lodcolumn=phename, chr=chr_id)
  m <- match(make.names(phename), sum_chr$pheno)
  sum_chr <- sum_chr[m,]
  sum_chr$pheno <- phename
  
  pos <- sum_chr$pos
  names(pos) <- phename
  wh <- apply(scan1_object,2,function(x) which.max(x)[1])
  sum_coef <- as.data.frame(matrix(NA,
                                   nrow(sum_chr),
                                   length(coef_names),
                                   dimnames=list(NULL,coef_names)))
  for(pheno_id in phename) {
    if(!is.null(object[[pheno_id]])) {
      ## Need to save these in data frame.
      tmp <- object[[pheno_id]][wh[pheno_id], seq_along(coef_names)]
      if(center)
        tmp <- tmp - mean(tmp)
      sum_coef[match(pheno_id, sum_chr$pheno),] <- tmp
    }
  }
  dplyr::bind_cols(sum_chr,sum_coef)
}
#' @method summary listof_scan1coef
#' @rdname listof_scan1coef
#' @export
summary.listof_scan1coef <- function(object, ...)
  summary_listof_scan1coef(object, ...)

#' Summary of object of class listof_scan1coeff
#'
#' Summary of object of class \code{\link{listof_scan1coeff}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param object object of class \code{listof_scan1coeff}
#' @param scan1_object object from \code{scan1}
#' @param ... arguments for \code{\link[qtl2ggplot]{ggplot_coef}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
summary_scan1coef <-
  function(object, scan1_object, map, ...) {
    if(!inherits(object, "listof_scan1coef")) {
      object <- list(object)
      names(object) <- dimnames(scan1_object)[[2]][1]
    }
    summary_listof_scan1coef(object, scan1_object, map, ...)
  }

#' @method summary scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export summary.scan1coef
#'
summary.scan1coef <- function(object, ...)
  summary_scan1coef(object, ...)

#' Subset of object of class listof_scan1coeff
#'
#' Subset of object of class \code{\link{listof_scan1coeff}}, which is a list of objects of class \code{scan1coef}.
#'
#' @param x object of class \code{listof_scan1coeff}
#' @param elements indexes or names of list elements in x
#' @param ... ignored
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @rdname listof_scan1coef
#' @export
subset_listof_scan1coef <- function(x, elements, ...) {
  x_class <- class(x)
  class(x) <- "list"
  x <- x[elements]
  class(x) <- x_class
  x
}

#' @method subset listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export subset.listof_scan1coef
#'
subset.listof_scan1coef <- function(x, ...)
  subset_listof_scan1coef(x, ...)

#' @method [ listof_scan1coef
#' @rdname listof_scan1coef
#' @export
#' @export [.listof_scan1coef
#'
`[.listof_scan1coef` <- function(x, ...)
  subset_listof_scan1coef(x, ...)
