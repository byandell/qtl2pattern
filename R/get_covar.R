#' Get covariates
#'
#' Get covariates using \code{analyses_tbl} filtered by \code{phename_output}.
#'
#' @param covar matrix of all covariates
#' @param analyses_tbl table of analyses setups
#'
#' @return matrix of covariates
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_covar(covar, datapath)}
#'
#' @export
get_covar <- function(covar, analyses_tbl) {
  ## Get covariate matrix covar.
  logical.col <- sapply(analyses_tbl, function(x) any(as.logical(x)))
  covar[, match(names(analyses_tbl)[logical.col],
                dimnames(covar)[[2]], nomatch=0),
        drop=FALSE]
}
