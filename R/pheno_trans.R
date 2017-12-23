#' Get phenotypes
#'
#' Get phenotypes using data frame of phenopypes filtered by \code{analyses_tbl}
#'
#' @param phe phenotypes in data frame
#' @param phename vector of phenotype names (subset of \code{colnames(phe)})
#' @param transform vector of function names (\code{NULL} for no transformations)
#' @param offset vector of offsets
#' @param winsor vector of winsorize values
#'
#' @return data frame of phenotypes
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{get_pheno(phe, analyses_tbl)}
#'
#' @export
#' @importFrom assertthat assert_that is.number
#'
pheno_trans <- function(phe, phename, transform = NULL, offset = 0,
                        winsor = 0.02) {
  # Get phenotype names. Make sure it is character, not factor.
  tmp <- !duplicated(phename)
  assertthat::assert_that(all(tmp))
  phename <- as.character(phename)
  phe <- phe[, match(phename, colnames(phe), nomatch=0), drop=FALSE]

  if(!is.null(transform)) {
    if(length(transform) == 1)
      transform <- rep_len(transform, length(phename))
    assertthat::assert_that(length(phename) == length(transform))
    
    ## Transform phenotype.
    not.id <- !(transform %in% c("id","identity"))
    if(any(not.id)) {
      if(length(offset == 1))
        offset <- rep_len(offset, length(phename))
      assertthat::assert_that(length(phename) == length(offset))
      for(i in which(not.id)) {
        tmp <- phe[, phename[i]] + offset[i]
        if(transform[i] %in% c("log","sqrt")) {
          # This will still yield -Inf if transform is log
          tmp[tmp < 0] <- NA
        }
        phe[,phename[i]] <- get(transform[i])(tmp)
      }
    }
    
    ## Parameter to winsorize will later be a value.
    if(any(is.logical(winsor)))
      winsor <- ifelse(winsor, 0.02, 0)
    if(length(winsor == 1))
      winsor <- rep_len(winsor, length(phename))
    assertthat::assert_that(length(phename) == length(winsor))
    wh <- which(winsor > 0)
    if(length(wh)) {
      for(i in wh)
        phe[,phename[i]] <- winsorize(unlist(phe[,phename[i]]), winsor[i])
    }
  }
  phe
}
