#' Collapse genoprob according to pattern
#'
#' @param snpprobs object of class \code{\link[qtl2]{calc_genoprob}}
#' @param action SNP gene action type
#'
#' @return object of class \code{\link[qtl2]{calc_genoprob}}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' \donttest{
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Download SNP info for DOex from web and read as RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Convert to snp probabilities
#' snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
#' snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
#' 
#' 
#' dim(snppr[[1]])
#' dim(snpprob_collapse(snppr, "additive")[[1]])
#' }
#' 
#' @export
snpprob_collapse <- function(snpprobs,
                             action = c("additive","add+dom","non-add",
                                        "recessive","dominant","basic")) {
  if(dim(snpprobs[[1]])[2] == 2)
    return(snpprobs)
  if(is.null(action))
    return(snpprobs)
  if(action %in% c("add+dom","basic")) # basic is null action
    return(snpprobs)
  if(action == "additive") {
    patprobs <- snpprobs[[1]][,c(1,3),,drop=FALSE]
    dimnames(patprobs)[[2]] <- as.character(0:1)
    patprobs <- patprobs + 0.5 * snpprobs[[1]][,c(2,2),,drop=FALSE]
  } else {
    snpsets <- switch(action,
                      "non-add" = c(0,1,0),
                      "recessive" = c(0,1,1),
                      "dominant" = c(0,0,1))
    
    ## Assume snpsets only has 0s and 1s
    if(all(snpsets==1) | all(snpsets==0) | length(snpsets) != 3)
      stop("need mix of 0s and 1s")
    
    snpsets <- as.character(snpsets)
    patprobs <- snpprobs[[1]][, !duplicated(snpsets),, drop = FALSE]
    dimnames(patprobs)[[2]] <- unique(snpsets)
    patprobsa <- snpsets[duplicated(snpsets)]
    patprobs[,patprobsa,] <- patprobs[, patprobsa,, drop = FALSE] + 
      snpprobs[[1]][, duplicated(snpsets),, drop = FALSE]
  }
  modify_aprobs(snpprobs, patprobs)
}
