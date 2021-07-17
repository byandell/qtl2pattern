#' Convert SNP info to map
#'  
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived from with
#'     \code{\link[qtl2]{index_snps}}):
#' \itemize{
#' \item \code{chr} - Character string or factor with chromosome
#' \item \code{pos} - Position (in same units as in the \code{"map"}
#'     attribute in \code{genoprobs}.
#' \item \code{sdp} - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' \item \code{snp} - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' \item \code{index} - Indices that indicate equivalent
#'     groups of SNPs.
#' \item \code{intervals} - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' \item \code{on_map} - Indicate whether SNP coincides with a marker
#'     in the \code{genoprobs}
#' }
#'
#' @return map as list of vectors of marker positions.
#' 
snpinfo_to_map <-
  function(snpinfo)
  {
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
      stop("snpinfo$index values outside of range [1, ",
           nrow(snpinfo), "]")
    
    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)
    
    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq(along=map)) {
      u <- unique(index[[i]])
      map[[i]] <- map[[i]][u]
      names(map[[i]]) <- snp[[i]][u]
    }
    
    names(map) <- uchr
    
    map
  }
