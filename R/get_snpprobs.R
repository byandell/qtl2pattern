#' Get SNP genotype probabilities in window around peak
#'
#' Get SNP information from SQLite database within \code{window_Mbp} of \code{peak_Mbp} on \code{chri_id}
#'
#' @param chr_id chromosome identifier
#' @param peak_Mbp position in Mbp of peak
#' @param window_Mbp half-width of \code{window} around \code{peak_Mbp}
#' @param phename names of phenotypes
#' @param probs_obj object of class \code{\link[qtl2]{calc_genoprob}} for \code{chr_id}
#' @param probs_map map of markers/pseudomarkers in \code{probs_obj}
#' @param snpinfo SNP information table from user supplied \code{query_variants}; see \code{\link[qtl2]{create_variant_query_func}}
#'
#' @return list with \code{snpprobs} and \code{snpinfo}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \donttest{get_snpprobs(chr_id, peak_Mbp, window_Mbp, scan_obj, probs_obj, datapath)}
#'
#' @export
#'
#' @importFrom dplyr bind_rows mutate
#' @importFrom qtl2 genoprob_to_snpprob index_snps
#'
get_snpprobs <- function(chr_id=NULL, peak_Mbp=NULL, window_Mbp=NULL,
                         phename, probs_obj, probs_map,
                         snpinfo = query_variants(chr_id,
                                                  peak_Mbp - window_Mbp,
                                                  peak_Mbp + window_Mbp)) {
  if(is.null(chr_id) | is.null(peak_Mbp) | is.null(window_Mbp))
    return(NULL)

  if(window_Mbp == 0) {
    window_Mbp <- 3
  }
  if(peak_Mbp == 0) {
    peak_Mbp <- mean(range(probs_map[[1]]))
  }

  snpinfo <- qtl2::index_snps(probs_map, snpinfo)
  
  list(snpprobs = qtl2::genoprob_to_snpprob(probs_obj, snpinfo),
       snpinfo = snpinfo)
}
