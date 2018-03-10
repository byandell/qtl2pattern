#' Pull mediator and SDP from mediate object
#' 
#' @param mediate_obj object from \code{\link{mediate1_test}}
#' @param med_ls object from \code{link{pheno_region}}
#' @param patterns data frame \code{summary} of \code{\link{tob_snps_all}} object
#' @param med_name name of mediator to pull
#' 
#' @export
#' @importFrom dplyr filter
#' 
pull_mediator <- function(mediate_obj, med_ls, patterns,
                         med_name = med_names[1]) {
  triad <- levels(mediate_obj$best$triad)[1]
  medID <- ifelse("symbol" %in% names(mediate_obj$best), "symbol", "longname")
  med_names <- dplyr::filter(mediate_obj$best, triad == triad)[[medID]]
  
  id <- med_ls[[2]]$id[med_ls[[2]][[medID]] == med_name]
  
  sdps <- unique(dplyr::filter(patterns, pheno == pheno_name)$sdp)
  pattern <- qtl2pattern::sdp_to_pattern(sdps, LETTERS[1:8])
  
  list(sdp = sdps[qtl2pattern::sdp_to_pattern(sdps, LETTERS[1:8]) == pattern][1],
       mediator = med_ls[[1]][, id, drop = FALSE])
}
