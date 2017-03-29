#' Allele plot for SNPs, alleles and allele pairs
#' 
#' Create table of alleles for various model fits.
#' 
#' @param scan_apr Object of class \code{\link[qtl2scan]{scan1}} with allele scan.
#' @param coefs Object of class \code{\link[qtl2scan]{scan1coef}}.
#' @param coefs36 Object of class \code{\link[qtl2scan]{scan1coef}} with allele pair coefficients.
#' @param scan_pat Object of class \code{\link{scan_pattern}}.
#' @param map Genome map.
#' @param haplo Haplotype allele letter to compare for number of copies.
#' @param trim If \code{TRUE}, trim extreme alleles.
#' 
#' @return Table with allele effects across sources.
#' 
#' @export
#' 
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows mutate
#' @importFrom stringr str_count str_detect str_split
#' 
allele1 <- function(scan_apr, coefs, coefs36, scan_pat, map, haplo,
                    trim = TRUE) {
  
  # Combine effects estimates.
  scan_pats <- dplyr::rename(
    tidyr::gather(summary(scan_pat, map),
                  allele, effect, -pheno, -chr, -pos, -lod),
    source = pheno)
  alleles <- dplyr::bind_rows(
    haplo = dplyr::mutate(
      tidyr::gather(
        summary(coefs, scan_apr, map), 
        allele, effect, -pheno, -chr, -pos, -lod),
      pheno = as.character(pheno),
      chr = as.character(chr)),
    diplo = dplyr::mutate(
      tidyr::gather(
        summary(coefs36, scan_apr, map), 
        allele, effect, -pheno, -chr, -pos, -lod),
      pheno = as.character(pheno),
      chr = as.character(chr)),
    .id = "source")
  
  alleles <- dplyr::bind_rows(alleles, scan_pats)
  alleles$source <- factor(alleles$source, c("haplo","diplo",
                                             colnames(scan_pat$scan)))
  if(trim)
    alleles <- trim_quant(alleles)
  
  tmpfn <- function(x) 
    sapply(stringr::str_split(x, ":"), 
           function(x) stringr::str_detect(x[2], haplo))
  alleles <- dplyr::mutate(alleles,
                           probe = stringr::str_count(allele, haplo),
                           probe = ifelse(tmpfn(source) & (allele == "het"), 1, probe),
                           probe = ifelse(tmpfn(source) & (allele == "alt"), 2, probe),
                           probe = factor(probe))
  attr(alleles, "probe") <- haplo
  class(alleles) <- c("allele1", class(alleles))
  alleles
}
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' 
summary.allele1 <- function(object, ...) {
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, source),
      pos = pos[1],
      min = min(effect),
      max = max(effect)))
}
trim_quant <- function(object, beyond = 3) {
  quant <- quantile(object$effect, c(.25,.75))
  range <- quant + c(-1,1) * beyond * diff(quant)
  object$effect <- pmin(pmax(object$effect, range[1]), range[2])
  object
}
trim_ends <- function(object, ends = 1, trim_val = NULL, effname = "effect") {
  ends <- rep(ends, length = 2)
  rk <- rank(object[[effname]])
  mintrim <- (rk <= ends[1])
  maxtrim <- (rk >= 1 + nrow(object) - ends[2])
  if(is.null(trim_val)) {
    trim_val <- c(object[[effname]][which(rk == min(rk[!mintrim]))[1]],
                  object[[effname]][which(rk == max(rk[!maxtrim]))[1]])
  }
  object[[effname]][mintrim] <- trim_val[1]
  object[[effname]][maxtrim] <- trim_val[2]
  object
}
