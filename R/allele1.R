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
#' 
#' @return Table with allele effects across sources.
#' 
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows mutate
#' @importFrom stringr str_count str_detect str_split
allele1 <- function(scan_apr, coefs, coefs36, scan_pat, map, haplo) {
  
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
    diplo_EF = dplyr::mutate(
      tidyr::gather(
        summary(coefs36, scan_pat$scan, map), 
        allele, effect, -pheno, -chr, -pos, -lod),
      pheno = as.character(pheno),
      chr = as.character(chr)),
    diplo_E = dplyr::mutate(
      tidyr::gather(
        summary(coefs36, subset(scan_pat$scan, lodcolumn=2), DOex$pmap), 
        allele, effect, -pheno, -chr, -pos, -lod),
      pheno = as.character(pheno),
      chr = as.character(chr)),
    .id = "source")
  
  alleles <- dplyr::bind_rows(alleles, scan_pats)
  
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
summary.allele1 <- function(object, ..., trim = TRUE) {
  if(trim)
    object <- trim_quant(object)
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, source),
      pos = pos[1],
      min = min(effect),
      max = max(effect)))
}
