#' Allele plot for SNPs, alleles and allele pairs
#' 
#' Create table of alleles for various model fits.
#' 
#' @param phe_df data frame with one phenotype
#' @param cov_mx covariate matrix
#' @param probD object of class \code{\link[qtl2geno]{calc_genoprob}}
#' @param map list of genome maps
#' @param K_chr kinship matrix
#' @param patterns data frame of pattern information
#' @param alt Haplotype allele letter(s) for alternative to reference.
#' @param trim If \code{TRUE}, trim extreme alleles.
#' @param ... additional parameters
#' 
#' @param scanH Object of class \code{\link[qtl2scan]{scan1}} with allele scan.
#' @param coefH Object of class \code{\link[qtl2scan]{scan1coef}}.
#' @param coefD Object of class \code{\link[qtl2scan]{scan1coef}} with allele pair coefficients.
#' @param scan_pat Object of class \code{\link{scan_pattern}}.
#' 
#' @return Table with allele effects across sources.
#' 
#' @export
#' 
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows filter mutate
#' @importFrom stringr str_count str_detect str_split
#' 
allele1 <- function(phe_df=NULL, cov_mx=NULL, probD=NULL, map=NULL, K_chr=NULL, patterns=NULL, 
                    alt=NULL, trim = TRUE, ...) {
  allele1_internal(phe_df, cov_mx, probD, map, K_chr, patterns,
                   alt, trim = TRUE, ...)
}
allele1_internal <- function(
  phe_df, cov_mx, probD, map, K_chr, patterns, alt, trim = TRUE, 
  probH = qtl2geno::genoprob_to_alleleprob(probD),
  scanH = qtl2scan::scan1(probH, phe_df, K_chr, cov_mx),
  coefH = qtl2scan::scan1coef(probH, phe_df, K_chr, cov_mx),
  coefD = qtl2scan::scan1coef(probD, phe_df, K_chr, cov_mx),
  scan_pat = NULL,
  ...) 
{
  if(!is.null(phe_df) && ncol(phe_df) > 1) {
    warning("only using first phenotype")
    phe_df <- phe_df[, 1, drop = FALSE]
  }
  pheno_name <- names(as.data.frame(phe_df))
  
  # Get patterns for pheno.
  if(!is.null(scan_pat)) {
    patterns <- dplyr::rename(scan_pat$patterns,
                              pattern = founders)
  } else {
    if(is.null(patterns))
      stop("need either patterns or scan_pat")
    patterns <- dplyr::filter(patterns,
                              pheno == pheno_name)
    scan_pat <- qtl2pattern::scan_pattern(probD, phe_df, K_chr, cov_mx,
                                          map, patterns)
  }
  if(!nrow(patterns))
    return(NULL)

  if(is.null(alt))
    alt <- paste0(
      stringr::str_split(
        stringr::str_replace(patterns$pattern[1], ".*:", ""),
        "")[[1]],
      collapse = "|")

  # Combine effects estimates.
  mar_df <- function(x, n) {
    mar <- rownames(x)
    x <- as.data.frame(x[, seq_len(n)])
    x$mar <- mar
    x
  }

  scan_pats <- tidyr::gather(
    dplyr::bind_rows(
      lapply(scan_pat$coef, mar_df, 3), 
      .id = "source"),
    allele, effect, -mar, -source)
  
  alleles <- dplyr::bind_rows(
    haplo = tidyr::gather(
      mar_df(coefH, 8),
      allele, effect, -mar),
    diplo = tidyr::gather(
        mar_df(coefD, 36), 
        allele, effect, -mar),
    .id = "source")
  
  alleles <- dplyr::bind_rows(alleles, scan_pats)
  map <- map[[1]]
  mar <- names(map)
  map <- data.frame(pos=map, mar=mar, stringsAsFactors = FALSE)
  alleles <- dplyr::inner_join(alleles, map, by = "mar")
  alleles$source <- factor(alleles$source, c("haplo","diplo",
                                             names(scan_pat$coef)))
  if(trim)
    alleles <- trim_quant(alleles)
  
  tmpfn <- function(x) 
    sapply(stringr::str_split(x, ":"), 
           function(x) stringr::str_detect(x[2], alt))
  alleles <- dplyr::mutate(alleles,
                           probe = stringr::str_count(allele, alt),
                           probe = ifelse(tmpfn(source) & (allele == "het"), 1, probe),
                           probe = ifelse(tmpfn(source) & (allele == "alt"), 2, probe),
                           probe = factor(probe))
  attr(alleles, "probe") <- alt
  class(alleles) <- c("allele1", class(alleles))
  alleles
}
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' 
summary.allele1 <- function(object, scan1_object=NULL, map=NULL, pos=NULL, ...) {
  if(is.null(pos)) {
    if(is.null(scan1_object))
      pos_Mbp <- median(object$pos)
    else
      pos_Mbp <- summary(scan1_object, map)$pos[1]
  } else {
    pos_Mbp <- pos
    if(pos_Mbp < min(object$pos) | pos_Mbp > max(object$pos))
      stop("position must be within range of allele positions")
  }
  
  tmpfn <- function(pos) {
    a <- abs(pos - pos_Mbp)
    which(a == min(a))
  }
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, source),
      min = min(effect[tmpfn(pos)]),
      max = max(effect[tmpfn(pos)]),
      pos = pos_Mbp))
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
