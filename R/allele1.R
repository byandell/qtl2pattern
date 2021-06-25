#' Allele plot for SNPs, alleles and allele pairs
#' 
#' Create table of alleles for various model fits.
#' 
#' @param probD object of class \code{\link[qtl2]{calc_genoprob}}
#' @param phe_df data frame with one phenotype
#' @param cov_mx covariate matrix
#' @param map list of genome maps
#' @param K_chr kinship matrix
#' @param patterns data frame of pattern information
#' @param alt Haplotype allele letter(s) for alternative to reference.
#' @param blups Create BLUPs if \code{TRUE}
#' @param ... additional parameters
#' 
#' @return Table with allele effects across sources.
#' 
#' @export
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr bind_rows filter mutate select
#' @importFrom stringr str_count str_detect str_split
#' @importFrom qtl2 genoprob_to_alleleprob scan1 scan1blup scan1coef
#' @importFrom rlang .data
#' 
allele1 <- function(probD, phe_df=NULL, cov_mx=NULL, map=NULL, K_chr=NULL, patterns=NULL, 
                    alt=NULL, blups = FALSE, ...) {
  allele1_internal(phe_df, cov_mx, probD, map, K_chr, patterns,
                   alt, blups, ...)
}
allele1_internal <- function(
  phe_df, cov_mx, probD, map, K_chr, patterns, alt, blups,
  probH = qtl2::genoprob_to_alleleprob(probD),
  scanH = qtl2::scan1(probH, phe_df, K_chr, cov_mx),
  coefH = scan1fn(probH, phe_df, K_chr, cov_mx),
  coefD = qtl2::scan1coef(probD, phe_df, K_chr, cov_mx),
  scan_pat = NULL,
  drop_covar = TRUE,
  ...) 
{
  if(!is.null(phe_df) && ncol(phe_df) > 1) {
    warning("only using first phenotype")
    phe_df <- phe_df[, 1, drop = FALSE]
  }
  pheno_name <- names(as.data.frame(phe_df))
  
  scan1fn <- ifelse(blups, 
                    qtl2::scan1blup, 
                    qtl2::scan1coef)
  
  # Get patterns for pheno.
  if(!is.null(scan_pat)) {
    patterns <- dplyr::rename(scan_pat$patterns,
                              pattern = .data$founders)
  } else {
    if(is.null(patterns))
      stop("need either patterns or scan_pat")
    patterns <- dplyr::filter(patterns,
                              .data$pheno == pheno_name)
    scan_pat <- scan1pattern(probD, phe_df, K_chr, cov_mx,
                             map, patterns, blups = blups)
  }
  if(!nrow(patterns))
    return(NULL)

  if(is.null(alt))
    alt <- paste0(
      stringr::str_split(
        stringr::str_replace(patterns$pattern[1], ".*:", ""),
        "")[[1]],
      collapse = "|")

  # Combine effects estimates. Object a has information about number of alleles.
  mar_df <- function(x, a) {
    mar <- rownames(x)
    n <- dim(a[[1]])[2]
    x <- as.data.frame(x[, seq_len(n)])
    x$mar <- mar
    x
  }

  scan_pats <- tidyr::pivot_longer(
    dplyr::select(
      dplyr::bind_rows(
        lapply(scan_pat$coef, mar_df, scan_pat$coef), 
        .id = "source"),
      -.data$intercept),
    .data$ref:.data$alt, names_to = "allele", values_to = "effect")
  
  # Identify covariates
  scan_pats <- dplyr::mutate(
    scan_pats,
    covar = !(.data$allele %in% c("alt","ref","het"))
  )
  if(drop_covar) {
    scan_pats <- dplyr::filter(
      scan_pats,
      !.data$covar)
  }
  
  # Combine alleles and allele pairs using codes inferred from probH and probD.
  alleles <- dplyr::bind_rows(
    haplotype = tidyr::pivot_longer(
      mar_df(coefH, probH),
      -.data$mar, names_to = "allele", values_to = "effect"),
    diplotype = tidyr::pivot_longer(
      mar_df(coefD, probD), 
      -.data$mar, names_to = "allele", values_to = "effect"),
    .id = "source")
  # Bind alleles with scan patterns (ref/het/alt)
  alleles <- dplyr::bind_rows(alleles, scan_pats)
  # Make sure covar has no missing values for alleles.
  alleles <- dplyr::mutate(
    alleles,
    covar = ifelse(is.na(.data$covar), FALSE, .data$covar))
  
  map <- map[[1]]
  mar <- names(map)
  map <- data.frame(pos=map, mar=mar, stringsAsFactors = FALSE)
  alleles <- dplyr::inner_join(alleles, map, by = "mar")
  alleles$source <- factor(alleles$source, c("haplotype","diplotype",
                                             names(scan_pat$coef)))
  
  # Set up probe for color in plot as number of copies of alternative allele.
  # ref = 0, het = 1, alt = 2; covariates = 5.
  tmpfn <- function(x) stringr::str_detect(x, ":")
  alleles <- dplyr::mutate(alleles,
                           probe = stringr::str_count(.data$allele, alt),
                           probe = ifelse(tmpfn(.data$source) & (.data$allele == "het"), 1, .data$probe),
                           probe = ifelse(tmpfn(.data$source) & (.data$allele == "alt"), 2, .data$probe),
                           probe = ifelse(.data$covar, 5, .data$probe),
                           probe = factor(.data$probe))
  attr(alleles, "probe") <- alt
  attr(alleles, "blups") <- blups
  class(alleles) <- c("allele1", class(alleles))
  alleles
}
#' @export
#' @importFrom dplyr filter
#' 
subset.allele1 <- function(x, sources = levels(x$source), ...) {
  dplyr::filter(x, source %in% sources)
}
#' @export
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom qtl2 find_peaks
#' 
summary.allele1 <- function(object, scan1_object=NULL, map=NULL, pos=NULL, ...) {
  if(is.null(pos)) {
    if(is.null(scan1_object))
      pos_center <- median(object$pos)
    else
      pos_center <- qtl2::find_peaks(scan1_object, map)$pos[1]
  } else {
    pos_center <- pos
    if(pos_center < min(object$pos) | pos_center > max(object$pos))
      stop("position must be within range of allele positions")
  }
  
  tmpfn <- function(pos, pos_center) {
    a <- abs(pos - pos_center)
    which(a == min(a))
  }
  dplyr::ungroup(
    dplyr::summarize(
      dplyr::group_by(object, .data$source),
      min = min(.data$effect[tmpfn(.data$pos, pos_center)]),
      lo_25 = quantile(.data$effect[tmpfn(.data$pos, pos_center)], 0.25),
      median = median(.data$effect[tmpfn(.data$pos, pos_center)]),
      hi_75 = quantile(.data$effect[tmpfn(.data$pos, pos_center)], 0.75),
      max = max(.data$effect[tmpfn(.data$pos, pos_center)]),
      pos = pos_center))
}
