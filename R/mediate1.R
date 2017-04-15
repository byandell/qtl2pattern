#' Mediate phenotype with nearby expression phenotypes
#'
#' @param chr_id Chromosome identifier.
#' @param pos_Mbp Position in megabase pairs.
#' @param window_Mbp Window width to scan for expression traits.
#' @param phe_df Data frame with phenotype.
#' @param cov_mx Covariate matrix
#' @param probs_obj Object from \code{\link[qtl2feather]{feather_genoprob}}
#' @param kinship Kinship matrix list.
#' @param map Physical map.
#' @param datapath Directory name for expression traits objects.
#' @param verbose Show expression trait number as calculated if \code{TRUE}
#' @param region Mediate over whole region if \code{TRUE}
#'
#' @export
#' @importFrom qtl2geno find_marker
#' @importFrom qtl2scan get_common_ids scan1
#' @importFrom DOread read_mrna
#'
mediate1 <- function(chr_id, pos_Mbp, window_Mbp, 
                     phe_df, cov_mx=NULL, probs_obj, kinship=NULL, map,
                     datapath, verbose = FALSE,
                     region = FALSE) {

  if(ncol(phe_df) > 1) {
    message("only using first phenotype")
    phe_df <- as.data.frame(phe_df[, 1, drop = FALSE])
  }

  scan_window <- pos_Mbp + c(-1,1) * window_Mbp
  
  ## RNAseq data
  # Get subset of expression data within scan window.
  indID <- rownames(phe_df)
  expr_mrna <- DOread::read_mrna(indID, chr_id, scan_window[1], scan_window[2], datapath)
  annot_mrna <- expr_mrna$annot
  expr_mrna <- expr_mrna$expr
  
  # Reduce to common data.
  if(!is.null(cov_tar)) {
    cov_mx <- as.matrix(cov_mx)
    if(is.null(colnames(cov_mx)))
      colnames(cov_mx) <- paste0("cov", seq_len(ncol(cov_mx)))
  }
  kinship <- kinship[[chr_id]]

  # Keep individuals with full records.
  ind2keep <- 
    qtl2scan::get_common_ids(phe_df, expr_mrna, cov_mx, kinship,
                             complete.cases = TRUE)
  phe_df <- phe_df[ind2keep,, drop = FALSE]
  expr_mrna <- expr_mrna[ind2keep,, drop = FALSE]
  if(!is.null(cov_mx))
    cov_mx <- cov_mx[ind2keep,, drop = FALSE]
  kinship <- kinship[ind2keep, ind2keep]

  # Decompose kinship
#  Ke <- qtl2scan::decomp_kinship(kinship)
  
  if(!region) {
    # Find marker at pos_Mbp on chromosome chr_id
    mar_id <- qtl2geno::find_marker(map, chr_id, pos_Mbp)

    # Subset genotype probabilities to this marker.
    probs_obj <- subset(probs_obj, chr = chr_id, mar = mar_id)
  }

  if(is.list(map))
    map <- map[[chr_id]]
  
  # Raw fit
  scan_max <- qtl2scan::scan1(probs_obj, phe_df, kinship, 
                              cov_mx)
  
  raw_lod <- max(scan_obj, map)
  raw_pos <- raw_lod$pos[1]
  raw_lod <- raw_lod[1, 3]
  
  n_mrna <- ncol(expr_mrna)
  scans <- matrix(NA, nrow(scan_obj), n_mrna)
  dimnames(scans) <- list(rownames(scan_obj), colnames(expr_mrna))
  scans[,1] <- scan_obj
  for(i in seq_len(n_mrna)) {
    if(verbose) cat(i, "")
    scans[,i] <- qtl2scan::scan1(probs_obj, phe_df, kinship,
                                 cbind(cov_mx,
                                       expr_mrna[, annot_mrna$id[i], drop = FALSE]))
  }
  scan_obj <- modify_object(scan_obj, scans)
  
  result <- list(scan = scan_obj, 
                 annot = annot_mrna,
                 map = map)
  attr(result, "pos") <- raw_pos
  attr(result, "lod") <- raw_lod
  attr(result, "pheno") <- names(phe_df)
  
  class(result) <- c("mediate1", class(result))
  result
}
