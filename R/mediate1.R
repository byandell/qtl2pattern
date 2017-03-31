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
#'
#' @export
#' @importFrom qtl2geno find_marker
#' @importFrom qtl2scan scan1
#' @importFrom DOread read_mrna
#'
mediate1 <- function(chr_id, pos_Mbp, window_Mbp, 
                     phe_df, cov_mx=NULL, probs_obj, kinship=NULL, map,
                     datapath, verbose = FALSE) {

  if(ncol(phe_df) > 1) {
    message("only using first phenotype")
    phe_df <- as.data.frame(phe_df[, 1, drop = FALSE])
  }

  # Find marker at pos_Mbp on chromosome chr_id
  mar_id <- qtl2geno::find_marker(map, chr_id, pos_Mbp)
  names(pos_Mbp) <- mar_id

  # Subset genotype probabilities to this marker.
  probs_max <- subset(probs_obj, chr = chr_id, mar = mar_id)

  # Raw fit
  scan_max <- qtl2scan::scan1(probs_max, phe_df, kinship[[chr_id]],
                              cov_mx)
  raw_lod <- max(scan_max, map)[, 3]

  scan_window <- pos_Mbp + c(-1,1) * window_Mbp
  ## RNAseq data
  # Get subset of expression data within scan window.
  indID <- rownames(phe_df)
  expr.mrna <- DOread::read_mrna(indID, chr_id, scan_window[1], scan_window[2], datapath)
  annot.mrna <- expr.mrna$annot
  expr.mrna <- expr.mrna$expr

  annot.mrna$lod_t_m <- rep(NA, nrow(annot.mrna))
  for(i in seq_len(nrow(annot.mrna))) {
    if(verbose) cat(i, "")
    scan_max <- qtl2scan::scan1(probs_max, phe_df, kinship[[chr_id]],
                                cbind(cov_mx,
                                      expr.mrna[, annot.mrna$id[i], drop = FALSE]))
    annot.mrna$lod_t_m[i] <- max(scan_max, pos_Mbp)[, 3]
  }

  attr(annot.mrna, "pos") <- pos_Mbp
  attr(annot.mrna, "lod") <- raw_lod
  attr(annot.mrna, "pheno") <- names(phe_df)

  class(annot.mrna) <- c("mediate1", class(annot.mrna))
  annot.mrna
}
