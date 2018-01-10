snpinfo_to_haplos <- function(snpinfo) {
  alleles <- dplyr::select(
    snpinfo,
    -(snp_id:alleles))
  # Would be better to have object that gives allele names rather than this opposite approach.
  infonames <- c("consequence","type","sdp","index","interval","on_map","pheno","lod","ensembl_gene")
  if(length(wh <- which(infonames %in% names(alleles)))) {
    alleles <- alleles[, -wh, drop = FALSE]
  }
  # Columns in between consequence and type should be alleles.
  # If not provided, assume we are in mouse with 8.
  if((nc <- ncol(alleles)) < 2) {
    warning("no alleles in snpinfo; assuming 8")
    nc <- 8
  }
  LETTERS[seq_len(nc)]
}
