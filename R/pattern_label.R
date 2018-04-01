#' Turn genotype probabilities into labels
#' 
#' @param genomat matrix of genotype probabilities at locus
#' @param sdp SNP distribution pattern for plot colors
#' @param allele Driver has alleles if \code{TRUE}, otherwise allele pairs.
#' 
#' @export
pattern_label <- function(genos, sdp = NULL, allele = TRUE) {
  geno_names <- colnames(genos)
  if(allele)
    haplos <- geno_names
  else
    haplos <- unique(unlist(stringr::str_split(geno_names, "")))
  if(is.null(sdp)) {
    label <- apply(genos, 1,
                   # Pick geno name with highest value; if >1, paste together.
                   function(x) paste(geno_names[which.max(x)], collapse = ""))
    color <- label
  } else {
    alt <- haplos[sdp_to_logical(sdp, haplos)]
    # Set up allele pair labels and count as color.
    if(allele) {
      tmpfn <- function(x, geno_names) {
        # Mostly genos are close to 0, 0.5, 1 and add to 1; so 2 char names.
        # order by decreasing value
        o <- order(-x)
        # find how many needed to exceed .67
        i <- 1 + sum(cumsum(x[o]) <= .67)
        paste(rep(geno_names[sort(o[seq_len(i)])], length = 2), collapse = "")
      }
      label <- apply(genos, 1, tmpfn, geno_names)
      # Color is 0, 1, 2 as number of copies of `alt` letters.
      color <- factor(round(2 * apply(genos[, alt, drop = FALSE], 1, sum)))
    } else {
      label <- apply(genos, 1,
                     # Pick geno name with highest value; if >1, paste together.
                     function(x) paste(geno_names[which.max(x)], collapse = ""))
      color <- sapply(stringr::str_split(label, ""), 
                        function(x, alt) sum(x %in% alt), alt)
    }
  }
  list(label = label, color = color)
}