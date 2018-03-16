#' @export
create_probs_query_func <- function(dbfile) {
  function(chr, start, stop,
           allele = TRUE,
           method = "feather",
           dirname = "genoprob") {
    read_probs(chr, start, stop, dbfile,
               allele,
               method = method,
               dirname = dirname)
  }
}
#' @export
create_mrna_query_func <- function(dbfile,
                                   fast = "feather",
                                   dirname = "RNAseq") {
  if(is.null(dbfile)) {
    # No mRNA data.
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE,
             fast = fast,
             dirname = dirname) 
      NULL
  } else {
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE,
             fast = fast,
             dirname = dirname)
      read_mrna(chr, start, stop, dbfile, local, qtl, fast, dirname)
  }
}
