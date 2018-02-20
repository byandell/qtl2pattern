#' @export
create_probs_query_func <- function(dbfile) {
  function(chr, start, stop, allele = TRUE) {
    read_probs(chr, start, stop, dbfile, allele, method = "feather")
  }
}
#' @export
create_mrna_query_func <- function(dbfile) {
  if(is.null(dbfile)) {
    # No mRNA data.
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE) 
      NULL
  } else {
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE)
      read_mrna(chr, start, stop, dbfile, local, qtl)
  }
}
