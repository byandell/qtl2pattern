#' @export
create_probs_query_func <- function(dbfile,
                                    method_val = "feather",
                                    probdir_val = "genoprob") {
  function(chr, start, stop,
           allele = TRUE,
           method = method_val,
           probdir = probdir_val) {
    read_probs(chr, start, stop, dbfile,
               allele,
               method = method,
               probdir = probdir)
  }
}
#' @export
create_mrna_query_func <- function(dbfile,
                                   fast_val = "feather",
                                   mrnadir_val = "RNAseq") {
  if(missing(dbfile) || is.null(dbfile)) {
    # No mRNA data.
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE,
             fast = fast_val,
             mrnadir = mrnadir_val) 
      NULL
  } else {
    function(chr, start, stop,
             local = TRUE,
             qtl = FALSE,
             fast = fast,
             mrnadir = mrnadir)
      read_mrna(chr, start, stop, dbfile, local, qtl, fast, mrnadir)
  }
}
