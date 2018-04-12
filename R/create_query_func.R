#' Create a function to query genotype probabilities
#'
#' Create a function that will connect to a database of genotype probability
#' information and return a list with `probs` object and a `map` object.
#'
#' @param dbfile Name of database file
#' @param method_val either \code{"fast"} or \code{"calc"} for type of genotype probabilities
#' @param probdir_val name of probability directory (default \code{"genoprob"})
#'
#' @return Function with six arguments, `chr`, `start`,
#'     `stop`, `allele`, `method` and `probdir`. It returns a list with `probs` and `map` objects
#'     spanning the region specified by the first three arguments.
#'     The `probs` element should be either a `calc_genoprob` or `fst_genoprob` object
#'     (see \code{\link[qtl2fst]{fst_genoprob}}). 
#'
#' @details Note that this function assumes that \code{probdir_val} has a file with the
#'     physical map with positions in Mbp and other files with genotype probabilities.
#'     See \code{\link{read_probs}} for details on how probabilities are read.
#'     See \code{\link[qtl2]{create_variants_query_func}} for original idea.
#'
#' @export
#'
create_probs_query_func <- function(dbfile,
                                    method_val = "fast",
                                    probdir_val = "genoprob") {
  function(chr = NULL, start = NULL, stop = NULL,
           allele = TRUE,
           method = method_val,
           probdir = probdir_val) {
    read_probs(chr, start, stop, dbfile,
               allele,
               method = method,
               probdir = probdir)
  }
}
#' Create a function to query mRNA data
#'
#' Create a function that will connect to a database of mRNA information
#' and return a list with `probs` object and a `map` object.
#'
#' @param dbfile Name of database file
#' @param fast_val either \code{"fst"} or \code{"feather"} for type of mRNA database
#' @param mrnadir_val name of mRNA data directory (default \code{"RNAseq"})
#'
#' @return Function with seven arguments, `chr`, `start`,
#'     `stop`, `local`, `qtl`, `fast` and `mrnadir`. It returns a list with `expr`, `annot` and `peaks` objects
#'     spanning the region specified by the first three arguments.
#'
#' @details Note that this function assumes positions are in Mbp.
#'     There are required columns for each element, to be detailed in time.
#'     See \code{\link{read_mrna}} for details on how mRNA data are read.
#'     See \code{\link[qtl2]{create_variants_query_func}} for original idea.
#'
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
