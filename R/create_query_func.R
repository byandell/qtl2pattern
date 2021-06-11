#' Create a function to query genotype probabilities
#'
#' Create a function that will connect to a database of genotype probability
#' information and return a list with `probs` object and a `map` object.
#'
#' @param dbfile Name of database file
#' @param method_val either \code{"fst"} or \code{"calc"} for type of genotype probabilities
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
#'     See \code{\link[qtl2]{create_variant_query_func}} for original idea.
#'
#' @examples 
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' create_qv <- function(dirpath) {
#'   # Download SNP info for DOex from web via RDS.
#'   # snpinfo is referenced internally in the created function.
#'   
#'   tmpfile <- tempfile()
#'   download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#'   snpinfo <- readRDS(tmpfile)
#'   unlink(tmpfile)
#'   snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#'   
#'   function(chr, start, end) {
#'     if(chr != "2") return(NULL)
#'     if(start < 96.5) start <- 96.5
#'     if(end > 98.5) end <- 98.5
#'     if(start >= end) return(NULL)
#'     dplyr::filter(snpinfo, .data$pos >= start, .data$pos <= end)
#'   }
#' }
#' query_variants <- create_qv(dirpath)
#' 
#' create_qg <- function(dirpath) {
#'   # Download Gene info for DOex from web via RDS
#'   # gene_tbl is referenced internally in the created function.
#'   
#'   tmpfile <- tempfile()
#'   download.file(file.path(dirpath, "c2_genes.rds"), tmpfile, quiet=TRUE)
#'   gene_tbl <- readRDS(tmpfile)
#'   unlink(tmpfile)
#'   
#'   function(chr, start, end) {
#'     if(chr != "2") return(NULL)
#'     if(start < 96.5) start <- 96.5
#'     if(end > 98.5) end <- 98.5
#'     if(start >= end) return(NULL)
#'     dplyr::filter(gene_tbl, .data$stop >= start, .data$start <= end)
#'   }
#' }
#' query_genes <- create_qg(dirpath)
#' 
#' # Examples for probs and mrna require either FST or RDS storage of data.
#' 
#' @export
create_probs_query_func <- function(dbfile,
                                    method_val = "fst",
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
#' @param mrnadir_val name of mRNA data directory (default \code{"RNAseq"})
#'
#' @return Function with seven arguments, `chr`, `start`,
#'     `stop`, `local`, `qtl`, `fast` and `mrnadir`. It returns a list with `expr`, `annot` and `peaks` objects
#'     spanning the region specified by the first three arguments.
#'
#' @details Note that this function assumes positions are in Mbp.
#'     There are required columns for each element, to be detailed in time.
#'     See \code{\link{read_mrna}} for details on how mRNA data are read.
#'     See \code{\link[qtl2]{create_variant_query_func}} for original idea.
#'
#' @export
create_mrna_query_func <- function(dbfile,
                                   mrnadir_val = "RNAseq") {
  if(missing(dbfile) || is.null(dbfile)) {
    # No mRNA data.
    function(chr = NULL, start = NULL, stop = NULL,
             local = TRUE,
             qtl = FALSE,
             mrnadir = mrnadir_val) 
      NULL
  } else {
    function(chr = NULL, start = NULL, stop = NULL,
             local = TRUE,
             qtl = FALSE,
             mrnadir = mrnadir_val)
      read_mrna(chr, start, stop, dbfile, local, qtl, mrnadir)
  }
}

# Dummy routines. See example above
query_genes <- function(...) {NULL}
query_variants <- function(...) {NULL}

