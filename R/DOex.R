#' Diversity Outbred Example
#' 
#' This is diversity outbred example data from 'qtl2data'.
#' 
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- qtl2::read_cross2(file.path(dirpath, "DOex.zip"))
#' DOex <- subset(DOex, chr = "2")
#' 
#' # Calculate genotype and allele probabilities
#' pr <- qtl2::calc_genoprob(DOex, error_prob=0.002)
#' 
#' # Download SNP info for DOex from web via RDS.
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_snpinfo.rds"), tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#' snpinfo <- dplyr::rename(snpinfo, pos = pos_Mbp)
#' 
#' # Convert to SNP probabilities
#' snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
#' snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
#' 
#' # Download Gene info for DOex from web via RDS
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_genes.rds"), tmpfile, quiet=TRUE)
#' genes <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' @return This has side actions of creating variables.
#'
DOex <- function() {
  example("DOex")
}


