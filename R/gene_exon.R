#' Get exons for set of genes
#'
#' Match up exon start,stop,strand with genes. Use \code{query_genes} to find features; see \code{\link[qtl2]{create_gene_query_func}}.
#'
#' @param top_snps_tbl table from \code{\link[qtl2]{top_snps}}
#' @param feature_tbl table of features from \code{query_genes}; see \code{\link[qtl2]{create_gene_query_func}}
#'
#' @return tbl of exon and gene features
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' \donttest{
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Download SNP info for DOex from web and read as RDS.
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
#' # Scan SNPs.
#' scan_snppr <- qtl2::scan1(snppr, DOex$pheno)
#' 
#' # Collect top SNPs
#' top_snps_tbl <- top_snps_pattern(scan_snppr, snpinfo)
#' 
#' # Download Gene info for DOex from web via RDS
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "c2_genes.rds"), tmpfile, quiet=TRUE)
#' gene_tbl <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' # Get Gene exon information.
#' out <- gene_exon(top_snps_tbl, gene_tbl)
#' summary(out, gene = out$gene[1])
#' }
#' 
#' @export
#' @rdname gene_exon
#' @importFrom dplyr arrange bind_rows desc distinct everything filter mutate select
#' @importFrom rlang .data
#' 
gene_exon <- function(top_snps_tbl,
                              feature_tbl = query_genes(chr_id, range_Mbp[1], range_Mbp[2])) {

  ## Only need distinct snp_id.
  top_snps_tbl <- dplyr::arrange(
    dplyr::select(
      dplyr::distinct(top_snps_tbl, .data$snp_id, .keep_all=TRUE),
      -.data$pheno),
    .data$pos)

  if(is.null(top_snps_tbl))
    return(NULL)
  if(!nrow(top_snps_tbl))
    return(NULL)

  chr_id <- as.character(unique(top_snps_tbl$chr))
  if(length(chr_id) != 1)
    stop("need exactly 1 chromosome in top_snps_tbl")
  range_Mbp <- range(top_snps_tbl$pos) + c(-1,1) * 0.005
  
  if(is.null(feature_tbl)) # Can happen if query_genes not supplied.
    return(NULL)
  
  gene_snp <- get_gene_snp(
    dplyr::select(
      top_snps_tbl, 
      .data$snp_id, .data$pos, .data$lod),
    feature_tbl)

  if(is.null(gene_snp))
    return(NULL)
  
  ## Need to get unique genes -- duplication with SNPs.
  gene_snp <- dplyr::distinct(gene_snp, .data$gene, .keep_all=TRUE)
  
  if(!nrow(gene_snp)) {
    return(NULL)
  }

  ## Use gene name as Name from feature_tbl.
  exons <- list()
  for(exoni in seq_len(nrow(gene_snp))) {
    genei <- gene_snp$gene[exoni]
    ## get genei and exons spanning genei
    exons[[genei]] <- dplyr::filter(
      dplyr::filter(feature_tbl,
                    .data$type %in% c("exon","gene"),
                    .data$start >= gene_snp$start[exoni],
                    .data$stop <= gene_snp$stop[exoni]),
      (!is.na(.data$Name) & .data$Name == genei) | .data$type == "exon")
    strandi <- gene_snp$strand[exoni]
    if(strandi != "." & !is.na(strandi))
      exons[[genei]] <- dplyr::filter(exons[[genei]], .data$strand == strandi)
  }
  out <- dplyr::distinct(
    dplyr::bind_rows(exons, .id = "gene"),
    .data$start, .data$stop, .data$strand, .keep_all=TRUE)
  class(out) <- c("gene_exon", class(out))

  out
}

#' Summary of exons for a gene with SNPs
#'
#' Returns table of gene and its exons.
#'
#' @param object Object of class \code{gene_exon}
#' @param ... additional parameters passed on to methods
#' @param gene_name name of gene as character string
#' @param top_snps_tbl table of top SNPs in region from \code{\link[qtl2]{top_snps}}
#' @param extra extra region beyond gene for SNPs (in Mbp)
#'
#' @return tbl of summary
#'
#' @keywords utilities
#'
#' @method summary gene_exon
#' @rdname gene_exon
#' @export
#' @importFrom dplyr arrange desc distinct filter group_by mutate n select summarize ungroup
#' 
summary.gene_exon <- function(object, gene_name=NULL,
                              top_snps_tbl = NULL,
                              extra = 0.005, ...) {
  ## Want to add columns for each phenotype
  ## with number of SNPs within extra of each gene.
  ## How to do this cleverly?
  if(!is.null(gene_name)) {
    out <- dplyr::distinct(
      dplyr::arrange(
        dplyr::select(
          dplyr::filter(object, .data$gene == gene_name),
          .data$gene, .data$source, .data$type, .data$start, .data$stop, .data$strand),
        dplyr::desc(.data$type)),
      .data$start, .data$stop, .data$strand, .keep_all=TRUE)
  } else {
    out <- dplyr::ungroup(
      dplyr::summarize(
        dplyr::group_by(
          dplyr::distinct(
            dplyr::filter(object, .data$type != "gene"),
            .data$start, .data$stop, .data$strand, .keep_all=TRUE),
          .data$gene),
        exons = dplyr::n(),
        min.len = min(.data$stop - .data$start),
        max.len = max(.data$stop - .data$start),
        sum.len = sum(.data$stop - .data$start),
        min_Mbp = min(.data$start),
        max_Mbp = max(.data$stop),
        strand = .data$strand[1]))
    out <- dplyr::select(
      out,
      .data$gene, .data$exons, .data$strand, .data$min_Mbp, .data$max_Mbp, dplyr::everything())
    
    if(!is.null(top_snps_tbl)) {
      ## Goal: add columns to out for each pheno in top_snps_tbl.
      ## Column should have number of SNPs within extra of gene.
      top_snps_tbl <- dplyr::select(top_snps_tbl, .data$pheno, .data$pos, .data$lod)
      pheno_names <- sort(unique(top_snps_tbl$pheno))
      outlim <- out[,c("min_Mbp","max_Mbp")]
      outlim[,1] <- outlim[,1] - extra
      outlim[,2] <- outlim[,2] + extra
      for(pheno_val in pheno_names) {
        out[[pheno_val]] <-
          apply(outlim, 1,
                function(x,y) {
                  in_region <- y$pos >= x[1] & y$pos <= x[2]
                  if(any(in_region))
                    max(y$lod[in_region])
                  else
                    NA
                },
                dplyr::filter(top_snps_tbl, .data$pheno == pheno_val))
      }
      tmp <- out[pheno_names]
      out$max_lod <- apply(tmp, 1, max)
      out <- dplyr::arrange(
        dplyr::select(
          out,
          .data$gene, .data$max_lod, dplyr::everything()),
        dplyr::desc(.data$max_lod))
    }
  }
  out
}
#' @param x Object of class \code{gene_exon}.
#' @param gene_val Name of gene from object \code{x}.
#' @param ... additional parameters passed on to methods
#' 
#' @rdname gene_exon
#' @export
subset.gene_exon <- function(x, gene_val, ...) {
  x <- dplyr::filter(x, .data$gene %in% gene_val)
  class(x) <- c("gene_exon", class(x))
  x
}

