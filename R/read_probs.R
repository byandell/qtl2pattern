#' Read genotype probability object from file
#'
#' Read object from file stored according to method.
#'
#' @param chr vector of chromosome identifiers
#' @param start_val, end_val start and end values in Mbp
#' @param datapath name of folder with Derived Data
#' @param allele read haplotype allele probabilities (if \code{TRUE}) or diplotype allele-pair probabilities (if \code{FALSE})
#' @param method method of genoprob storage
#'
#' @return list with \code{probs} = large object of class \code{\link[qtl2geno]{calc_genoprob}} and \code{map} = physical map for selected \code{chr}
#'
#' @author Brian S Yandell, \email{brian.yandell@@wisc.edu}
#' @keywords utilities
#'
#' @examples
#' \dontrun{read_probs(chr, datapath)}
#'
#' @export
#' @importFrom qtl2feather subset_feather_genoprob
#' @importFrom stringr str_replace
#'
read_probs <- function(chr=NULL, start_val=NULL, end_val=NULL, datapath,
                       allele = TRUE, method = c("feather","calc")) {

  method <- match.arg(method)

  map <- readRDS(file.path(datapath, "pmap.rds"))

  if(is.null(chr)) {
    if(!allele & method == "calc")
      stop("must supply a chr")

    chr <- names(map)
  }

  map <- map[chr]

  probs <- switch(method,
                  feather = read_probs_feather(chr, datapath, allele),
                  calc    = read_probs_calc   (chr, datapath, allele))

  # Map may have extra markers. Trim.
  dnames <- dimnames(probs)$mar[chr]
  for(chri in chr) {
    wh <- match(dnames[[chri]], names(map[[chri]]))
    map[[chri]] <- map[[chri]][wh]
  }

  if(length(chr) == 1) {
    if(!is.null(start_val) & !is.null(end_val)) {

      ## Reduce to region of interest.
      wh <- which(map[[chr]] >= start_val &
                    map[[chr]] <= end_val)
      if(!length(wh))
        return(NULL)
      map[[chr]] <- map[[chr]][wh]

      switch(method,
             calc = {
               pr <- probs
               pr[[chr]] <- pr[[chr]][,,wh, drop = FALSE]
               probs <- modify_object(probs, pr)
             },
             feather = {
               probs <- qtl2feather::subset_feather_genoprob(probs, mar = wh)
             })
    }
  }

  list(probs = probs,
       map = map)
}
