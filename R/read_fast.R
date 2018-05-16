#' Read fast database with possible rownames
#' 
#' Read fast database with format \code{feather} or \code{fst}. Use first column
#' of database as rownames if desired. R/qtl2 routines assume data frames have
#' rownames to use to align individuals.
#' 
#' @param datapath character string path to database
#' @param columns names or indexes for columns to be extracted
#' @param rownames use first column of rownames if \code{TRUE} (can supply column number)
#' @param fast one of \code{"feather"} or \code{"fst"}
#' 
#' @export
#' @importFrom feather read_feather
#' @importFrom fst read_fst
#' @seealso \code{\link[fst]{read_fst}}, \code{\link[feather]{read_feather}}
#' 
read_fast <- function(datapath, columns = NULL, rownames = TRUE, fast = c("feather","fst")) {

  # Set up fast specific stuff.
  fast <- match.arg(fast)
  readfn <- switch(fast,
                   feather = feather::read_feather,
                   fst     = fst::read_fst)
  
  out <- as.data.frame(readfn(datapath, columns),
                       stringsAsFactors = FALSE)
  if(nrow(out) == 0)
    return(NULL)
  
  if(rownames) {
    rownames <- as.integer(rownames)
    if(fast== "fst") {
      # read_fst requires column name, not index.
      rownames <- colnames(readfn(datapath, from = 1, to = 1))[1]
    }
    # Row names (IDs) must be in first column of database (or use rownames as integer)
    # Row names not applied if any duplication.
    rowId <- unlist(readfn(datapath, rownames))
    if(length(rowId) == nrow(out) & !any(duplicated(rowId)))
      rownames(out) <- rowId
  }
  out
}