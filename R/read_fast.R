#' Read fast database with possible rownames
#' 
#' Read fast database with format \code{fst}. Use first column
#' of database (must be named `ind`) as rownames if desired. R/qtl2 routines assume data frames have
#' rownames to use to align individuals.
#' 
#' @param datapath character string path to database
#' @param columns names or indexes for columns to be extracted
#' @param rownames use first column of rownames if \code{TRUE} (can supply column number)
#' 
#' @return extracted data frame with appropriate rows and columns.
#' 
#' @export
#' @importFrom fst read_fst
#' @seealso \code{\link[fst]{read_fst}}
#' 
read_fast <- function(datapath, columns = NULL, rownames = TRUE) {

  out <- as.data.frame(fst::read_fst(datapath, columns),
                       stringsAsFactors = FALSE)
  if(nrow(out) == 0)
    return(NULL)
  
  if(rownames) {
    # Row names (IDs) must be in database column named "ind".
    # Row names not applied if any duplication.
    rowId <- unlist(fst::read_fst(datapath, "ind"))
    if(length(rowId) == nrow(out) & !any(duplicated(rowId)))
      rownames(out) <- rowId
  }
  out
}