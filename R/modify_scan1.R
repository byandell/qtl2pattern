#' Modify scan1 lod matrix.
#'
#' @param scan1_object object of class \code{\link[qtl2scan]{scan1}}
#' @param newlod matrix of new LOD values
#' 
#' @export
#' 
modify_scan1 <- function(scan1_object, newlod) {
  x_attr <- attributes(scan1_object)
  x_attrnam <- names(x_attr)
  x_class <- class(scan1_object)
  
  for(obj in c("sample_size","SE", "hsq")) {
    if(obj %in% x_attrnam) {
      attr(newlod, obj) <- x_attr[[obj]]
    }
  }
  class(newlod) <- x_class
  newlod
}