#' Compute distance matrix
#' 
#' @description
#' `compute_distmat()` computes and returns a distance matrix.
#'
#' @param data a numeric matrix or data frame of spatial coordinates.
#' @param method a string of the method. Either `"geodesic"` or `"euclidean"` is
#' available. Default is `"geodesic"`.
#'
#' @details
#' km is the unit when `"geodesic"` is specified.
#'
#' @returns a numeric distance matrix.
#' 
#' @seealso [geosphere::distm()]
#'
#' @export

compute_distmat <- function(data, method = "geodesic") {
  method <- match.arg(method, choices = c("geodesic", "euclidean"))
  data <- as.matrix(data)
  
  if (method == "geodesic") {
    distmat <- as.matrix(geosphere::distm(data, fun = geosphere::distGeo)) / 1000
  }
  if (method == "euclidean") {
    distmat <- as.matrix(stats::dist(data, method = "euclidean"))
  }
  return(distmat)
}
