#' Compute distance matrix
#' 
#' @description
#' `compute_distmat()` computes and returns a distance matrix.
#'
#' @param data a numeric matrix or data frame of spatial coordinates.
#' @param method a character string of method `"geodesic"` or `"euclidean"`. Default is `"geodesic"`.
#' @param distGeo_a a numeric value of the major (equatorial) radius of
#' the ellipsoid for `"geodesic"` passed onto [geosphere::distGeo()].
#' Default value is for WGS84.
#' @param distGeo_f a numeric value of ellipsoid flattening for
#' `"geodesic"` passed onto [geosphere::distGeo()].
#' Default value is for WGS84.
#'
#' @details
#' [geosphere::distGeo()] is used internally when `"geodesic"` is specified.
#' Output unit is converted to km. 
#'
#' @returns a numeric distance matrix.
#' 
#' @seealso [geosphere::distm()],[geosphere::distGeo()]
#'
#' @export

compute_distmat <- function(data, method = "geodesic", distGeo_a = 6378137, distGeo_f = 1/298.257223563) {
  method <- match.arg(method, choices = c("geodesic", "euclidean"))
  data <- as.matrix(data)
  
  if (method == "geodesic") {
    distmat <- as.matrix(geosphere::distm(data, fun = function(x,y) geosphere::distGeo(p1 = x, p2 = y, a = distGeo_a, f = distGeo_f))) / 1000
  }
  if (method == "euclidean") {
    distmat <- as.matrix(stats::dist(data, method = "euclidean"))
  }
  return(distmat)
}
