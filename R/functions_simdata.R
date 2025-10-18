#------------------------------------------------------------------------------#
# Simulate data ####
#------------------------------------------------------------------------------#
#' Generate coordinates on a circle outline
#' 
#' @description
#' This function generates coordinates of equal distances on a circle outline centered at (0,0).
#'
#' @param radius a numeric value of radius from the center.
#' @param n an integer of the number of vertices. 
#'
#' @returns a numeric matrix of coordinates.
#'
#' @export

gen_circle_coords <- function(radius, n) {
  # Generate angles for each vertex
  angles <- seq(0, 2 * pi, length.out = n + 1)[-1]  # Exclude the last angle (2*pi)
  
  # Calculate x and y coordinates
  x_coords <- radius * sin(angles)
  y_coords <- radius * cos(angles)
  
  # Combine coordinates into a data frame
  coords <- as.matrix(data.frame(x = x_coords, y = y_coords))
  
  return(coords)
}

#------------------------------------------------------------------------------#
#' Generate Gaussian data clusters
#' 
#' @description
#' This function generates Gaussian data clusters.
#'
#' @param n an integer of the number of clusters
#' @param size an integer vector of sizes of a cluster. If an integer value is supplied, all clusters will be of this size.
#' @param center a data matrix or data frame of coordinates.
#' @param ... the optional argument includes `sigma`.
#' \itemize{
#'   \item \code{sigma}: a covariance matrix or a list of covariance matrices to specify different covariance matrices. See [mvtnorm::rmvnorm()].
#' }
#'
#' @returns a data frame of coordinates with clusters as a column.
#'
#' @seealso [mvtnorm::rmvnorm()]
#'
#' @export

gen_gaussian_data <- function(n, size, center = NULL, ...) {
  args <- list(...)
  sigma <- args$sigma
  
  if (length(size) == 1) {
    size <- rep(size, n)
  }
  
  if (is.null(center)) { center <- matrix(sample(x = seq(-10,10,0.01), size = 2*n, replace = T), nrow = n) }
  
  data <- list()
  for (i in 1:n) {
    if (is.null(sigma)) {
      data[[i]] <- mvtnorm::rmvnorm(n = size[i], mean = center[i,])
    } else {
      if (!is.list(sigma)) {
        data[[i]] <- mvtnorm::rmvnorm(n = size[i], mean = center[i,], sigma = sigma)
      } else {
        data[[i]] <- mvtnorm::rmvnorm(n = size[i], mean = center[i,], sigma = sigma[[i]])
      }
    }
    # label cluster
    data[[i]] <- cbind(data[[i]],i)
  }
  
  data <- do.call(rbind, data)
  colnames(data) <- c("x","y","clust")
  data <- as.data.frame(data)
  data$clust <- as.factor(data$clust)
  return(data)
}
