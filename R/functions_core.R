require("dplyr")
require("magrittr")
require("sf") # for geometry
require("geosphere") # for calculating geodesic distance matrix
require("lhs") # for lhs sampling
require("moocore") # for multiobjective optimisation functions

#------------------------------------------------------------------------------#
# Algorithm ####
#------------------------------------------------------------------------------#
#' Geodesic distance computation
#' 
#' @description
#' This function calculates and returns the shortest distance in radians between a single (x,y) and all points (xa,ya) using spherical law of cosines.
#' @param x a numeric value of longitude.
#' @param y a numeric value of latitude.
#' @param ax a numeric vector of longitudes.
#' @param ay a numeric vector of latitudes.
#' @return a numeric vector of distances in radians.
#' @export

slc <- function(x, y, ax, ay, input = 'rad'){
  if(length(ax) != length(ay)) stop('ax must be same length as ay')
  if(length(ax) == 0) return(0); if(length(ax)>0){
    if(input=='deg'){x <- x*pi/180; y <- y*pi/180; ax <- ax*pi/180; ay <- ay*pi/180;}
    step.1 = sin(y) * sin (ay) + cos(y) * cos(ay) * cos(ax - x)
    
    # floating point bullshit, as sometimes 1 is greater than 1 (Rinferno!) 
    step.1[step.1 > 1] = 1
    step.1[step.1 < (-1)] = -1
    dist <- acos(step.1)	
    return(dist)}
}
#------------------------------------------------------------------------------#
#' Distance computation
#' 
#' @description
#' This function computes and returns either geodesic or Euclidean distances.
#' @param x a numeric value of longitude.
#' @param y a numeric value of latitude.
#' @param ax a numeric vector of longitudes.
#' @param ay a numeric vector latitudes.
#' @param method a string of method, This must be one of "geodesic" or "euclidean". Default is "geodesic".
#' @return a numeric vector of distances. Radian is the unit when "geodesic" is specified.
#' @export

compute_dist <- function(x, y, ax, ay, method = "geodesic") {
  if (method == "geodesic") {
    dist <- slc(x = x, y = y, ax = ax, ay = ay, input = 'deg')
  }
  if (method == "euclidean") {
    x <- c(x, ax)
    y <- c(y, ay)
    mat <- matrix(c(x,y), ncol = 2)
    dist <- as.numeric(as.matrix(dist(mat))[-1,1])
  }
  return(dist)
}
#------------------------------------------------------------------------------#
#' Distance matrix computation
#' 
#' @description
#' This function computes and returns the distance matrix computed using the specified method
#' to compute the distances between the rows of a data matrix. Km is the unit when "geodesic" is specified.
#' @param data a numeric matrix or data frame of two columns. The first and second columns are longitude and latitude respectively.
#' @param method a string of method. It must be one of "geodesic" or "euclidean". Default is "geodesic".
#' @return a numeric distance matrix. Km is the unit when "geodesic" is specified.
#' @export

compute_distmat <- function(data, method = "geodesic") {
  data <- as.matrix(data)
  
  if (method == "geodesic") {
    distmat <- as.matrix(geosphere::distm(data, fun = geosphere::distGeo))
  }
  if (method == "euclidean") {
    distmat <- as.matrix(dist(data, method = "euclidean"))
  }
  return(distmat)
}
#------------------------------------------------------------------------------#
#' Radial basis function
#' 
#' @description
#' This function kernelises the distance using the radial basis function. It is also known as the Gaussian kernel.
#' @param distance a numeric or complex vector of distances.
#' @param sigma a numeric value of sigma.
#' @return a numeric or complex vector of distances \eqn([0,1]) in the feature space.
#' @export

rbf <- function(distance, sigma) {
  exp(-(distance^2)/(2*sigma^2)) 
}
#------------------------------------------------------------------------------#
#' Spatial cost computation
#'
#' @description
#' This function computes and returns the statistic of the spatial cost function for a given cluster.
#' @param space_distmat a numeric spatial distance matrix.
#' @param clust_points a numeric vector of point indices in the targeted cluster.
#' @param sigma a numeric value of sigma passed on to [rbf()].
#' @return a numeric vector of spatial statistics for a given cluster.

compute_spacestat <- function(space_distmat, clust_points, sigma) {
  # transform the spatial distance matrix into a kernel matrix
  k_matrix <- rbf(space_distmat, sigma)
  
  # calculate the distance to the centroid in the feature space
  k_ii <- diag(k_matrix)
  k_ik <- k_matrix[ , clust_points, drop = FALSE]
  k_kk <- k_matrix[clust_points, clust_points, drop = FALSE]
  space_stat <- k_ii - 2*rowMeans(k_ik) + mean(k_kk)
  
  return(space_stat)
}

#------------------------------------------------------------------------------#
#' Intersect check
#' 
#' @description
#' This function checks and indicates if any clusters are intersecting based on their convex hulls.
#' @param data a numeric data matrix or data frame.
#' @param clust a numeric vector of cluster assignment. Default is NULL when it is already present as a column in the data.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()]. Default is 4326.
#' @return TRUE or FALSE
#' @export

intersects_bool <- function(data, clust = NULL, crs = 4326) {
  data <- as.data.frame(data)
  
  if (is.null(crs)) { crs <- NA } # Euclidean
  if (!is.null(clust)) { data$clust <- clust }
  
  data_sf <- sf::st_as_sf(data, coords = c(1,2), crs = crs)
  
  # compute concave hulls
  hulls <- data_sf %>%
    dplyr::group_by(clust) %>%
    dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
    sf::st_concave_hull(ratio = 1) # convex hulls
  
  # indicate T/F if there is any intersects
  intersects <- sf::st_intersects(hulls$geometry, sparse = F)
  diag(intersects) <- NA
  bool <- if (any(intersects == T, na.rm = T)) T else F
  
  return(bool)
}
#------------------------------------------------------------------------------#
#' Reorder clusters
#' 
#' @description
#' This function reorders and returns a vector of cluster assignment in ascending order.
#' @param clust a numeric vector of cluster assignment.
#' @return a numeric vector of reordered cluster assignment.

reorder_clust <- function(clust) {
  # e.g. c(1,4,4,2) will become c(1,2,2,3)
  sets <- list()
  for (i in 1:length(clust)) {
    sets[[i]] <- which(clust == clust[i])
    if (length(sets[[i]]) == 0) sets[[i]] <- NA
  }
  
  max_set_length <- max(sapply(sets, length))
  sets <- lapply(sets, function(x) if (length(x) < max_set_length) {c(x, rep(0, max_set_length - length(x)))} else {x})
  sets <- do.call(rbind, sets)
  sets <- unique(sets)
  sets <- as.matrix(sets[rowSums(sets) != 0, ]) # remove empty clusters
  
  # put NA to the back
  # count the number of rows with NA
  n_row <- nrow(sets)
  n_row_na <- length(unique(which(is.na(sets), arr.ind = T)[,"row"]))
  if (n_row_na > 0) {
    sets <- as.matrix(na.omit(sets))
    sets <- rbind(sets,NA)
  }
  
  # reorder cluster labels
  for (i in 1:nrow(sets)) {
    for (j in 1:ncol(sets)) {
      clust[sets[i,j]] <- i
    }
  }
  return(clust)
}

#------------------------------------------------------------------------------#
#' Clustering evaluation
#' 
#' @description
#' This function evaluates clustering performance and returns a summary and clusters that are below the critical size.
#' @param data a data frame.
#' @param space_distmat a numeric spatial distance matrix. Default is NULL.
#' @param space_distmethod a string of method used for distance computation. It must be one of "geodesic" or "euclidean". Default is "geodesic".
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()]. Default is 4326.
#' @param ... these arguments include `sigma`.
#' @param sigma a numeric value of sigma passed on to `rbf()`. It must be specified when `space_distmat` is supplied.
#' @details
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' 
#' The critical size of a cluster is defined as \eqn{\frac{N}{2K}} where \eqn{N} is the number of data points and \eqn{k} is the number of clusters.
#' @return a list of the following objects.
#' \itemize{
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{clust_below_size}: a numeric vector of clusters below the critical size.
#' }
#' @export

eval_blobs <- function(data, space_distmat = NULL, space_distmethod = "geodesic", crs = 4326, ...) {
  args <- list(...)
  sigma <- args$sigma
  
  # total number of points
  N <- nrow(data)
  # total number of clusters
  k <- length(unique(na.omit(data$clust))) # NA is excluded
  # initialise empty vectors 
  space_ss <- time_range <- time_evenness <- n <- numeric(k)
  
  # loop over k to obtain within cluster statistics
  for (j in 1:k) {
    clust_points <- which(data$clust == j)
    # if (length(clust_points) == 0) next ####
    data_k <- subset(data, data$clust == j)

    # spatial objectives
    # if "on the fly" mode is used
    if (is.null(space_distmat)) {
      space_centroids <- aggregate(data, list(data$clust), FUN = mean)[ ,c(2,3,5)] # 2,3 are spatial coordinates and 5 is clust
      space_centroids_k <- subset(space_centroids, space_centroids$clust == j)
      space_sq_distances <- vapply(1:nrow(data_k), function(i) compute_dist(x = data_k[i,1],
                                                                        y = data_k[i,2],
                                                                        ax = space_centroids_k[ ,1],
                                                                        ay = space_centroids_k[ ,2],
                                                                        method = space_distmethod),
                             numeric(1)) ^ 2
      space_ss[j] <- sum(space_sq_distances)
    
    # else if spatial distance matrix is supplied
    } else {
      # kernel k means cost function
      space_ss[j] <- sum(compute_spacestat(space_distmat = space_distmat, clust_points = clust_points, sigma = sigma)[clust_points])
    }
    # temporal objectives
    time_range[j] <- max(data_k[ ,3], na.rm = T) - min(data_k[ ,3], na.rm = T)
    time_evenness[j] <- 1 / (1 + var(diff(sort(data_k[ ,3])))) # NA if there are fewer than 3 data points
    # other constraint related statistics
    n[j] <- length(clust_points)
  }
  
  # Calculate the summary statistics
  space_wcss <- sum(space_ss)
  time_range_mean <- mean(time_range, na.rm = TRUE)
  time_range_sd <- sd(time_range, na.rm = TRUE)
  time_evenness_mean <- mean(time_evenness, na.rm = TRUE) # na.rm = TRUE to acess the overall performance as we normally do not care clusters with fewer than 3 points
  time_evenness_sd <- sd(time_evenness, na.rm = TRUE) # na.rm = TRUE to acess the overall performance as we normally do not care clusters with fewer than 3 points
  size_diff <- max(n, na.rm = TRUE) - min(n, na.rm = TRUE)
  size_mean <- mean(n, na.rm = TRUE)
  size_sd <- sd(n, na.rm = TRUE)
  
  # evaluate failed clust
  # n.points must be > N / k / 2
  clust_below_size <- which(n < N/k/2)
  n_fail <- sum(n[clust_below_size])
  
  # evaluate if blobs are intersecting in space
  intersects <- intersects_bool(data = data, crs = crs)
  
  # return a data frame of all the statistics
  summary <- data.frame(k = k,
                        space_wcss = space_wcss,
                        time_range_mean = time_range_mean,
                        time_range_sd = time_range_sd,
                        time_evenness_mean = time_evenness_mean,
                        time_evenness_sd = time_evenness_sd,
                        size_diff = size_diff,
                        size_mean = size_mean,
                        size_sd = size_sd,
                        n_fail = n_fail,
                        intersects = intersects)
  
  return(list(summary = summary, clust_below_size = clust_below_size))
}

#------------------------------------------------------------------------------#
#' Assign the starting cluster members - spatial distances computed on the fly
#' 
#' @description
#' Assign the starting cluster members by approximating the maximum spatial distant points. The distances are computed on the fly.
#' @param data a data matrix or data frame.
#' @param k an integer of the number of clusters
#' @param space_distmethod a string of the method for distance computation. It must be either "geodesic" or "euclidean". Default is "geodesic".
#' @return a data frame with assigned starting clusters as a column.

start_blobs_otf <- function(data, k, space_distmethod = "geodesic") {
  data <- as.data.frame(data)
  # assign an index to put data back into correct order at the end
  data$order <- 1:nrow(data)
  # initialise one point each to k blobs
  data$clust <- NA
  # a bit faster to handle data as a matrix
  data <- as.matrix(data)
  
  # start from k roughly equally spaced random locations 
  # (just permute a bunch and pick the one with the least smallest distance)
  N <- 100
  mat <- matrix( , N, k)
  mindist <- numeric(N)
  for(n in 1:N){
    i <- sort(sample(1:nrow(data), size = k))
    mat[n, ] <- i
    points <- data[i, 1:2]
    cb <- combn(k, 2)
    NC <- ncol(cb)
    distances <- numeric(NC)
    for(c in 1:NC) {
      distances[c] <- compute_dist(x = points[cb[1, c], 1],
                                   y = points[cb[1, c], 2],
                                   ax = points[cb[2, c], 1],
                                   ay = points[cb[2, c], 2],
                                   method = space_distmethod)
    } 
    mindist[n] <- min(distances)
  }
  # pick the set with largest minimum distance between two points to ensure maximum separation between k points;
  # pick the first one if two are tied 
  start <- mat[which(mindist == max(mindist))[1],]
  # Assign cluster memberships to the starting points
  data[start, "clust"] <- 1:k
  data <- as.data.frame(data)
  
  return(data)
}
#------------------------------------------------------------------------------#
#' Assign the starting cluster members - spatial distances supplied by a distance matrix
#' 
#' @description
#' Assign the starting cluster members by approximating the maximum spatial distant points. The distances are supplied by a distance matrix.
#' @param data a data matrix or data frame.
#' @param k an integer of the number of clusters.
#' @param space_distmat a numeric spatial distance matrix.
#' @return a data frame with assigned starting clusters as a column.

start_blobs_distmat <- function(data, k, space_distmat) {
  data <- as.data.frame(data)
  # assign an index to put data back into correct order at the end
  data$order <- 1:nrow(data)
  # initialise one point each to k blobs
  data$clust <- NA

  # a bit faster to handle data as a matrix
  data <- as.matrix(data)
  
  # start from k roughly equally spaced random locations 
  # (just permute a bunch and pick the one with the least smallest distance)
  N <- 100
  mat <- matrix(, N, k)
  mindist <- numeric(N)
  for(n in 1:N){
    i <- sort(sample(1:nrow(data), size = k)) # sample k points from the data, sort them 
    mat[n, ] <- i # store sorted index in the matrix by row
    cb <- combn(k, 2) # all combinations of k chooses 2 by column
    NC <- ncol(cb) # NC number of combinations
    distances <- numeric(NC) # a vector of distances of NC long
    for(c in 1:NC) distances[c] <- space_distmat[ i[cb[1, c]], i[cb[2, c]] ] # extract from distance matrix the distances for all combinations of points
    mindist[n] <- min(distances) # find the minimum
  }
  # Pick the set with largest minimum distance between two points to ensure maximum separation between k points;
  # pick the first one if two are tied 
  start <- mat[which(mindist == max(mindist))[1], ]
  # Assign cluster memberships to the starting points
  data[start,"clust"] <- 1:k
  data <- as.data.frame(data)
  
  return(data)
}
#------------------------------------------------------------------------------#
#' Assign the starting cluster members
#' 
#' @description
#' Assign the starting cluster members by approximating the maximum spatial distant points. A wrapper of `start_blobs_otf` and `start_blobs_distmat`.
#' @param data a data matrix or data frame.
#' @param k an integer of the number of clusters
#' @param space_distmat a numeric spatial distance matrix. Default is NULL.
#' @param space_distmethod a string of the method for distance computation. It must be either "geodesic" or "euclidean". Default is "geodesic".
#' @details
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' @return a data frame with assigned starting clusters as a column.
#' @export

start_blobs <- function(data, k, space_distmat = NULL, space_distmethod = "geodesic") {
  if (is.null(space_distmat)) {
    data <- start_blobs_otf(data = data, k = k, space_distmethod = space_distmethod)
  } else {
    data <- start_blobs_distmat(data = data, k = k, space_distmat = space_distmat)
  }
  
  return(data)
}
#------------------------------------------------------------------------------#
#' Assign clusters - spatial distances computed on the fly
#' 
#' @description
#' This function assigns clusters for a given k and r. Spatial distances are computed on the fly.
#' 
#' @param data a data matrix or data frame with assigned (starting) clusters as a column.
#' @param k an integer of the number of clusters.
#' @param r an numeric value of spatial relative weight. It must be \eqn{[0,1]}.
#' @param space_distmethod a string of the method for distance computation. It must be either "geodesic" or "euclidean". Default is "geodesic".
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' @return a data frame with assigned clusters as a column.

find_blobs_otf <- function(data, k, r, space_distmethod = "geodesic") {
  # a bit faster to handle data as a matrix
  data <- as.matrix(data)
  # for indexing speed, keep those assigned at the top
  data <- data[order(data[ ,"clust"]), ]
  # randomise the order of the points to be assigned (tba)
  tba_points <- which(is.na(data[ ,"clust"])) # to be assigned
  # for when iteration is implemented there will be no NA
  tba_points <- if (length(tba_points) == 0) 1:nrow(data) else tba_points
  tba_points <- sample(tba_points)
  a_points <- which(!1:nrow(data) %in% tba_points) # assigned
  # assigned at the top, randomised unassigned rows that follow
  data <- data[c(a_points,tba_points), ] 
  
  # reorder space_distmat to match the reordered data
  order <- data[ ,"order"]
  # extract clust to make the code cleaner
  clust <- data[,"clust"]
  
  # initialise stat matrices for the unassigned points
  start <- if (length(a_points) > nrow(data)) 1 else length(a_points) + 1
  N <- length(start:nrow(data))
  space_stat <- time_stat <- n <- numeric(k)
  
  # loop through every point (incremental updating)
  for (i in start:nrow(data)) {
    # loop through every k
    for (j in 1:k) {
      clust_points <- which(clust == j)
      n[j] <- length(clust_points)
      if (n[j] == 0) next
      data_k <- data[clust_points, , drop = FALSE]
      # spatial centroid
      space_centroid <- c(mean(data_k[ ,1]), mean(data_k[ ,2]))
      # compute the spatial cost
      space_stat[j] <- compute_dist(x = data[i,1], y = data[i,2], ax = space_centroid[1], ay = space_centroid[2], method = space_distmethod)
      # compute the temporal cost
      clust_points_tmp <- if (i %in% clust_points) clust_points[clust_points != i] else clust_points
      
      if (length(clust_points_tmp) == 0) {
        time_stat[j] <- 0  # only i was in the cluster
      } else {
        time_stat[j] <- min(abs(data[i, 3] - data[clust_points_tmp, 3]))
      }
    }
    
    # normalise the cost
    space_stat_norm <- (space_stat - min(space_stat, na.rm = T))/(max(space_stat, na.rm = T) - min(space_stat, na.rm = T))
    time_stat_norm <- (time_stat - min(time_stat, na.rm = T))/(max(time_stat, na.rm = T) - min(time_stat, na.rm = T))
    space_stat_norm[is.na(space_stat_norm)] <- 0
    time_stat_norm[is.na(time_stat_norm)] <- 0
    
    # scalarising (into a single objective statistic)
    obj <- space_stat_norm*r - time_stat_norm*(1-r)
    
    # assign cluster
    # check if there are tied clusters
    clust_tmp <- which(obj == min(obj, na.rm = T))
    
    if (length(clust_tmp) > 1) {
      # first check if one of those have fewer points
      n_clust_tmp <- n[clust_tmp]
      clust_tmp_minn_idx <- which(n_clust_tmp == min(n_clust_tmp))
      
      if (length(clust_tmp_minn_idx) > 1) {
        # randomly pick one if tied
        clust_tmp_minn_idx <- sample(clust_tmp_minn_idx, 1)
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      } else {
        # else pick the cluster with fewer assigned points
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      }
    } else {
      clust[i] <- clust_tmp
    }
  }
  
  # Put the data back to order
  data[ ,"clust"] <- clust
  data <- data[order(data[ ,"order"]), ]
  data <- as.data.frame(data)
  return(data)
}
#------------------------------------------------------------------------------#
#' Assign clusters - spatial distances supplied by a distance matrix
#' 
#' @description
#' This function assigns clusters for a given k and r. Spatial distances are supplied by a distance matrix.
#' 
#' @param data a data matrix or data frame with assigned (starting) clusters as a column.
#' @param k an integer of the number of clusters.
#' @param r an numeric value of spatial relative weight. It must be \eqn{[0,1]}.
#' @param space_distmat a numeric spatial distance matrix. Default is NULL.
#' @details
#' When distances are supplied by a distance matrix, kernel method is applied to compute cluster centroids in space.
#' @return a data frame with assigned clusters as a column.

find_blobs_distmat <- function(data, k, r, space_distmat, sigma) {
  # a bit faster to handle data as a matrix
  data <- as.matrix(data)

  # for indexing speed, keep those assigned at the top
  data <- data[order(data[ ,"clust"]), ]
  # randomise the order of the points to be assigned (tba)
  tba_points <- which(is.na(data[ ,"clust"])) # to be assigned
  # for when iteration is implemented there will be no NA
  tba_points <- if (length(tba_points) == 0) 1:nrow(data) else tba_points
  tba_points <- sample(tba_points)
  a_points <- which(!1:nrow(data) %in% tba_points) # assigned
  # assigned at the top, randomised unassigned rows that follow
  data <- data[c(a_points,tba_points), ] 
  
  # Reorder space_distmat/k.matrix to match the reordered data
  order <- data[ ,"order"]
  space_distmat <- space_distmat[order, order] 
  
  # transform the spatial distance matrix into a kernel matrix
  if (is.null(sigma)) stop("sigma is missing") # error message for missing sigma
  k_matrix <- rbf(space_distmat, sigma) 

  # Extract clust to make the code cleaner
  clust <- data[,"clust"]
  
  # initialise stat matrices for the unassigned points
  start <- if (length(a_points) > nrow(data)) 1 else length(a_points) + 1
  N <- length(start:nrow(data))
  space_stat <- time_stat <- n <- numeric(k)
  
  # loop through every point (incremental updating)
  for (i in start:nrow(data)) {
    # loop through every k
    for (j in 1:k) {
      clust_points <- which(clust == j)
      n[j] <- length(clust_points)
      if (n[j] == 0) next
      # compute the spatial cost == j
      # calculate the distance to the centroid in Hilbert space
      # O(c^2) is unavoidable 
      k_ii <- k_matrix[i, i] # included so the sigma makes more sense
      k_ik <- k_matrix[i, clust_points]
      k_kk <- k_matrix[clust_points, clust_points, drop = FALSE] # O(c^2)
      space_stat[j] <- k_ii - 2*mean(k_ik) + mean(k_kk) 

      # compute the temporal cost
      clust_points_tmp <- if (i %in% clust_points) clust_points[clust_points != i] else clust_points
      
      if (length(clust_points_tmp) == 0) {
        time_stat[j] <- 0  # only i was in the cluster
      } else {
        time_stat[j] <- min(abs(data[i, 3] - data[clust_points_tmp, 3]))
      }
      
    }
    
    # normalise the cost
    space_stat_norm <- (space_stat - min(space_stat, na.rm = T))/(max(space_stat, na.rm = T) - min(space_stat, na.rm = T))
    time_stat_norm <- (time_stat - min(time_stat, na.rm = T))/(max(time_stat, na.rm = T) - min(time_stat, na.rm = T))
    space_stat_norm[is.na(space_stat_norm)] <- 0
    time_stat_norm[is.na(time_stat_norm)] <- 0
    
    # scalarising (into a single objective statistic)
    obj <- space_stat_norm*r - time_stat_norm*(1-r)
    
    # assign cluster
    # check if there are tied clusters
    clust_tmp <- which(obj == min(obj, na.rm = T))
    
    if (length(clust_tmp) > 1) {
      # first check if one of those have fewer points
      n_clust_tmp <- n[clust_tmp]
      clust_tmp_minn_idx <- which(n_clust_tmp == min(n_clust_tmp))
      
      if (length(clust_tmp_minn_idx) > 1) {
        # randomly pick one if tied
        clust_tmp_minn_idx <- sample(clust_tmp_minn_idx, 1)
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      } else {
        # else pick the cluster with fewer assigned points
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      }
    } else {
      clust[i] <- clust_tmp
    }
  }
  
  # Put the data back to order
  data[ ,"clust"] <- clust
  data <- data[order(data[ ,"order"]), ]
  data <- as.data.frame(data)
  return(data)
}

#------------------------------------------------------------------------------#
#' Assign clusters
#' 
#' @description
#' This function assigns clusters for a given k and r.
#' 
#' @param data a data matrix or data frame with assigned (starting) clusters as a column.
#' @param k an integer of the number of clusters.
#' @param r an numeric value of spatial relative weight. It must be \eqn{[0,1]}.
#' @param space_distmat a numeric spatial distance matrix. Default is NULL.
#' @param space_distmethod a string of the method for distance computation. It must be either "geodesic" or "euclidean". Default is "geodesic".
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' 
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' In that case, kernel method is applied to compute cluster centroids in space.
#' @return a data frame with assigned clusters as a column.
#' @export

find_blobs <- function(data, k, r, space_distmat = NULL, space_distmethod = "geodesic", ...) {
  args <- list(...)
  sigma <- args$sigma
  
  if (is.null(space_distmat)) {
    data <- find_blobs_otf(data = data, k = k, r = r, space_distmethod = space_distmethod)
  } else {
    data <- find_blobs_distmat(data = data, k = k, r = r, space_distmat = space_distmat, sigma = sigma)
  }
  return(data)
}
#------------------------------------------------------------------------------#
#' Compare two cluster assignments
#' 
#' @description
#' This function compares two cluster assignments and retains the shared members.
#' 
#' @param b1 a data frame with assigned clusters as a column.
#' @param b2 a data frame with assigned clusters as a column.
#' @return a data frame with shared assigned clusters as a column.

compare_blobs <- function(b1, b2){
  # compare blobs and retain the commonality
  # b1 and b2: two outputs of find_blobs() using the same dataset
  
  k1 <- length(unique(na.omit(b1$clust))) # account for non-assigned points, i.e. clust == NA
  k2 <- length(unique(na.omit(b2$clust)))
  k <- max(k1,k2)
  
  # all the possible pairs of cluster assignment
  pairs <- expand.grid(b1=1:k, b2=1:k)
  # Number of pairs
  N <- nrow(pairs)
  # initialise a empty vector of match rate
  match <- numeric(N)
  # calculate match rate for each pair
  for(n in 1:N){
    v1 <- subset(b1, clust == pairs$b1[n])$order # get the index
    v2 <- subset(b2, clust == pairs$b2[n])$order # get the index
    t1 <- v1 %in% v2 
    t2 <- v2 %in% v1
    match[n] <- (sum(t1) + sum(t2)) / (length(t1) + length(t2))
  }
  pairs$match <- match
  # remove pairs where match rate is 0
  pairs <- subset(pairs, match != 0)
  # order pairs by match rate
  pairs <- pairs[order(pairs$match, decreasing=TRUE), ]
  # ensure unique pairs of length k by removing duplicated label at k+1 
  pairs <- pairs[!duplicated(pairs$b1) & !duplicated(pairs$b2), ]
  N <- nrow(pairs)
  # update labels of b2
  if (N > 0) {
    for(n in 1:N) b2$newclust[b2$clust == pairs$b2[n]] <- pairs$b1[n]
  }
  
  # for NAs, copy the original labels
  b2$newclust[is.na(b2$newclust)] <- b2$clust[is.na(b2$newclust)]
  
  # store matches between blobs
  same <- b1$clust == b2$newclust
  same[is.na(same)] <- FALSE
  b <- b1
  b$clust[!same] <- NA
  
  # ensure all clusters have at least one point assigned
  missing <- (1:k)[!(1:k) %in% unique(b$clust)]
  M <- length(missing)
  L_na <- length(which(is.na(b$clust)))
  if(M > 0 & L_na > 0) {
    for(m in 1:M){
      i <- sample(which(is.na(b$clust)), size=1)
      b$clust[i] <- missing[m]
    }
  }
  return(b)
}

#------------------------------------------------------------------------------#
#' Core local-search algorithm - "iterative" method for convergence
#' 
#' @description
#' This function performs the core bi-objective optimisation algorithm to assign clusters for a given k and r.
#' Convergence is directed by the "iterative" method.
#' 
#' @param data a data matrix or data frame with assigned (starting) clusters as a column.
#' @param k an integer of the number of clusters.
#' @param r an numeric value of spatial relative weight. It must be \eqn{[0,1]}.
#' @param iter an integer of the number of iterations. Default is 3.
#' @param space_distmat a numeric spatial distance matrix. Default is NULL.
#' @param space_distmethod a string of the method for distance computation. It must be either "geodesic" or "euclidean". Default is "geodesic".
#' @param converge_ari a numeric value of Adjusted Rand Index (ARI) to indicate convergence threshold between two searches. It must be \eqn{[0,1]}. Default is NULL.
#' @param random_start a logical operator. Should random start or [start_blobs()] be used? Default is F.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()]. Default is 4326.
#' @param ... these arguments include `sigma`.
#' @param sigma a numeric value of sigma passed on to [rbf()]. It must be specified when `space_distmat` is supplied.
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' 
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' In that case, kernel method is applied to compute cluster centroids in space.
#' 
#' When the "iterative" method is specified, one search is performed in every iteration. It runs until the last iteration or convergence. 
#' When `converge_ari` is specified, convergence is defined and activated when ARI between the latest and the previous search is
#' above the specified threshold and at least three iterations are run.
#' @return a list of the following objects.
#' \itemize{
#'   \item \code{data}: a data frame of the input data with assigned clusters as a column.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{clust_below_size}: a numeric vector of clusters below the critical size.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#' }

blob_search_iter <- function(data, k, r, iter = 3, space_distmat = NULL,
                             space_distmethod = "geodesic", converge_ari = NULL, random_start = F, crs = 4326, ...) {
  
  if (random_start == T) {
    data <- as.data.frame(data)
    data$clust <- NA
    start <- sample(1:nrow(data), k)
    data$clust[start] <- 1:k
  } else {
    # start_blobs() to pick centroids
    data <- start_blobs(data = data, k = k, space_distmat = space_distmat, space_distmethod = space_distmethod)
  }
  # initialise counter counting find_blobs()
  t <- 0
  
  # intialise trace table
  trace_df <- data.frame()
  
  for (i in 1:iter) {
    
    data_old <- data
    # find_blobs()
    data <- find_blobs(data = data, k = k, r = r, space_distmat = space_distmat, space_distmethod = space_distmethod, ...)
    # count find_blobs() executed
    t <- t + 1
    
    if (t > 0) {
      # check convergence
      ari <- mclust::adjustedRandIndex(data$clust, data_old$clust)
      data_common <- compare_blobs(data, data_old)
      n_common <- sum(!is.na(data_common$clust))
      
      # eval_blobs()
      eval_out <- eval_blobs(data = data, space_distmat = space_distmat, space_distmethod = space_distmethod, crs = crs, ...)
      clust_below_size <- eval_out$clust_below_size
      trace_df_newrow <- eval_out$summary
      trace_df_newrow$iter <- t
      trace_df_newrow$ari <- ari
      trace_df <- rbind(trace_df, trace_df_newrow)
      
      # if converged between t and t-1, break
      if (!is.null(converge_ari)) {
        if (all(ari >= converge_ari) == TRUE & t >= 3) break
      }
    }
  }
  
  trace_df$r <- r
  summary <- trace_df[nrow(trace_df),]
  summary$iter_method <- "iterative"
  
  # order is not useful in the result
  # data <- subset(data, select = -order)
  data$order <- NULL
  
  return(list(data = data, summary = summary, clust_below_size = clust_below_size, trace = trace_df))
}

#------------------------------------------------------------------------------#
#' Core local-search algorithm - "commonality" method for convergence
#' 
#' @description
#' This function performs the core bi-objective optimisation algorithm to assign clusters for a given k and r.
#' Convergence is directed by the "commonality" method.
#' 
#' @inheritParams blob_search_iter
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' 
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' In that case, kernel method is applied to compute cluster centroids in space.
#' 
#' When the "commonality" method is specified, in every iteration, three searches are run and the shared members are retained for the next iteration.
#' After the last iteration or when searches converge, the search with the most shared members will be used to assign clusters.
#' When `converge_ari` is specified, convergence is defined and activated when all ARIs between all two pairs of searches is
#' above the specified threshold and at least three iterations are performed. 
#' @inherit blob_search_iter returns

blob_search_common <- function(data, k, r, iter = 3, space_distmat = NULL,
                               space_distmethod = "geodesic", converge_ari = NULL, random_start = F, crs = 4326, ...) {
  
  if (random_start == T) {
    data_list <- list()
    for (i in 1:3) {
      data_list[[i]] <- as.data.frame(data)
      data_list[[i]]$clust <- NA
      start <- sample(1:nrow(data), k)
      data_list[[i]]$clust[start] <- 1:k
    }
  } else {
    # starting 3 blobs
    # start_blobs() to pick centroids
    data_list <- lapply(1:3, function(i) { start_blobs(data = data, k = k, space_distmat = space_distmat, space_distmethod = space_distmethod) })
  }
  
  # pairs to compare
  pairs <- combn(3, 2, simplify = FALSE)
  # initialise counter counting find.blobs()
  t <- 0
  
  # initialise empty list for trace 
  trace_list <- list()
  
  for (i in 1:iter) {
    # find_blobs()
    data_list <- lapply(data_list, function(x) { find_blobs(data = x, k = k, r = r, space_distmat = space_distmat, space_distmethod = space_distmethod, ...) })
    t <- t + 1
    
    # check convergence
    # Adjusted Rand Index [0 (random),1 (identical)]
    ari <- sapply(pairs, function(x) {
      mclust::adjustedRandIndex(data_list[[ x[1] ]][ , "clust"], data_list[[ x[2] ]][ , "clust"])
    })
    
    # eval_blobs()
    if (trace == T) {
      trace_list[[i]] <- lapply(data_list, function(x) eval_blobs(data = x, space_distmat = space_distmat, space_distmethod = space_distmethod, ...)$summary)
      trace_list[[i]] <- do.call(rbind, trace_list[[i]])
      trace_list[[i]]$iter <- t
    }
    
    # compare_blobs(), retain commonality
    data_list <- lapply(pairs, function(x) { compare_blobs(data_list[[ x[1] ]], data_list[[ x[2] ]]) })
    
    # extract the number of common points
    n_common <- sapply(data_list, function(x) sum(!is.na(x$clust)))
    
    # further fields to be added to the trace table
    trace_list[[i]]$ari <- paste(ari, collapse = ",")
    # trace_list[[i]]$n_common <- paste(n_common, collapse = ",")
    
    # if converged between candidates, break
    if (!is.null(converge_ari)) {
      # ari >= a certain similarity and iteration >= 3
      if (all(ari >= converge_ari) == TRUE & t >= 3) break
    }
  }
  
  best_idx <- which.max(n_common)
  data <- find_blobs(data_list[[best_idx]], k = k, r = r, space_distmat = space_distmat, space_distmethod = space_distmethod, ...)
  n_common <- n_common[best_idx]
  
  eval_out <- eval_blobs(data = data, space_distmat = space_distmat, space_distmethod = space_distmethod, crs = crs, ...)
  clust_below_size <- eval_out$clust_below_size
  summary <- eval_out$summary
  summary$iter <- t
  # summary$n_common <- n_common
  summary$ari <- paste(ari, collapse = ",")
  summary$iter_method <- "commonality"
  summary$r <- r
  
  trace_df <- do.call(rbind, trace_list)
  trace_df$r <- r
  
  # order is not useful in the result
  # data <- subset(data, select = -order)
  data$order <- NULL
  
  return(list(data = data, summary = summary, clust_below_size = clust_below_size, trace = trace_df))
}

#------------------------------------------------------------------------------#
#' Core local-search algorithm
#' 
#' @description
#' This function performs the core bi-objective optimisation algorithm to assign clusters for a given k and r.
#' 
#' @inheritParams blob_search_iter
#' @param filter_intersects a logical operator. Should an assignment with intersects in space be removed? Default is T.
#' @param filter_clustsize a logical operator. Should a cluster below a critical size be assigned NA? Default is T.
#' @param max_na a numeric value of the maximum proportion of NAs allowed. It must be \eqn{[0,1]}. Default is 0.05.
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' 
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' In that case, kernel method is applied to compute cluster centroids in space.
#' 
#' When the `"iterative"` method is specified in `iter_method`, one search is performed in every iteration. It runs until the last iteration or convergence. 
#' When `converge_ari` is specified, convergence is defined and activated when ARI between the latest and the previous search is
#' above the specified threshold and at least three iterations are run.
#' 
#' When the `"commonality"` method is specified in `iter_method`, in every iteration, three searches are run and the shared members are retained for the next iteration.
#' After the last iteration or when searches converge, the search with the most shared members will be used to assign clusters.
#' When `converge_ari` is specified, convergence is defined and activated when all ARIs between all two pairs of searches is
#' above the specified threshold and at least three iterations are performed. 
#' 
#' The critical size of a cluster is defined as \eqn{\frac{N}{2K}} where \eqn{N} is the number of data point and \eqn{k} is the number of clusters.
#' @return a list of the following objects.
#' \itemize{
#'   \item \code{data}: a data frame of the input data with assigned clusters in a column.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#' }
#' @seealso [sf::st_as_sf()]
#' @export

blob_search <- function(data, k, r, iter = 3, space_distmat = NULL, space_distmethod = "geodesic",
                        iter_method = "iterative", converge_ari = NULL, random_start = F, crs = 4326,
                        filter_intersects = T, filter_clustsize = T, max_na = 0.05, ...) {
  N <- nrow(data)
  
  if (iter_method == "iterative") {
    blob <- blob_search_iter(data = data, k = k, r = r, iter = iter,
                             space_distmat = space_distmat, space_distmethod = space_distmethod,
                             converge_ari = converge_ari, random_start = random_start, crs = crs, ...)
  }
  
  if (iter_method == "commonality") {
    blob <- blob_search_common(data = data, k = k, r = r, iter = iter,
                               space_distmat = space_distmat, space_distmethod = space_distmethod,
                               converge_ari = converge_ari, random_start = random_start, crs = crs, ...)
  }
  
  # filter for constraints
  if (filter_intersects == T) {
    intersects <- blob$summary$intersects
    if (intersects == T) {
      # message("At least two clusters intersect. Return NULL.")  
      return(1)
    }
  }
  
  # initialise n.points.removed column. If no point is removed, the entries will all be 0
  blob$summary$n_removed <- 0
  
  if (filter_clustsize == T) {
    clust_below_size <- blob$clust_below_size
    if (length(clust_below_size) > 0) {
      # assign NA to removed clusters
      n_removed <- blob$summary$n_fail
      blob$data$clust[which(blob$data$clust %in% clust_below_size)] <- NA
      # reassign the cluster number
      # make sure there is no gap in the sequence of k
      blob$data$clust <- reorder_clust(blob$data$clust)
      # update eval.res$summary
      eval_out_updated <- eval_blobs(data = blob$data, space_distmat = space_distmat, space_distmethod = space_distmethod, ...)
      updated_cols <- intersect(names(blob$summary), names(eval_out_updated$summary))
      blob$summary[ , updated_cols] <- eval_out_updated$summary[ , updated_cols]
      blob$summary$n_removed <- n_removed
    }
  } else {
    # reassign the cluster number
    # make sure there is no gap in the sequence of k
    blob$data$clust <- reorder_clust(blob$data$clust)
  }
  
  # return NULL if too many points are removed
  if (blob$summary$n_removed > N * max_na) return(2)
  if (blob$summary$k < 2) return(3)
  
  blob$summary$n_fail <- NULL # confusing to have both n_removed and this in the output
  blob$trace$n_fail <- NULL # confusing to have both n_removed and this in the output
  
  # remove clust.below.size as filter has been applied so it is irrelevant to the downstream
  blob$clust_below_size <- NULL

  return(blob)
}

#------------------------------------------------------------------------------#
# Populate solutions ####
#------------------------------------------------------------------------------#
#' Find duplicates
#' 
#' @description
#' This function find duplicates in a set of cluster assignments.
#' 
#' @param clust a numeric matrix of cluster assignments. Each row is an assignment.
#' @param ari a numeric value of Adjusted Rand Index (ARI) to indicate duplication threshold between two assignments. It must be \eqn{[0,1]}. Default is 1.
#' @return a list of the following objects.
#' \itemize{
#'   \item \code{idx}: a numeric vector of duplicate indices.
#'   \item \code{pairwise_ari}: a data frame of summary statistics.
#'   \item \code{freq}: a contingency table of the frequency of duplicates, an object of class "table", an array of integer values. See also [base::table()].
#'   \item \code{pairs_dup}: a numeric matrix of pairs of duplicates. The second row stores indices of duplicates of the first row.
#' }
#' @seealso [base::table()]
#' @export

find_dup <- function (clust, ari = 1) {
  # As NA is ignored by mclust::adjustedRandIndex()
  clust[is.na(clust)] <- 0
  
  # number of solutions
  n_sol <- nrow(clust)
  # get all combination of pairs
  pairs <- combn(n_sol, 2)
  # calculate the pairwise ARI for all pairs of solutions
  pairwise_ari <- future.apply::future_vapply(1:ncol(pairs),
                                              function(x) mclust::adjustedRandIndex(clust[pairs[1,x], ], clust[pairs[2,x], ]),
                                              numeric(1),
                                              future.seed = T)
  
  # index the pairs with ari >= ari
  pairs_dup_idx <- which(pairwise_ari >= ari)
  # subset the columns of duplicated pairs
  pairs_dup <- pairs[, pairs_dup_idx, drop = F]
  
  # if a ~ b and b ~ c, then not necessarily a ~ c
  # the following steps makes sure they are the pairs with ari >= ari.dup 
  if (ncol(pairs_dup) > 1) {
    # dependence boolean
    pairs_dup_depend_bool <- logical()
    for (i in 2:ncol(pairs_dup)) {
      pairs_dup_depend_bool[i] <- pairs_dup[1, i] %in% pairs_dup[2, 1:(i-1)]
    }
    pairs_dup_depend_idx <- which(pairs_dup_depend_bool == T)
    # remove the dependent columns 
    if (length(pairs_dup_depend_idx) > 0) {
      pairs_dup <- pairs_dup[ , -pairs_dup_depend_idx, drop = F]
    }
  }
  
  # note it is sensitive to order
  # second row is considered the duplicate of the first row
  idx <- unique(pairs_dup[2, ])
  
  # frequency of duplicates
  freq <- table(pairs_dup[1,])
  
  return(list(idx = idx, pairwise_ari = pairwise_ari, freq = freq, pairs_dup = pairs_dup))
}
#------------------------------------------------------------------------------#
#' Populate solutions by weighted sum scalarisation
#' 
#' @description
#' This function populate solutions by weighted sum scalarisation of the bi-objective function in [blob_search()] for a given k. 
#' 
#' @inheritParams blob_search_iter
#' @inheritParams blob_search
#' @params r_range a numeric vector of length 2 indicating the lower and upper bounds of the relative spatial weight. They must be \eqn{[0,1]}.
#' @params run an integer of the number of runs.
#' @details
#' When distances are computed on the fly, kernel method cannot be applied to compute cluster centroids in space.
#' 
#' When `space_distmat` is specified. Distance matrix implementation will ignore `space_distmethod`.
#' In that case, kernel method is applied to compute cluster centroids in space.
#' 
#' When the `"iterative"` method is specified in `iter_method`, one search is performed in every iteration. It runs until the last iteration or convergence. 
#' When `converge_ari` is specified, convergence is defined and activated when ARI between the latest and the previous search is
#' above the specified threshold and at least three iterations are run.
#' 
#' When the `"commonality"` method is specified in `iter_method`, in every iteration, three searches are run and the shared members are retained for the next iteration.
#' After the last iteration or when searches converge, the search with the most shared members will be used to assign clusters.
#' When `converge_ari` is specified, convergence is defined and activated when all ARIs between all two pairs of searches is
#' above the specified threshold and at least three iterations are performed. 
#' 
#' The critical size of a cluster is defined as \eqn{\frac{N}{2K}} where \eqn{N} is the number of data point and \eqn{k} is the number of clusters.
#' 
#' Scalarisation is achieved by varying the relative spatial weight generated by Latin hypercube sampling using [lhs::randomLHS()].
#' 
#' To parallelise runs, [future.apply::future_lapply()] is implemented. See [future::future] and [future.apply::future.apply] for more information.
#' 
#' @return a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of cluster assignments. Each row is a solution.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#'   \item \code{n_filtered}: a data frame of numbers of filtered solutions.
#' }
#' @seealso [sf::st_as_sf()], [lhs::randomLHS()], [future::future], [future.apply::future.apply]
#' @export
blob_populate <- function (data, k, r_range = c(0.5,1), iter = 3, run = 10, space_distmat = NULL,
                           space_distmethod = "geodesic", iter_method = "iterative", converge_ari = NULL,
                           random_start = F, crs = 4326,
                           filter_intersects = T, filter_clustsize = T, max_na = 0.05, ...) {
  # r_range: a vector of 2L indicating the lower and upper limit of relative weight
  
  if (length(r_range) > 1) {
    # LHS sampling for more evenly distributed parameters
    lhs_samples <- lhs::randomLHS(run,1)
    # Scale to the range
    r_samples <- sort(as.vector(min(r_range) + lhs_samples * (max(r_range) - min(r_range))))
  } else {
    r_samples <- rep(r_range, run)
  }
  
  pop <- future.apply::future_lapply(r_samples, function (r) {
    blob_search(data = data, k = k, r = r, iter = iter, space_distmat = space_distmat,
                space_distmethod = space_distmethod, iter_method = iter_method,
                converge_ari = converge_ari, random_start = random_start,
                filter_intersects = filter_intersects, filter_clustsize = filter_clustsize, max_na = max_na, ...)
  }, future.seed = T)
  
  # count the runs that are filtered out, informative about the range of r to be searched
  filtered_intersects <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 1)
  filtered_intersects <- sapply(filtered_intersects, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_intersects <- sum(filtered_intersects)
  
  filtered_size <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 2)
  filtered_size <- sapply(filtered_size, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_size <- sum(filtered_size)
  
  filtered_k1 <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 3)
  filtered_k1 <- sapply(filtered_k1, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_k1 <- sum(filtered_k1)
  
  # if blob returns 1 or 2, set NULL so rbind in the following step will omit them
  pop <- lapply(pop, function(x) if(!is.list(x)) NULL else x)
  # extract each element and add a column to indicate the run
  pop <- do.call(rbind, pop)
  data_list <- lapply(seq_along(pop[, "data"]), function (i) cbind(pop[, "data"][[i]], run = i))
  summary_list <- lapply(seq_along(pop[, "summary"]), function (i) cbind(pop[, "summary"][[i]], run = i))
  trace_list <- lapply(seq_along(pop[, "trace"]), function (i) cbind(pop[, "trace"][[i]], run = i))
  
  # extract clust and append the run in the beginning of the vector
  clust_list <- lapply(seq_along(data_list), function (i) data_list[[i]]$clust)
  
  # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
  clust <- do.call(rbind, clust_list)
  summary <- do.call(rbind, summary_list)
  trace <- do.call(rbind, trace_list)
  
  # find and filter exact duplicates
  filtered_dup <- 0
  
  if (!is.null(clust)) {
    summary$dup <- 0
    # Here, clust <- do.call(rbind, clust_list) must return a matrix even if there is only one solution
    if (nrow(clust) > 1) {
      dup <- find_dup(clust, ari = 1)
      
      if (length(dup$idx) > 0) {
        # record the duplicate freq
        summary$dup[as.numeric(names(dup$freq))] <- as.vector(dup$freq)
        # filter the duplicates
        clust <- clust[-dup$idx, , drop = F]
        summary <- summary[-dup$idx, ]
        trace <- trace[-which(trace$run %in% dup$idx),]
        # record the total no. of dup filtered
        filtered_dup <- length(dup$idx)
      }
    }
  }
  
  n_filtered <- data.frame(intersects = filtered_intersects, size = filtered_size, k1 = filtered_k1, dup = filtered_dup)

  pop <- list(clust = clust, summary = summary, trace = trace, n_filtered = n_filtered)

  return(pop)
}
#------------------------------------------------------------------------------#
# Optimising number of k ####
#------------------------------------------------------------------------------#
# It works by looping from k = 2 to sqrt(N) 
blob_kpopulate <- function(data, k_range = NULL, r_range = c(0.5,1), iter = 3, run = 10, 
                           space_distmat = NULL, space_distmethod = "geodesic",
                           iter_method = "iterative", converge_ari = NULL,
                           random_start = F, filter_intersects = T, filter_clustsize = T, max_na = 0.05, ...) {
  
  # iter: iteration for each k
  # run: run for each k
  
  # number of samples
  N <- nrow(data)
  if (is.null(k_range)) {
  	k_range <- integer(2)
  	k_range[1] <- 2 # lower bound
    k_range[2] <- floor(sqrt(N)) # upper bound, an heuristic to maximise information e.g. 100 points: 10 blobs, 10 points each 
  }
  
  pop_list <- list()
  for (k in k_range[1]:k_range[2]) {
    pop_list[[k]] <- blob_populate(data = data, k = k, r_range = r_range, iter = iter, run = run, 
                              space_distmat = space_distmat, space_distmethod = space_distmethod,
                              iter_method = iter_method, converge_ari = converge_ari,
                              random_start = random_start, filter_intersects = filter_intersects, filter_clustsize = filter_clustsize, max_na = max_na, ...)
    
    # if no solution for a given k, stop the loop
    if (is.null(pop_list[[k]]$clust)) {
      pop_list[k] <- NULL
      next # may be next vs break? assuming too overlap to pursue further
    }
  }
  
  pop_list <- pop_list[-1] # as the loop starts at 2
  
  # extract each element and add a column to indicate the initial k
  pop <- do.call(rbind, pop_list)
  
  summary_list <- lapply(seq_along(pop[ , "summary"]), function (i) cbind(pop[ , "summary"][[i]], k_o = i + 1))
  trace_list <- lapply(seq_along(pop[ , "trace"]), function (i) cbind(pop[ , "trace"][[i]], k_o = i + 1)) #####c
  clust_list <- pop[ ,"clust"]
  n_filtered_list <- pop[ , "n_filtered"]

  # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
  clust <- do.call(rbind, clust_list)
  summary <- do.call(rbind, summary_list)
  trace <- do.call(rbind, trace_list)
  n_filtered <- do.call(rbind, n_filtered_list)
  n_filtered$k_o <- 2:(nrow(n_filtered) + 1)
  
  pop <- list(clust = clust, summary = summary, trace = trace, n_filtered = n_filtered)
  return(pop)
}

#------------------------------------------------------------------------------#
# Multi-objective optimsiation (kind of) ####
#------------------------------------------------------------------------------#
pop_update <- function(pop_old, pop) {
  # append pop
  pop <- Map(rbind, pop_old, pop)
  
  # index summary and trace for easy filtering
  # include Pareto front solution indices for easy filtering
  pop$summary$idx <- 1:nrow(pop$summary)
  pop$trace$idx <- NULL # remove column for the new
  
  if("k_o" %in% names(pop$trace)) {
    pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[ , c("run","batch","k_o","idx")]), by = dplyr::join_by(run, batch, k_o))
  } else {
    pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[ , c("run","batch","idx")]), by = dplyr::join_by(run, batch))
  }
  
  # find and filter exact duplicates
  if (!is.null(pop$clust)) {
    if (nrow(pop$clust) > 1) {
      dup <- find_dup(pop$clust, ari = 1)
      
      if (length(dup$idx) > 0) {
        # record the duplicate freq
        pop$summary$dup[as.numeric(names(dup$freq))] <- pop$summary$dup[as.numeric(names(dup$freq))] + as.vector(dup$freq)
        # filter the duplicates
        pop$clust <- pop$clust[-dup$idx, , drop = F]
        pop$summary <- pop$summary[-dup$idx, ]
        pop$trace <- subset(pop$trace, !idx %in% dup$idx)
        # append the dup count
        pop$n_filtered$dup[nrow(pop$n_filtered)] <- pop$n_filtered$dup[nrow(pop$n_filtered)] + length(dup$idx)
        
        # reindex the solution
        pop$summary$idx <- 1:nrow(pop$summary)
        pop$trace$idx <- NULL # remove column for the new
        
        if("k_o" %in% names(pop$trace)) {
          pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[ , c("run","batch","k_o","idx")]), by = dplyr::join_by(run, batch, k_o))
        } else {
          pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[ , c("run","batch","idx")]), by = dplyr::join_by(run, batch))
        }
      }
    }
  }
  
  return(pop)
}
#------------------------------------------------------------------------------#
eval_moo <- function(pop_pareto_objspace_list) {
  n_batch <- length(pop_pareto_objspace_list)
  # last pareto objspace to set upper and lower limits
  pareto_objspace <- pop_pareto_objspace_list[[n_batch]]

  # In case NULL for some batches
  for(i in 1:n_batch) {
    n_obj <- ncol(pop_pareto_objspace_list[[i]])
    obj <- names(pop_pareto_objspace_list[[i]])
    if (!is.null(n_obj)) break
  }

  # obtain upper and lower bounds from the total Pareto front
  ub <- vapply(1:n_obj, function(i) max(pareto_objspace[[obj[i]]]), numeric(1))
  lb <- vapply(1:n_obj, function(i) min(pareto_objspace[[obj[i]]]), numeric(1))
  
  # normalise the objectives
  pop_pareto_objspace_list <- lapply(pop_pareto_objspace_list, function(x) moocore::normalise(x, to_range = c(0,1), lower = lb, upper = ub))
  
  # compute IGD, IGD plus and HV
  igd <- vapply(1:(n_batch - 1),
                function(i)
                  moocore::igd(x = pop_pareto_objspace_list[[i]], reference = pop_pareto_objspace_list[[i+1]]),
                numeric(1))
  
  igd_plus <- vapply(1:(n_batch - 1),
                     function(i)
                       moocore::igd_plus(x = pop_pareto_objspace_list[[i]], reference = pop_pareto_objspace_list[[i+1]]),
                     numeric(1))
  
  hv <- vapply(1:n_batch,
               function(i)
                 moocore::hypervolume(x = pop_pareto_objspace_list[[i]], reference = rep(1, n_obj)),
               numeric(1))
  
  igd <- append(NA, igd)
  igd_plus <- append(NA, igd_plus)
  
  moo_quality <- data.frame(igd = igd, igd_plus = igd_plus, hv = hv)
  
  return(moo_quality)
}
#------------------------------------------------------------------------------#
parse_obj <- function(obj) {
  maximise_obj_idx <- grep(pattern = "^-", obj)
  obj <- sub("^-", "", obj)
  out <- list(obj = obj, maximise_obj_idx = maximise_obj_idx)
  return(out)
}

#------------------------------------------------------------------------------#
append_moo_out <- function(pop, obj, pareto_idx, moo_quality) {
  # create a column in summary to indicate Pareto front solutions
  pop$summary$pareto <- 0
  pop$summary$pareto[pareto_idx] <- 1
  pop$summary$pareto <- pop$summary$pareto
  
  # output the Pareto clust as clust 
  pop$clust <- pop$clust[pareto_idx, , drop = F]
  # output orignal idx for id in the summary table
  pop$clust_idx <- pareto_idx
  
  # output eval_moo() output
  pop$moo_quality <- moo_quality
  
  # output obj input parameter
  pop$obj <- obj
  
  return(pop)
}

#------------------------------------------------------------------------------#
filter_pareto_similar <- function(pop, ari) {
  
  # how not to rely on pareto_objspace but only pop$summary? so we can simplify the output
  # parse parameter obj to get the column indices and parsed obj names
  parse_obj_out <- parse_obj(pop$obj)
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  # use the parsed obj names in the following steps
  obj <- parse_obj_out$obj
  # select the columns of the objective space.
  pareto_objspace <- subset(pop$summary, pareto == 1)[ , obj]
  pareto_idx <- subset(pop$summary, pareto == 1)$idx
  clust <- pop$clust
  
  # multiply obj to be maximised by -1
  if (length(maximise_obj_idx) > 0) {
    for(j in maximise_obj_idx) {
      pareto_objspace[[ obj[j] ]] <- -pareto_objspace[[ obj[j] ]]
    }
  }
  
  # find_dup() is sensitive to order, therefore when it comes to similarity we want to prioritise based on the objectives
  # order the solutions by the objectives
  pareto_objspace$idx <- 1:nrow(pareto_objspace)
  ordered_idx <- eval(
    parse(
      text = paste0("pareto_objspace[order(", paste0(paste0("pareto_objspace$", obj), collapse = ","), "), ]$idx")
    ) 
  )
  
  ordered_clust <- clust[ordered_idx, ]
  similar <- find_dup(clust = ordered_clust, ari = ari)
  
  # reorder the idx according to the original idx
  # names(similar$freq) <- ordered_idx[as.numeric(names(similar$freq))]
  similar$idx <- ordered_idx[similar$idx]
  # similar$pairs_dup <- apply(similar$pairs_dup, c(1,2), function(x) ordered_idx[x])
  
  # filter similar if any
  if (length(similar$idx) > 0) {
    pop$clust <- clust[-similar$idx, ] # should i filter already or just update the summary?
    pop$clust_idx <- pop$clust_idx[-similar$idx]
    pop$summary$pareto_similar <- 0
    pop$summary$pareto_similar[pareto_idx[similar$idx]] <- 1 #############

    if (length(unique(pop$summary$k_o)) > 1) {
      counts <- subset(pop$summary[ , c("pareto_similar","k_o","batch")], pareto_similar == 1)
      if (length(counts$pareto_similar) == 0) {
        pop$n_filtered$pareto_similar <- 0
      } else {
        counts <- dplyr::group_by(counts, k_o, batch)
        counts <- as.data.frame(dplyr::summarise(counts, pareto_similar = sum(pareto_similar), .groups = "drop"))
        pop$n_filtered <- dplyr::left_join(pop$n_filtered, counts, by = dplyr::join_by(k_o, batch))
        pop$n_filtered$pareto_similar[is.na(pop$n_filtered$pareto_similar)] <- 0
      }
    } else {
      counts <- subset(pop$summary[,c("pareto_similar","batch")], pareto_similar == 1)
      if (length(counts$pareto_similar) == 0) {
        pop$n_filtered$pareto_similar <- 0
      } else {
        counts <- dplyr::group_by(counts, batch)
        counts <- as.data.frame(dplyr::summarise(counts, pareto_similar = sum(pareto_similar), .groups = "drop"))
        pop$n_filtered <- dplyr::left_join(pop$n_filtered, counts, by = dplyr::join_by(batch))
        pop$n_filtered$pareto_similar[is.na(pop$n_filtered$pareto_similar)] <- 0
      }
    }
    
    # rearrange the columns
    if ("k_o" %in% names(pop$n_filtered)) {
      pop$n_filtered <- pop$n_filtered[ , c("intersects", "size", "k1", "dup", "pareto_similar", "k_o", "batch")]
    } else {
      pop$n_filtered <- pop$n_filtered[ , c("intersects", "size", "k1", "dup", "pareto_similar", "batch")]
    }
  }
  
  return(pop)
}

#------------------------------------------------------------------------------#
blob_moo <- function (data, k = NULL, k_range = NULL, r_range = c(0.5,1), iter = 3, max_run = 500, run = 100, space_distmethod = "geodesic",
                      space_distmat = NULL, iter_method = "iterative", converge_ari = NULL,
                      filter_intersects = T, filter_clustsize = T, max_na = 0.05,
                      obj = NULL, pareto_similar_ari = NULL, ...) {
  
  #-------------------------------------------------------------------#
  # sort the batches and runs from the input ####
  if (max_run < run) stop("max_run < run")
  
  # number of standard batches as set for the parameter
  n_batch <- floor(max_run/run)
  # remaining run to match the max number of runs in total
  remain_run <- max_run - n_batch*run 
  # a vector of runs per batch 
  RUN <- c(rep(run, n_batch), remain_run)
  # remove batch of 0 run
  RUN <- RUN[RUN != 0]
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # initialise MOO required object ####
  # pareto indices for every batch
  pop_pareto_idx_list <- list()
  # pareto objspace for every batch
  pop_pareto_objspace_list <- list()
  
  if(is.null(obj)) {
    if (!is.null(k)) { # if blob_populate() is called
      if (filter_clustsize == F) {
        obj <- c("space_wcss","-time_range_mean","-time_evenness_mean","size_diff")
      } else {
        obj <- c("-k","n_removed","space_wcss","-time_range_mean","-time_evenness_mean","size_diff")
      }
    } else { # if blob_kpopulate() is called
      obj <- c("-k","n_removed","space_wcss","-time_range_mean","-time_evenness_mean","size_diff")
    }
  }
  
  # parse parameter obj to get the column indices and parsed obj names
  parse_obj_out <- parse_obj(obj)
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  # store parameter obj input for output
  obj_input <- obj
  # use the parsed obj names in the following steps
  obj <- parse_obj_out$obj
  
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # main loop to compare appending batches ####
  for (i in 1:length(RUN)) {
    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    # populate a bunch of solutions
    if (!is.null(k)) {
      pop <- blob_populate(data = data, k = k, r_range = r_range, iter = iter, run = RUN[i], space_distmethod = space_distmethod,
                           space_distmat = space_distmat, iter_method = iter_method, converge_ari = converge_ari,
                           filter_intersects = filter_intersects, filter_clustsize = filter_clustsize, max_na = max_na, ...)
    } else {
      pop <- blob_kpopulate(data = data, k_range = k_range, r_range = r_range, iter = iter, run = RUN[i], space_distmethod = space_distmethod,
                            space_distmat = space_distmat, iter_method = iter_method, converge_ari = converge_ari,
                            filter_intersects = filter_intersects, filter_clustsize = filter_clustsize, max_na = max_na, ...)
    }
    
    
    # index the batch
    pop$summary$batch <- i
    pop$trace$batch <- i
    pop$n_filtered$batch <- i
    
    # if there is at least a single solution, index the summary and trace tables
    if(!is.null(pop$clust)) {
      # index summary and trace for easy filtering
      # include Pareto front solution indices for easy filtering
      pop$summary$idx <- 1:nrow(pop$summary)
      if("k_o" %in% names(pop$trace)) {
        pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[,c("run", "k_o", "batch", "idx")]), by = dplyr::join_by(run, k_o, batch))
      } else {
        pop$trace <- dplyr::left_join(pop$trace, dplyr::distinct(pop$summary[,c("run", "batch", "idx")]), by = dplyr::join_by(run, batch)) 
      }
    } 
    
    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    # update pop as number of solutions increase from batch to batch ####
    if (i > 1) {
      pop <- pop_update(pop_old = pop_old, pop = pop)
    }
    
    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    # if there are more than one solution, find the parteo optimal solution(s)
    if (!is.null(pop$clust) & nrow(pop$clust) > 1) {
      # select Pareto front solutions ####
      # objective space
      pop_objspace <- pop$summary[obj]
      
      # multiply obj to be maximised by -1
      if (length(maximise_obj_idx) > 0) {
        for(j in maximise_obj_idx) {
          pop_objspace[[ obj[j] ]] <- -pop_objspace[[ obj[j] ]]
        }
      }
      
      # extract the idx of Pareto optimal solutions
      pop_pareto_idx_list[[i]] <- which(moocore::is_nondominated(pop_objspace) == TRUE)
      # filter for the Pareto optimal objective space
      pop_pareto_objspace_list[[i]] <- pop_objspace[pop_pareto_idx_list[[i]], obj]
    }
    
    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    # store the current as old to be appended in the next loop
    pop_old <- pop
  }
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # if there there is no solution at all, return NULL
  if (is.null(pop_pareto_idx_list[[n_batch]])) return(NULL)
  
  # Pareto optimal idx and objspace from the last appended pop
  pareto_idx <- pop_pareto_idx_list[[n_batch]]
  
  # evaluate MOO performance
  moo_quality <- eval_moo(pop_pareto_objspace_list)
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # append MOO output ####
  pop <- append_moo_out(pop = pop, obj = obj_input, pareto_idx = pareto_idx, moo_quality = moo_quality)
  
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  # filter similar solutions on the Pareto front ####
  if (!is.null(pareto_similar_ari)) {
    print(pareto_idx)
    if (length(pareto_idx) > 1) {
      pop <- filter_pareto_similar(pop, ari = pareto_similar_ari)
    }
  }
  #-------------------------------------------------------------------#
  #-------------------------------------------------------------------#
  return(pop)
}

#------------------------------------------------------------------------------#
# Simulate data ####
#------------------------------------------------------------------------------#
# Generate coordinates on a circle line
gen_circle_coords <- function(radius, n) {
  # radius: a numeric of the radius
  # n: an integer of the number of vertices
  
  # Generate angles for each vertex
  angles <- seq(0, 2 * pi, length.out = n + 1)[-1]  # Exclude the last angle (2*pi)
  
  # Calculate x and y coordinates
  x_coords <- radius * cos(angles)
  y_coords <- radius * sin(angles)
  
  # Combine coordinates into a data frame
  coords <- as.matrix(data.frame(x = x_coords, y = y_coords))
  
  return(coords)
}
#------------------------------------------------------------------------------#
# Generate clusters of Guassian data
gen_gaussian_data <- function(n, size, center = NULL, ...) {
  # k: an integer indicating the number of clusters
  # size: an integer indicating the size of the cluster. To generate unequal sized clusters, supply a vector of integers
  # center: a data frame of coordinates. If none is given, it will be sampled randomly
  
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
    # Label cluster
    data[[i]] <- cbind(data[[i]],i)
  }
  
  data <- do.call(rbind, data)
  colnames(data) <- c("x","y","clust")
  data <- as.data.frame(data)
  data$clust <- as.factor(data$clust)
  return(data)
}
