#------------------------------------------------------------------------------#
# Functions for the local search algorithm ####
#------------------------------------------------------------------------------#
#' Compute distance matrix
#' 
#' @description
#' This function computes and returns a distance matrix.
#'
#' @param data a numeric matrix or data frame of spatial coordinates.
#' @param method a string of the method. Either "geodesic" or "euclidean" is available. Default is "geodesic".
#'
#' @details
#' km is the unit when "geodesic" is specified.
#'
#' @returns a numeric distance matrix.
#'
#' @export

compute_distmat <- function(data, method = "geodesic") {
  # checks
  method <- match.arg(method, choices = c("geodesic", "euclidean"))
  stopifnot("data must be a numeric matrix or data frame of spatial coordinates." = is.matrix(data) | is.data.frame(data))
  # convert data to matrix
  data <- as.matrix(data)
  # checks
  stopifnot("data must be a numeric matrix or data frame of spatial coordinates." = is.numeric(data))
  #-------------------------------------------------------------------#
  if (method == "geodesic") {
    distmat <- as.matrix(geosphere::distm(data, fun = geosphere::distGeo)) / 1000
  }
  if (method == "euclidean") {
    distmat <- as.matrix(stats::dist(data, method = "euclidean"))
  }
  return(distmat)
}

#------------------------------------------------------------------------------#
#' Compute spatial costs
#'
#' @description
#' This function computes and returns the spatial costs of all points to the medoid of a given cluster.
#'
#' @param space_distmat a spatial distance matrix.
#' @param clust_points a numeric vector of point indices in the cluster.
#' @param weights a numeric vector of weights for each data point. Default is NULL.
#' 
#' @returns a numeric vector of spatial costs of a given cluster.

compute_spacecost <- function(space_distmat, clust_points, weights = NULL) {
  # checks
  stopifnot("weights must be a numeric vector of the same length as the number of data points" = (is.vector(weights) & is.numeric(weights) & (dim(space_distmat)[1]==length(weights))) | is.null(weights))
  #-------------------------------------------------------------------#
  w <- weights
  if (is.null(weights)) {
    w <- rep(1,nrow(space_distmat))
  }
  
  W <- outer(w, w, "+")/2
  space_distmat <- W * space_distmat
  
  clust_mat <- space_distmat[clust_points, clust_points, drop = FALSE]
  distsum <- rowSums(clust_mat)
  medoid <- clust_points[which.min(distsum)[1]]
  
  space_cost <- (space_distmat[ , medoid])
  return(space_cost)
}

#------------------------------------------------------------------------------#
#' Check for spatial intersects
#' 
#' @description
#' This function checks and returns a Boolean if any clusters intersect in space.
#'
#' @param data a numeric matrix or data frame of spatial coordinates.
#' @param clust a numeric vector of the cluster assignment. Default is NULL.
#' @param coords a vector of strings or numeric values indicating the columns of spatial coordinates (e.g. "longitude" and "latitude"). Default is the first two columns.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()]. Default is 4326.
#' @param hull_convex_ratio a numeric value indicating the convexity of the hulls passed onto [sf::st_concave_hull()]. 1 returns convex and 0 maximally concave hulls. Default is 0.5.
#'
#' @returns TRUE or FALSE

check_intersects <- function(data, clust = NULL, coords = c(1,2), crs = 4326, hull_convex_ratio = 0.5) {
  # checks
  stopifnot("data must be a numeric matrix or data frame of spatial coordinates." = (is.matrix(data) | is.data.frame(data)) & is.numeric(as.matrix(data)))
  #-------------------------------------------------------------------#
  # st_union breaks when s2 is on
  suppressMessages(sf::sf_use_s2(FALSE))
  data <- as.data.frame(data)
  #-------------------------------------------------------------------#
  if (!is.null(clust)) { data$clust <- clust }
  #-------------------------------------------------------------------#
  data_sf <- sf::st_as_sf(data, coords = coords, crs = crs)
  #-------------------------------------------------------------------#
  # obtain convex/concave hulls
  suppressMessages({
    hulls <- stats::aggregate(data_sf$geometry, by = list(clust = data_sf$clust), function(x){
      x <- sf::st_combine(x)
      x <- sf::st_union(x, by_feature = TRUE)
      # If by_feature is TRUE each feature geometry is unioned individually. This can for instance be used to resolve internal boundaries after polygons were combined using st_combine. https://r-spatial.github.io/sf/reference/geos_combine.html
      return(x)
    }) 
    hulls <- sf::st_as_sf(hulls)
    hulls <- sf::st_concave_hull(hulls, ratio = hull_convex_ratio) # convex hulls
    #-------------------------------------------------------------------#
    # mark TRUE/FALSE for each comparison
    hulls <- sf::st_make_valid(hulls)
    intersects <- sf::st_intersects(hulls$geometry, sparse = F)
  })
  diag(intersects) <- NA
  bool <- if (any(intersects == TRUE, na.rm = TRUE)) TRUE else FALSE
  #-------------------------------------------------------------------#
  # reset s2
  suppressMessages(sf::sf_use_s2(TRUE))
  return(bool)
}

#------------------------------------------------------------------------------#
#' Evaluate blobs
#' 
#' @description
#' This function evaluates the performance of clustering.
#'
#' @inheritParams check_intersects
#' @inheritParams compute_spacecost
#' @param data a data frame with spatial coordinates, age and cluster of the data 
#' @param age a string or numeric value indicating the column of age of the data. Default is the third column.
#' @param space_distmat a spatial distance matrix. Default is NULL.
#' @param space_distmethod a string of the method, used when `space_distmat` is not specified. Either "geodesic" or "euclidean" is available. Default is "geodesic".
#' 
#' @details
#' The size threshold of a cluster is defined as \eqn{\frac{n}{2k}} where \eqn{n} is the number of data points and \eqn{k} is the number of clusters.
#'
#' @returns a data frame of summary statistics.

eval_blobs <- function(data,
                       coords = c(1,2),
                       age = 3,
                       crs = 4326,
                       hull_convex_ratio = 0.5,
                       space_distmat = NULL,
                       space_distmethod = NULL,
                       weights = NULL) {
  
  # checks
  stopifnot("data must be a data frame with spatial coordinates, age and cluster." = is.data.frame(data) & !is.null(data[ , coords]) & !is.null(data[ , age]) & !is.null(data$clust))
  #-------------------------------------------------------------------#
  # total number of points
  N <- nrow(data)
  # total number of clusters
  k <- length(unique(stats::na.omit(data$clust))) # NA is excluded
  # initialise empty vectors 
  space_d <- time_r <- time_e <- n <- numeric(k)
  #-------------------------------------------------------------------#
  # check if is.null(space_distmat)
  space_distmat <- check_space_distmat(space_distmat = space_distmat, space_distmethod = space_distmethod)
  #-------------------------------------------------------------------#
  # loop over K to obtain within cluster statistics
  for (j in 1:k) {
    clust_points <- which(data$clust == j)
    if (length(clust_points) == 0) next
    data_k <- subset(data, data$clust == j)
    #-------------------------------------------------------------------#
    # spatial objective
    space_d[j] <- sum(compute_spacecost(space_distmat = space_distmat,
                                        clust_points = clust_points,
                                        weights = weights)[clust_points])
    #-------------------------------------------------------------------#
    # temporal objectives
    age <- data_k[ , age] 
    time_r[j] <- max(age, na.rm = T) - min(age, na.rm = T)
    time_e[j] <- 1 / (1 + stats::var(diff(sort(age)))) # NA if there are fewer than 3 data points
    #-------------------------------------------------------------------#
    # clust size
    n[j] <- length(clust_points)
  }
  #-------------------------------------------------------------------#
  # calculate the summary statistics
  space_wcd <- sum(space_d)
  time_wcr <- mean(time_r, na.rm = TRUE)
  time_wce <- mean(time_e, na.rm = TRUE)
  # na.rm = TRUE to assess the overall performance as we do not care clusters with fewer than 3 points
  #-------------------------------------------------------------------#
  # flag the number of clusters below the threshold
  clust_sizef <- length(which(n < N/k/2))
  #-------------------------------------------------------------------#
  # evaluate if blobs are intersecting in space
  intersects <- check_intersects(data = data, coords = coords, crs = crs, hull_convex_ratio = hull_convex_ratio)
  #-------------------------------------------------------------------#
  # return a data frame of all the statistics
  summary <- data.frame(k = k,
                        space_wcd = space_wcd,
                        time_wcr = time_wcr,
                        time_wce = time_wce,
                        intersects = intersects,
                        clust_sizef = clust_sizef)
  
  return(summary)
}

#------------------------------------------------------------------------------#
#' Assign the starting cluster members
#' 
#' @description
#' This function assigns the initial cluster points by approximating the maximum spatial separation between points or by random assignment.
#'
#' @param data a data matrix or data frame with spatial coordinates and age of the data.
#' @param k an integer of the number of clusters.
#' @param random_init a Boolean to randomly select initial cluster points. Default is FALSE.
#' @inheritParams compute_spacecost
#'
#' @returns a data frame with spatial coordinates, age, cluster and order of the data.

init_blobs <- function(data, k, space_distmat, random_init = FALSE) {
  if (random_init == TRUE) {
    data <- as.data.frame(data)
    data$clust <- NA
    init <- sample(1:nrow(data), k)
    data$clust[init] <- 1:k
    data$order <- 1:nrow(data)
    return(data)
  }
  #-------------------------------------------------------------------#
  data <- as.data.frame(data)
  # assign an index to put data back into correct order at the end
  data$order <- 1:nrow(data)
  # initialise one point each to k blobs
  data$clust <- NA
  #-------------------------------------------------------------------#
  # a bit faster to handle data as a matrix
  data <- as.matrix(data)
  #-------------------------------------------------------------------#
  # start from k roughly equally spaced random locations 
  # (just permute a bunch and pick the one with the least smallest distance)
  N <- 100
  mat <- matrix( , N, k)
  min_dist <- numeric(N)
  #-------------------------------------------------------------------#
  for(n in 1:N){
    i <- sort(sample(1:nrow(data), size = k)) # sample k points from the data, sort them 
    mat[n, ] <- i # store sorted index in the matrix by row
    cb <- utils::combn(k, 2) # all combinations of k chooses 2 by column
    NC <- ncol(cb) # NC number of combinations
    dist <- numeric(NC)
    for(c in 1:NC) dist[c] <- space_distmat[ i[cb[1, c]], i[cb[2, c]] ] # extract from distance matrix the distances for all combinations of points
    min_dist[n] <- min(dist)
  }
  #-------------------------------------------------------------------#
  # pick the set with largest minimum distance between two points to ensure maximum separation between k points;
  # pick the first one if two are tied 
  init <- mat[which(min_dist == max(min_dist))[1], ]
  # Assign cluster memberships to the starting points
  data[init, "clust"] <- 1:k
  data <- as.data.frame(data)
  #-------------------------------------------------------------------#
  return(data)
}

#------------------------------------------------------------------------------#
#' Assign clusters
#' 
#' @description
#' This function assigns clusters for a given k and r.
#' 
#' @param data a data frame output from [init_blobs()].
#' @param k an integer of the number of clusters.
#' @param r an numeric value of the spatial relative weight of range \eqn{[0,1]}.
#' @param age a string or numeric value indicating the column of age. Default is the third column.
#' @inheritParams compute_spacecost
#'
#' @returns a data frame with spatial coordinates, age, cluster and order of the data.

find_blobs <- function(data, k, r, space_distmat, age = 3, weights = NULL) {
  # a bit faster to handle data as a matrix
  data <- as.matrix(data)
  #-------------------------------------------------------------------#
  # for indexing speed, keep those assigned at the top
  data <- data[order(data[ , "clust"]), ]
  #-------------------------------------------------------------------#
  # randomise the order of the points to be assigned (tba)
  tba_points <- which(is.na(data[ , "clust"])) # to be assigned
  # for when iteration is implemented there will be no NA
  tba_points <- if (length(tba_points) == 0) 1:nrow(data) else tba_points
  tba_points <- sample(tba_points)
  a_points <- which(!1:nrow(data) %in% tba_points) # assigned
  # assigned at the top, randomised unassigned rows that follow
  data <- data[c(a_points, tba_points), ] 
  #-------------------------------------------------------------------#
  # Reorder distmat to match the reordered data
  order <- data[ , "order"]
  space_distmat <- space_distmat[order, order]
  #-------------------------------------------------------------------#
  # Extract clust to make the code cleaner
  clust <- data[ , "clust"]
  #-------------------------------------------------------------------#
  # initialise stat matrices for the unassigned points
  start <- if (length(a_points) > nrow(data)) 1 else length(a_points) + 1
  N <- length(start:nrow(data))
  space_cost <- time_cost <- n <- numeric(k)
  #-------------------------------------------------------------------#
  # precompute spacecost
  space_costmat <- vapply(1:k, function(j) {
    clust_points <- which(clust == j)
    compute_spacecost(space_distmat = space_distmat,
                      clust_points = clust_points,
                      method = method,
                      weights = weights)
  }, FUN.VALUE = numeric(nrow(data)))
  
  # prefetch age column for time computation
  age <- data[ , age]
  #-------------------------------------------------------------------#
  # loop through every point (incremental updating)
  for (i in start:nrow(data)) {
    #-------------------------------------------------------------------#
    # loop through every k
    for (j in 1:k) {
      # O(c^2) is unavoidable 
      clust_points <- which(clust == j)
      n[j] <- length(clust_points)
      if (n[j] == 0) next
      #-------------------------------------------------------------------#
      # compute the spatial cost == j
      space_cost[j] <- space_costmat[i, j]
      #-------------------------------------------------------------------#
      # compute the temporal cost
      clust_points_tmp <- if (i %in% clust_points) clust_points[clust_points != i] else clust_points
      #-------------------------------------------------------------------#
      if (length(clust_points_tmp) == 0) {
        time_cost[j] <- 0  # only i was in the cluster
      } else {
        time_cost[j] <- min(abs(age[i] - age[clust_points_tmp]))
      }
    }
    #-------------------------------------------------------------------#
    # normalise the cost
    space_cost_norm <- (space_cost - min(space_cost, na.rm = T)) / (max(space_cost, na.rm = T) - min(space_cost, na.rm = T))
    time_cost_norm <- (time_cost - min(time_cost, na.rm = T)) / (max(time_cost, na.rm = T) - min(time_cost, na.rm = T))
    space_cost_norm[is.na(space_cost_norm)] <- 0
    time_cost_norm[is.na(time_cost_norm)] <- 0
    #-------------------------------------------------------------------#
    # weighted scalarising (into a single cost)
    cost <- space_cost_norm*r - time_cost_norm*(1-r)
    #-------------------------------------------------------------------#
    # assign cluster
    # check if there are tied clusters
    clust_tmp <- which(cost == min(cost, na.rm = T))
    #-------------------------------------------------------------------#
    if (length(clust_tmp) > 1) {
      # check if one of those has fewer points
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
  #-------------------------------------------------------------------#
  # Put the data back to order
  data[ , "clust"] <- clust
  data <- data[order(data[ , "order"]), ]
  data <- as.data.frame(data)
  #-------------------------------------------------------------------#
  return(data)
}

#------------------------------------------------------------------------------#
#' Core local-search algorithm #### come back here
#' 
#' @description
#' This function performs an iterative bi-objective local search algorithm to assign clusters for a given k and r.
#' 
#' @inheritParams init_blobs
#' @inheritParams find_blobs
#' @inheritParams eval_blobs
#' @param iter an integer of the number of iterations. Default is 10.
#' @param converge_ari a numeric value of the Adjusted Rand Index (ARI) that sets the convergence threshold between iterations. It must be \eqn{[0,1]}. Default is 1.
#' @param filter_intersects a Boolean to remove an assignment with intersects in space. Default is TRUE.
#' @param filter_clustsize a Boolean to remove an assignment with clusters below the expected size. Default is TRUE.
#' @param max_na a numeric value of the maximum proportion of NAs allowed. It must be \eqn{[0,1]}. Default is 0.05.
#'
#' @details
#' Local search completes when consecutive iterations reach ARI specified by `converge_ari`, or else 1 (identical) , after at least 2 iterations.
#'
#' The size threshold of a cluster is defined as \eqn{\frac{n}{2k}} where \eqn{n} is the number of data points and \eqn{k} is the number of clusters.
#'
#' @returns 
#' a numeric value of 1 when the output is invalid,
#' otherwise a data frame with the following elements
#' \itemize{
#'   \item \code{data}: a data frame of the input data with assigned clusters as a column.
#'   \item \code{summary}: a data frame of the summary statistics. 
#'   \item \code{trace}: a data frame of the summary statistics per iteration.
#' }
#'
#' @seealso [sf::st_as_sf()]
#'
#' @export

blob_search <- function(data,
                        k,
                        r,
                        iter = 10,
                        converge_ari = 1,
                        coords = c(1,2),
                        age = 3,
                        crs = 4326,
                        hull_convex_ratio = 0.5,
                        random_init = FALSE,
                        filter_intersects = TRUE,
                        filter_clustsize = TRUE,
                        space_distmat = NULL,
                        space_distmethod = NULL,
                        weights = NULL) {
  
  # checks
  
  #-------------------------------------------------------------------#
  # select the relevant columns
  data <- data[, c(coords,age)]
  #-------------------------------------------------------------------#
  space_distmat <- check_space_distmat(space_distmat, space_distmethod)
  #-------------------------------------------------------------------#
  # search algorithm
  # init_blobs() to pick medoids
  data <- init_blobs(data = data, k = k, space_distmat = space_distmat, random_init = random_init)
  #-------------------------------------------------------------------#
  # initialise counter counting find_blobs()
  t <- 0
  # intialise trace table
  trace <- data.frame()
  #-------------------------------------------------------------------#
  for (i in 1:iter) {
    data_old <- data
    #-------------------------------------------------------------------#
    # find_blobs()
    data <- find_blobs(data = data, k = k, r = r, space_distmat = space_distmat, weights = weights)
    #-------------------------------------------------------------------#
    # count find_blobs() executed
    t <- t + 1
    #-------------------------------------------------------------------#
    if (t > 0) {
      #-------------------------------------------------------------------#
      # check convergence
      ari <- mclust::adjustedRandIndex(data$clust, data_old$clust)
      #-------------------------------------------------------------------#
      # eval_blobs()
      eval_out <- eval_blobs(data,
                             crs = crs,
                             hull_convex_ratio = hull_convex_ratio,
                             space_distmat = space_distmat,
                             weights = weights)
      clust_sizef <- eval_out$clust_sizef
      trace_row <- eval_out$summary
      trace_row$iter <- t
      trace_row$ari <- ari
      trace <- rbind(trace, trace_row)
      #-------------------------------------------------------------------#
      # if converged between t and t-1, break
      if (!is.null(converge_ari)) {
        if (all(ari >= converge_ari) == TRUE & t >= 3) break
      }
    }
  }
  #-------------------------------------------------------------------#
  trace$r <- r
  summary <- trace[nrow(trace), ]
  # order is not useful for output, remove the column
  data$order <- NULL
  #-------------------------------------------------------------------#
  # filters / constraints
  N <- nrow(data)
  #-------------------------------------------------------------------#
  # 1. intersects
  if (filter_intersects == T) {
    intersects <- summary$intersects
    if (intersects == T) {
      return(1)
    }
  }
  #-------------------------------------------------------------------#
  # 2. cluster size
  # initialise n_removed column. If no point is removed, the entries will all be 0
  summary$n_removed <- 0
  if (filter_clustsize == T) {
    if (length(clust_below_size) > 0) {
      #-------------------------------------------------------------------#
      # assign NA to removed clusters
      n_removed <- summary$n_fail
      data$clust[which(data$clust %in% clust_below_size)] <- NA
      #-------------------------------------------------------------------#
      # reassign the cluster number
      # make sure there is no gap in the sequence of k
      data$clust <- reorder_clust(data$clust)
      #-------------------------------------------------------------------#
      # update eval_out$summary
      eval_out_updated <- eval_blobs(data,
                                     crs = crs,
                                     hull_convex_ratio = hull_convex_ratio,
                                     space_distmat = space_distmat,
                                     weights = weights)
      updated_cols <- intersect(names(summary), names(eval_out_updated$summary))
      summary[ , updated_cols] <- eval_out_updated$summary[ , updated_cols]
      summary$n_removed <- n_removed
    }
  } else {
    #-------------------------------------------------------------------#
    # reassign the cluster number
    # make sure there is no gap in the sequence of k
    data$clust <- reorder_clust(data$clust)
  }
  #-------------------------------------------------------------------#
  # 3. filter too many NAs and k < 2
  # return NULL if too many points are removed
  if (summary$n_removed > N * max_na) return(2)
  if (summary$k < 2) return(3)
  #-------------------------------------------------------------------#
  # 4. remove some columns
  summary$n_fail <- NULL # confusing to have both n_removed and this in the output
  trace$n_fail <- NULL # confusing to have both n_removed and this in the output
  #-------------------------------------------------------------------#
  # 5. remove clust.below.size as filter has been applied so it is irrelevant to the downstream
  clust_below_size <- NULL
  #-------------------------------------------------------------------#
  summary <- summary[, c("k", "r", "space_wcd", "time_wcr", "time_wce", "intersects", "n_removed", "iter", "ari")]
  trace <- trace[, c("k", "r", "space_wcd", "time_wcr", "time_wce", "iter", "ari")]
  rownames(summary) <- NULL
  rownames(trace) <- NULL
  #-------------------------------------------------------------------#
  blob <- list(data = data,
               summary = summary,
               trace = trace,
               clust_below_size = clust_below_size,
               space_kmat_optim_out = space_kmat_optim_out)
  class(blob) <- "blob"
  return(blob)
}

#------------------------------------------------------------------------------#
# Populate solutions ####
#------------------------------------------------------------------------------#
#' Find duplicates
#' 
#' @description
#' This function finds duplicates in a set of cluster assignments. It uses [mclust::adjustedRandIndex()] to measure similarity.
#' 
#' @param clust a numeric matrix of cluster assignments. Each row is an assignment.
#' @param ari a numeric value of the Adjusted Rand Index (ARI) that sets duplication threshold between two assignments. It must be \eqn{[0,1]}. Default is 1. See also [mclust::adjustedRandIndex()].
#' 
#' @returns a list of the following objects.
#' \itemize{
#'   \item \code{idx}: a numeric vector of duplicate indices.
#'   \item \code{pairwise_ari}: a data frame of summary statistics.
#'   \item \code{freq}: a contingency table of the frequency of duplicates, an object of class "table", an array of integer values. See also [base::table()].
#'   \item \code{pairs_dup}: a numeric matrix of pairs of duplicates. The second row stores indices of duplicates of the first row.
#' }
#'
#' @seealso [base::table()], [mclust::adjustedRandIndex()]
#'
#' @export

find_dup <- function (clust, ari = 1) {
  # as NA is ignored by mclust::adjustedRandIndex()
  clust[is.na(clust)] <- 0
  #-------------------------------------------------------------------#
  # number of solutions
  n_sol <- nrow(clust)
  # get all combination of pairs
  pairs <- utils::combn(n_sol, 2)
  # calculate the pairwise ARI for all pairs of solutions
  pairwise_ari <- future.apply::future_vapply(1:ncol(pairs),
                                              function(x) mclust::adjustedRandIndex(clust[pairs[1,x], ], clust[pairs[2,x], ]),
                                              numeric(1),
                                              future.seed = TRUE)
  #-------------------------------------------------------------------#
  # index the pairs with ari >= ari
  pairs_dup_idx <- which(pairwise_ari >= ari)
  # subset the columns of duplicated pairs
  pairs_dup <- pairs[, pairs_dup_idx, drop = FALSE]
  #-------------------------------------------------------------------#
  # if a ~ b and b ~ c, then not necessarily a ~ c
  # the following steps makes sure they are the pairs with ari >= ari.dup 
  if (ncol(pairs_dup) > 1) {
    #-------------------------------------------------------------------#
    # dependence boolean
    pairs_dup_depend_bool <- logical()
    pairs_dup_tmp <- pairs_dup
    for (i in 2:ncol(pairs_dup)) {
      # if the reference is in one of the previous duplicates
      if(pairs_dup_tmp[1, i] %in% pairs_dup_tmp[2, 1:(i-1)]) {
        # mark the dependence
        pairs_dup_depend_bool[i] <- TRUE
        # break the chain by assigning 0 to the duplicate
        pairs_dup_tmp[2,i] <- 0
      } else {
        pairs_dup_depend_bool[i] <- FALSE
      }
    }
    pairs_dup_depend_idx <- which(pairs_dup_depend_bool == TRUE)
    #-------------------------------------------------------------------#
    # remove the dependent columns 
    if (length(pairs_dup_depend_idx) > 0) {
      pairs_dup <- pairs_dup[ , -pairs_dup_depend_idx, drop = FALSE]
      pairwise_ari <- pairwise_ari[pairs_dup_idx[-pairs_dup_depend_idx]]
    }
  }
  #-------------------------------------------------------------------#
  # note it is sensitive to order
  # second row is considered the duplicate of the first row
  idx <- unique(pairs_dup[2, ])
  #-------------------------------------------------------------------#
  # frequency of duplicates
  freq <- table(pairs_dup[1,])
  #-------------------------------------------------------------------#
  dup_out <- list(idx = idx,
                  pairwise_ari = pairwise_ari,
                  freq = freq,
                  pairs_dup = pairs_dup)
  return(dup_out)
}
#------------------------------------------------------------------------------#
#' Convert a list of blob objects to a pop object
#' 
#' @description
#' This function converts a list of `blob` objects to a `pop` object
#' and the counts of filtered solutions.
#' 
#' @param blob_list a list of output from [blob_search()]
#' 
#' @returns a `pop` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of the cluster assignments. Each row is a solution.
#'   \item \code{summary}: a data frame of the summary statistics.
#'   \item \code{trace}: a data frame of the summary statistics per iteration.
#'   \item \code{n_filtered}: a data frame of the numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: a list of [stats::optim()] output from the optimisation of \eqn{\beta} in [distmat_to_kmat()].
#' }

convert_to_pop <- function(blob_list) {
  pop <- blob_list
  #-------------------------------------------------------------------#
  # count the runs that are filtered out, informative about the range of r to be searched
  filtered_intersects <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 1)
  filtered_intersects <- sapply(filtered_intersects, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_intersects <- sum(filtered_intersects)
  #-------------------------------------------------------------------#
  filtered_size <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 2)
  filtered_size <- sapply(filtered_size, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_size <- sum(filtered_size)
  #-------------------------------------------------------------------#
  filtered_k1 <- lapply(pop, function(x) if (is.list(x)) FALSE else x == 3)
  filtered_k1 <- sapply(filtered_k1, function(x) if (length(x) != 1) x <- FALSE else x) 
  filtered_k1 <- sum(filtered_k1)
  #-------------------------------------------------------------------#
  # if blob returns 1 or 2, set NULL so rbind in the following step will omit them
  pop <- lapply(pop, function(x) if(!is.list(x)) NULL else x)
  #-------------------------------------------------------------------#
  # extract each element and add a column to indicate the run
  pop <- do.call(rbind, pop) 
  data_list <- lapply(seq_along(pop[, "data"]), function (i) cbind(pop[, "data"][[i]], run = i))
  summary_list <- lapply(seq_along(pop[, "summary"]), function (i) {
    if (!is.null(pop[, "summary"][[i]])) { cbind(pop[, "summary"][[i]], run = i) }
  })
  trace_list <- lapply(seq_along(pop[, "trace"]), function (i) {
    if (!is.null(pop[, "trace"][[i]])) { cbind(pop[, "trace"][[i]], run = i) }
  })
  #-------------------------------------------------------------------#
  # extract clust
  clust_list <- lapply(seq_along(data_list), function (i) data_list[[i]]$clust)
  #-------------------------------------------------------------------#
  # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
  clust <- do.call(rbind, clust_list)
  summary <- do.call(rbind, summary_list)
  trace <- do.call(rbind, trace_list)
  #-------------------------------------------------------------------#
  # find and filter exact duplicates
  filtered_dup <- 0
  if (!is.null(clust)) {
    summary$dup <- 0
    # here, clust <- do.call(rbind, clust_list) must return a matrix even if there is only one solution
    if (nrow(clust) > 1) {
      dup <- find_dup(clust, ari = 1)
      if (length(dup$idx) > 0) {
        #-------------------------------------------------------------------#
        # record the duplicate freq
        summary$dup[as.numeric(names(dup$freq))] <- as.vector(dup$freq)
        #-------------------------------------------------------------------#
        # filter the duplicates
        clust <- clust[-dup$idx, , drop = F]
        summary <- summary[-dup$idx, ]
        trace <- trace[-which(trace$run %in% dup$idx),]
        #-------------------------------------------------------------------#
        # record the total no. of dup filtered
        filtered_dup <- length(dup$idx)
      }
    }
    #-------------------------------------------------------------------#
    # index the summary and trace
    summary$idx <- 1:nrow(summary)
    trace <- merge(trace, summary[, c("idx","run")], by = "run", all.x = TRUE)
  }
  #-------------------------------------------------------------------#
  n_filtered <- data.frame(intersects = filtered_intersects,
                           size = filtered_size,
                           k1 = filtered_k1,
                           dup = filtered_dup)
  #-------------------------------------------------------------------#
  summary <- summary[, c("idx", "k", "r", "run", "space_wcd", "time_wcr", "time_wce", "intersects", "n_removed", "iter", "ari", "dup")]
  trace <- trace[, c("idx", "k", "r", "run", "space_wcd", "time_wcr", "time_wce", "iter", "ari")]
  
  rownames(summary) <- NULL
  rownames(trace) <- NULL
  #-------------------------------------------------------------------#
  pop <- list(clust = clust,
              summary = summary,
              trace = trace,
              n_filtered = n_filtered)
  return(pop)
}

#------------------------------------------------------------------------------#
#' Populate solutions by weighted sum scalarisation
#' 
#' @description
#' This function populates solutions by weighted sum scalarisation of the bi-objective
#' local search algorithm in [blob_search()] for a given k. 
#' 
#' @inheritParams start_blobs
#' @inheritParams find_blobs
#' @inheritParams blob_search
#' @param k an integer value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the number of clusters.
#' @param r a numeric value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the relative spatial weight.
#' They must be \eqn{[0,1]}. Default is c(0.5,1).
#' @param run an integer of the number of runs. Default is 100.
#' 
#' @details
#' When "kkmeans" method is specified, a diffusion kernel is applied to compute the distance to the centroid in space for
#' kernel k-means clustering in space. 
#' See [distmat_to_kmat()] for more details for converting a distance matrix to a kernel matrix.
#'
#' Clusters are assigned in every iteration. It iterates until the set length or convergence. 
#'
#' When `converge_ari` is specified, convergence between iterations is defined.
#' The search will complete when ARI meets the threshold and least 3 iterations are run.
#'
#' The critical size of a cluster is defined as \eqn{\frac{n}{2k}} where \eqn{n} is the number of data point and \eqn{k} is the number of clusters.
#' 
#' Scalarisation is achieved by varying the relative spatial weight generated by Latin hypercube sampling using [lhs::randomLHS()].
#' 
#' To parallelise runs, [future.apply::future_lapply()] is implemented. See [future::future] and [future.apply::future.apply] for more information.
#' 
#' `k` is set to \eqn{[2, \lfloor\sqrt{n}\rfloor]} when NULL.
#' 
#' @returns a `pop` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of the cluster assignments. Each row is a solution.
#'   \item \code{summary}: a data frame of the summary statistics.
#'   \item \code{trace}: a data frame of summary statistics per iteration.
#'   \item \code{n_filtered}: a data frame of the numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} in [distmat_to_kmat()] when `space_kmat` is not supplied.
#' }
#' 
#' @seealso [distmat_to_kmat()], [sf::st_as_sf()], [lhs::randomLHS()], [mclust::adjustedRandIndex()], [future::future], [future.apply::future.apply]
#'
#' @export

blob_populate <- function(data,
                          k,
                          r = c(0.5,1),
                          iter = 10,
                          run = 100,
                          converge_ari = 1,
                          coords = c(1,2),
                          age = 3,
                          space_clustmethod = "kmedoids",
                          crs = 4326,
                          hull_convex_ratio = 0,
                          random_start = FALSE,
                          filter_intersects = TRUE,
                          filter_clustsize = TRUE,
                          max_na = 0.05,
                          space_kmat = NULL,
                          space_distmat = NULL,
                          space_distmethod = NULL,
                          w_knn = NULL,
                          l_normalise = NULL,
                          beta_par = NULL,
                          weights = NULL) {
  
  space_clustmethod <- match.arg(space_clustmethod, choices = c("kmedoids", "kkmeans"))
  #-------------------------------------------------------------------#
  # select the relevant columns
  data <- data[,c(coords,age)]
  #-------------------------------------------------------------------#
  if (space_clustmethod == "kkmeans") {
    # compute space_kmat
    if (is.null(space_kmat)) {
      if (is.null(w_knn)) w_knn <- 7
      if (is.null(l_normalise)) l_normalise <- TRUE
      if (is.null(beta_par)) beta_par <- 10
      #-------------------------------------------------------------------#
      if(is.null(space_distmat)) {
        if (is.null(space_distmethod)) {
          space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
          message(paste0(space_distmethod," is used to compute space_distmat"))
        } else {
          space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
        }
        space_kmat_out <- compute_kmat(data = data[, c(1,2)],
                                       method = space_distmethod,
                                       k = max(k),
                                       w_knn = w_knn,
                                       l_normalise = l_normalise,
                                       beta_par = beta_par)
      } else {
        space_kmat_out <- distmat_to_kmat(distmat = space_distmat,
                                          k = max(k),
                                          w_knn = w_knn,
                                          l_normalise = l_normalise,
                                          beta_par = beta_par)
      }
      space_kmat <- space_kmat_out$kmat
      space_kmat_optim_out <- space_kmat_out$optim_out
    } else {
      space_kmat_optim_out <- NULL # As it is one of the returned items
    }
  }
  
  if (space_clustmethod == "kmedoids"){
    if(is.null(space_distmat)) {
      if (is.null(space_distmethod)) {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
        message(paste0(space_distmethod," is used to compute space_distmat"))
      } else {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      }
      space_distmat <- compute_distmat(data = data[, c(1,2)],
                                       method = space_distmethod)
    }
    space_kmat_optim_out <- NULL # As it is one of the returned items
  }
  #-------------------------------------------------------------------#
  # sample r
  if (length(r) == 2) {
    # LHS sampling for more evenly distributed parameters
    lhs_samples <- lhs::randomLHS(run, 1)
    # scale to the range
    r_samples <- sort(as.vector(min(r) + lhs_samples * (max(r) - min(r))))
  } else {
    r_samples <- rep(r, run)
  }
  #-------------------------------------------------------------------#
  # assign k's lb and ub when NULL
  N <- nrow(data)
  if (is.null(k)) {
    k <- integer(2L)
    # lower bound
    k[1] <- 2L 
    # upper bound, an heuristic to maximise information e.g. 100 points: 10 blobs, 10 points each 
    k[2] <- as.integer(floor(sqrt(N)))
  }
  #-------------------------------------------------------------------#
  # for single k
  if (length(k) == 1) {
    # run blob_search() in parallel
    blob_list <- future.apply::future_lapply(r_samples, function (r) {
      blob_search(data = data,
                  k = k,
                  r = r,
                  iter = iter,
                  converge_ari = converge_ari,
                  coords = c(1,2),
                  age = 3,
                  space_clustmethod = space_clustmethod,
                  crs = crs,
                  hull_convex_ratio = hull_convex_ratio,
                  random_start = random_start,
                  filter_intersects = filter_intersects,
                  filter_clustsize = filter_clustsize,
                  max_na = max_na,
                  space_kmat = space_kmat,
                  space_distmat = space_distmat,
                  weights = weights) 
    }, future.seed = T)
    #-------------------------------------------------------------------#
    pop <- convert_to_pop(blob_list)
    pop$space_kmat_optim_out <- space_kmat_optim_out
  } else {
    # for a range of k
    if (length(k) == 2) {
      # create a grid of all combinations
      k_vec <- min(k):max(k)
      r_vec <- r_samples
      grid <- expand.grid(k_vec = k_vec, r_vec = r_vec)
      #-------------------------------------------------------------------#
      # run blob_search() in parallel for all combinations of parameters
      blob_list <- future.apply::future_Map(function(k, r) {
        list(
          blob = blob_search(data = data,
                             k = k,
                             r = r,
                             iter = iter,
                             converge_ari = converge_ari,
                             coords = coords,
                             age = age,
                             space_clustmethod = space_clustmethod,
                             crs = crs,
                             hull_convex_ratio = hull_convex_ratio,
                             random_start = random_start,
                             filter_intersects = filter_intersects,
                             filter_clustsize = filter_clustsize,
                             max_na = max_na,
                             space_kmat = space_kmat,
                             space_distmat = space_distmat,
                             weights = weights),
          k = k
        )
      }, grid$k_vec, grid$r_vec, future.seed = TRUE)
      #-------------------------------------------------------------------#
      # group runs of same k
      pop_list <- vector("list", k[2])
      for(i in 1:length(blob_list)) {
        m <- blob_list[[i]][["k"]]
        # append to pop_list
        pop_list[[m]] <- append(pop_list[[m]], list(blob_list[[i]][["blob"]]))
      }
      #-------------------------------------------------------------------#
      pop_list <- pop_list[-(1:(min(k)-1))] # as the loop starts at 2 and minus k [2,min(k)]
      # convert to pop object for each k group
      pop_list <- lapply(pop_list, convert_to_pop) 
      #-------------------------------------------------------------------#
      # extract each element and add a column to indicate the initial k
      pop <- do.call(rbind, pop_list)
      summary_list <- lapply(seq_along(pop[ , "summary"]), function (i) {
        if (!is.null(pop[, "summary"][[i]])) { cbind(pop[ , "summary"][[i]], k_o = i + min(k) - 1) }
      })
      trace_list <- lapply(seq_along(pop[ , "trace"]), function (i) {
        if (!is.null(pop[, "trace"][[i]])) { cbind(pop[ , "trace"][[i]], k_o = i + min(k) - 1) }
      }) 
      clust_list <- pop[ ,"clust"]
      n_filtered_list <- pop[ , "n_filtered"]
      #-------------------------------------------------------------------#
      # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
      clust <- do.call(rbind, clust_list)
      summary <- do.call(rbind, summary_list)
      trace <- do.call(rbind, trace_list)
      n_filtered <- do.call(rbind, n_filtered_list)
      n_filtered <- cbind(k_o = min(k):max(k), n_filtered)
      #-------------------------------------------------------------------#
      # reindex the solutions
      if (!is.null(clust)) {
        summary$idx <- NULL
        summary$idx <- 1:nrow(summary)
        trace$idx <- NULL
        trace <- merge(trace, summary[, c("idx","k_o","run")], by = c("k_o","run"), all.x = TRUE)
      } else {
        summary <- NULL
        trace <- NULL
      }
      #-------------------------------------------------------------------#
      summary <- summary[, c("idx", "k_o", "k", "r", "run", "space_wcd", "time_wcr", "time_wce", "intersects", "n_removed", "iter", "ari", "dup")]
      trace <- trace[, c("idx", "k_o", "k", "r", "run", "space_wcd", "time_wcr", "time_wce", "iter", "ari")]
      
      rownames(summary) <- NULL
      rownames(trace) <- NULL
      
      pop <- list(clust = clust,
                  summary = summary,
                  trace = trace,
                  n_filtered = n_filtered,
                  space_kmat_optim_out = space_kmat_optim_out)
    }
  }
  #-------------------------------------------------------------------#
  class(pop) <- "pop"
  return(pop)
}
