#' Core local-search algorithm
#' 
#' @description
#' `blob_search()` performs an iterative bi-objective local-search algorithm to
#' assign clusters for a given number of clusters (`k`) and spatial relative
#' weight (`r`).
#' 
#' @param data a data frame or matrix with spatial coordinates and age of the data.
#' @param k number of clusters.
#' @param r spatial relative weight of range \eqn{[0,1]}.
#' @param iter number of iterations. Default is `10L`.
#' @param converge_ari a numeric value of the Adjusted Rand Index (ARI) that
#' sets the convergence threshold between iterations. It must be of range
#' \eqn{[0,1]}. Default is `1`.
#' @param coords a vector of strings or integers indicating the columns of
#' spatial coordinates (e.g. `c("longitude, latitude)`). Default is `c(1L,2L)`
#' (the first and second columns).
#' @param age a string or an integer indicating the age of the data (e.g.
#' `age`). Default is `3L` (the 3rd column).
#' @param crs coordinate reference system passed on to [sf::st_as_sf()].
#' Default is `4326`.
#' @param space_distmat a spatial distance matrix or `dist` object.
#' Default is `NULL`.
#' @param space_distmethod spatial distance method used when `space_distmat`
#' is not specified. Either `"geodesic"` or `"euclidean"` is available (see
#' [compute_distmat()]). Default is `"geodesic"`.
#' @param filter_intersects a logical value to remove a solution with intersects
#' in space. Default is `TRUE`.
#' @param hull_convex_ratio a numeric value indicating the convexity of the
#' hulls passed onto [sf::st_concave_hull()]. `1` returns convex and `0`
#' maximally concave hulls. Default is `0.5`.
#' @param filter_clustsize a logical value to remove a solution with clusters
#' below the expected size. Default is `TRUE`.
#' @param outlier_removal (testing) a logical value to remove outliers. Default is `FALSE`.
#' @param outlier_iqrm (testing) multiplier of IQR for detecting outliers. Default is `1.5`.
#' @param random_init a logical value to choose random initial cluster points
#' instead of using a heuristic to choose points with approximately greatest 
#' separations. Default is `FALSE`.
#' @param weights a numeric vector of weights for each data point. Default is `NULL`.
#' @param sf_use_s2 a logical value to control spherical geometry.
#' See [sf::sf_use_s2()]. Default is `TRUE`.
#' 
#' @details
#' The core local-search algorithm searches for clusters that optimise
#' within-cluster spatial proximity and temporal coverage simultaneously.
#' Spatially, it uses k-medoid clustering approach (Kaufman & Rousseeuw, 1987),
#' implemented similarly to the standard k-means clustering algorithm 
#' (Lloyd, 1982). Temporally, a greedy algorithm is used to optimise range and
#' evenness of time points simultaneously. The spatial and temporal cost
#' functions are converted into a single cost function using a weighted sum
#' approach where the costs are normalised with respect to the costs across
#' \eqn{k} clusters by min-max normalisation. The relative spatial weight is
#' specified by the parameter `r`.
#' 
#' A partition is evaluated based on three objectives:
#' * `space_wcd` Within-cluster spatial sum of distances to their corresponding cluster
#'  medoids (total sum)
#' * `time_wcr` Within-cluster temporal range (average across clusters)
#' * `time_wce` Within-cluster temporal evenness (average across clusters)
#' 
#' The algorithm runs in an iterative fashion. The search is complete when the
#' adjusted Rand Index (ARI) reach the threshold specified by the parameter
#' `converge_ari` (`1` by default) after at least 2 iteration or until the set
#' length specified by the parameter `iter`.
#' 
#' An ideal partition would also results in clusters of roughly equal sizes. The
#' size threshold of a cluster is defined as \eqn{\frac{n}{2k}} where \eqn{n} is
#' the number of data points and \eqn{k} the number of clusters.
#' `filter_clustsize = T` validates the solution accordingly. Furthermore,
#' `filter_intersects = T` constrains feasible solutions to
#' partitions with non-overlapping boundaries. The boundaries are constructed by
#' concave/convex hulls, of which the convexity is controlled by
#' `hull_convex_ratio` passed onto the `ratio` of [sf::st_concave_hull()].
#'
#' `outlier_removal` is being tested at the moment. The idea is to remove points
#' whose distance to their medoid is considered an outlier based on the
#' interquantile range (IQR) of distances to medoid across \eqn{k} clusters.
#'
#' @returns
#' Either an S3 object of class `blob` with the following components:
#'  * `data`: a data frame of the input data with assigned clusters in column `clust`.
#'  * `summary`: a data frame of summary statistics. They include the returning
#'  number of clusters (`k`), original `k` parameter (`k_o`), three
#'  objective values (`space_wcd`, `time_wcr` and `time_wce`), number of
#'  iterations (`iter`), ARI with the previous iteration (`ari`), presence
#'  of intersecting clusters (`intersects`) and number of clusters flagged
#'  for expected size (`clustsize_f`).
#'  * `trace`: a data frame of summary statistics across iterations.
#'  * `params`: a list of parameter values.
#' or a status code if
#'  * `1`: degenerate case when there is only 1 returning cluster.
#'  * `2`: the solution has intersecting clusters with `filter_intersects == T`.
#'  * `3`: the solution has at least a cluster with lower than expected size
#'  with `filter_clustsize == T`.
#'  * `4`: both cases `2` and `3`.
#'
#' @seealso [compute_distmat()], [sf::st_as_sf()], [sf::st_concave_hull()],
#' [sf::sf_use_s2()], [mclust::adjustedRandIndex()]
#' 
#' @references
#' S Lloyd. “Least squares quantization in PCM”. en. In: IEEE Trans. Inf. Theory
#' 28.2 (Mar.1982), pp. 129–137.
#' 
#' L Kaufman and P Rousseeuw. “Clustering by means of medoids”. In:
#' International Conference on Statistical Data Analysis Based on the L1-Norm
#' and Related Methods (1987), pp. 405–416.
#'
#' @export

blob_search <- function(data,
                        k,
                        r,
                        iter = 10L,
                        converge_ari = 1,
                        coords = c(1L,2L),
                        age = 3L,
                        crs = 4326,
                        space_distmat = NULL,
                        space_distmethod = NULL,
                        filter_intersects = TRUE,
                        hull_convex_ratio = 0.5,
                        filter_clustsize = TRUE,
                        outlier_removal = FALSE,
                        outlier_iqrm = 3,
                        random_init = FALSE,
                        weights = NULL,
                        sf_use_s2 = TRUE) {
  
  # checks
  check_input_bs(data = data,
                 k = k,
                 r = r,
                 iter = iter,
                 converge_ari = converge_ari,
                 coords = coords,
                 age = age,
                 space_distmat = space_distmat,
                 filter_intersects = filter_intersects,
                 hull_convex_ratio = hull_convex_ratio,
                 filter_clustsize = filter_clustsize,
                 random_init = random_init,
                 weights = weights,
                 sf_use_s2 = sf_use_s2)
  
  params <- mget(c("k", "r", "iter", "converge_ari", "filter_intersects",
                   "hull_convex_ratio", "filter_clustsize",
                   "outlier_removal", "weights"))
  
  # select the relevant columns
  data <- as.data.frame(data)
  # check if space_distmat is supplied
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod)
  
  # init_blobs() to initialise medoids
  data <- init_blobs(data = data, k = k, space_distmat = space_distmat, random_init = random_init)
  # initialise counter counting find_blobs()
  t <- 0
  # initialise trace table
  trace <- data.frame()
  
  for (i in 1:iter) {
    data_old <- data
    
    # find_blobs()
    data <- find_blobs(data = data, k = k, r = r,
                       space_distmat = space_distmat,
                       age = age,
                       weights = weights)
    
    # count find_blobs() executed
    t <- t + 1
    
    if (t > 0) {
      # check convergence by ARI
      ari <- mclust::adjustedRandIndex(data$clust, data_old$clust)
      
      # evaluate the output
      summary <- eval_blob(data = data,
                           crs = crs,
                           hull_convex_ratio = hull_convex_ratio,
                           space_distmat = space_distmat,
                           weights = weights,
                           sf_use_s2 = sf_use_s2)
      
      # other summary columns
      summary$k_o <- as.integer(k)
      summary$r <- r
      summary$ari <- ari
      summary$iter <- as.integer(t)
      # append trace row
      trace_row <- summary
      trace <- rbind(trace, trace_row)
     
      # if converged between t and t-1, break
      if (!is.null(converge_ari)) {
        if (all(ari >= converge_ari) == TRUE & t > 2L) break
      }
    }
  }

  # remove outliers
  if (outlier_removal == TRUE) {
    out <- remove_outliers(data = data,
                           space_distmat = space_distmat,
                           space_distmethod = space_distmethod,
                           weights = weights,
                           iqrm = outlier_iqrm)
    data <- out$data
    n_outliers <- out$n_outliers
  } else {
    n_outliers <- 0
  }

  summary$n_outliers <- n_outliers
  
  # filter solution with intersects and clustsize_f
  status_nclust <- status_intersects <- status_clustsizef <- 0
  
  if (length(unique(data$clust)) < 2L) {
    message("Status code 1: Only 1 returning cluster.")
    return (1L) 
  }
  
  if (filter_intersects == TRUE) {
    if (summary$intersects == T) status_intersects <- 1
  }
  if (filter_intersects == TRUE) {
    if (summary$intersects == T) status_intersects <- 1
  }
  if (filter_clustsize == TRUE) {
    if (summary$clustsize_f > 0) status_clustsizef <- 1
  }
  
  if (sum(status_intersects, status_clustsizef) == 1) {
    if (status_intersects == 1) {
      message("Status code 2: Intersecting clusters.")
      return(2L)
    }
    
    if (status_clustsizef == 1) {
      message("Status code 3: Cluster size below threshold.")
      return(3L)
    }
  } 
  
  if (sum(status_intersects, status_clustsizef) == 2) {
    message("Status code 4: Intersecting clusters & cluster size below threshold.")
    return(4L)
  }
  
  clust <- reorder_clust(data$clust)
  data$clust <- NULL
  # reorder the columns
  summary <- select_summary_bs(summary)
  trace <- select_trace_bs(trace)
  
  return(new_blob(clust = clust, summary = summary, trace = trace, data = data, params = params))
}

# init_blobs ------------------------------------------------------------------
init_blobs <- function(data, space_distmat, k, random_init = FALSE) {
  # initialise clust
  N <- nrow(data)
  clust <- rep(NA, N)

  if (random_init == TRUE) {
    init <- sample(1:N, k)
    clust[init] <- 1:k
    return(data)
  }

  # start from k roughly equally spaced random locations 
  # (just permute a bunch and pick the one with the least smallest distance)
  m <- 100
  mat <- matrix( , m, k)
  min_dist <- numeric(m)
  
  # loop through 100 randomly sampled points
  for(i in 1:m){
    # sample k points from the data, sort them 
    j <- sort(sample(1:N, size = k))
    mat[i, ] <- j
    # all combinations of k chooses 2 
    cb <- utils::combn(k, 2)
    nc <- ncol(cb)
    dist <- numeric(nc)
    # extract from distance matrix the distances for all combinations of points
    for (c in 1:nc) dist[c] <- space_distmat[j[cb[1, c]], j[cb[2, c]]]
    min_dist[i] <- min(dist)
  }
  
  # pick the set with the largest minimum distance between two points
  # pick the first one if two are tied 
  init <- mat[which(min_dist == max(min_dist))[1], ]
  
  # assign cluster to the initial points
  clust[init] <- 1:k
  data$clust <- clust

  return(data)
}

# find_blobs ------------------------------------------------------------------
find_blobs <- function(data, k, r, space_distmat, age = 3L, weights = NULL) {
  N <- nrow(data)
  clust <- data$clust
  ages <- data[ , age]
  
  # initialise vectors for the loop
  space_costs <- time_costs <- n <- rep(NA, k)
  
  # precompute space_cost
  space_costmat <- vapply(1:k, function(j) {
    clust_points <- which(clust == j)
    compute_space_costs(space_distmat = space_distmat,
                        clust_points = clust_points,
                        weights = weights)
  }, FUN.VALUE = numeric(N))
  
  # loop through every point (incremental updating)
  for (i in sample(1:N)) {
    # loop through clusters
    for (j in 1:k) {
      # points in clust j
      clust_points <- which(clust == j)
      # number of points in clust j
      n[j] <- length(clust_points)
      # next cluster if there is no point
      if (n[j] == 0) next
      
      # store space_cost for point i in clust j
      space_costs[j] <- space_costmat[i, j]
      
      # get the set of clust_points excluding point i
      clust_points_tmp <- if (i %in% clust_points) clust_points[clust_points != i] else clust_points
      
      # compute the temporal cost for point i in clust j
      if (length(clust_points_tmp) == 0) {
        # if only point i is in the cluster
        time_costs[j] <- 0  
      } else {
        time_costs[j] <- min(abs(ages[i] - ages[clust_points_tmp]))
      }
    }
    
    # normalise the cost by min-max normalisation
    space_costs_norm <- (space_costs - min(space_costs, na.rm = T)) / (max(space_costs, na.rm = T) - min(space_costs, na.rm = T))
    time_costs_norm <- (time_costs - min(time_costs, na.rm = T)) / (max(time_costs, na.rm = T) - min(time_costs, na.rm = T))
    space_costs_norm[is.na(space_costs_norm)] <- 0
    time_costs_norm[is.na(time_costs_norm)] <- 0
    
    # single cost by weighted-sum scalarisation
    costs <- space_costs_norm*r - time_costs_norm*(1-r)
    # cluster(s) with the minimum cost
    clust_tmp <- which(costs == min(costs, na.rm = T))
    # if there is more than one candidate cluster
    if (length(clust_tmp) > 1) {
      # check if one has fewer points
      n_clust_tmp <- n[clust_tmp]
      clust_tmp_minn_idx <- which(n_clust_tmp == min(n_clust_tmp))
      
      if (length(clust_tmp_minn_idx) > 1) {
        # randomly pick one if tied
        clust_tmp_minn_idx <- sample(clust_tmp_minn_idx, 1)
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      } else {
        # else pick the cluster with fewer points in it 
        clust[i] <- clust_tmp[clust_tmp_minn_idx]
      }
    } else {
      # else just assign the cluster
      clust[i] <- clust_tmp
    }
  }
  
  # update clust in the output
  data$clust <- clust
  
  return(data)
}

# eval_blob ------------------------------------------------------------------
eval_blob <- function(data,
                      coords = c(1L,2L),
                      age = 3L,
                      crs = 4326,
                      hull_convex_ratio = 0.5,
                      space_distmat = NULL,
                      space_distmethod = NULL,
                      weights = NULL,
                      sf_use_s2 = TRUE) {
  # total number of points
  N <- sum(!is.na(data$clust))
  # clust
  clust <- data$clust
  # total number of clusters
  k <- length(unique(stats::na.omit(data$clust))) # NA is excluded
  # initialise empty vectors 
  space_d <- time_r <- time_e <- n <- rep(NA, k)
  
  # check if is.null(space_distmat)
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod)
  ages <- data[ , age]
  
  # for weighted k medoids
  weights <- weights %||% rep(1, nrow(data))
  
  # loop through k to obtain within cluster statistics
  for (j in 1:k) {
    clust_points <- which(clust == j)
    # skip if there is no point in the cluster
    if (length(clust_points) == 0) next
    
    ages_j <- ages[clust_points]
    
    # spatial objective
    space_costs_j <- compute_space_costs(space_distmat = space_distmat,
                                         clust_points = clust_points,
                                         weights = weights)[clust_points]
    
    space_d[j] <- sum(weights[clust_points] * space_costs_j)
    
    # temporal objectives
    time_r[j] <- max(ages_j) - min(ages_j) 
    time_e[j] <- 1 / (1 + stats::var(diff(sort(ages_j)))) # NA if there are fewer than 3 data points
    
    # clust size
    n[j] <- length(clust_points)
  }
  
  # calculate the summary statistics
  space_wcd <- sum(space_d)
  time_wcr <- mean(time_r, na.rm = TRUE)
  time_wce <- mean(time_e, na.rm = TRUE)
  # na.rm = TRUE to assess the overall performance
  # as we do not care clusters with fewer than 3 points
  
  # flag the number of clusters below the threshold
  clustsize_f <- length(which(n < N/k/2))
  
  # evaluate if blobs are intersecting in space
  intersects <- check_intersects(data = data, coords = coords, crs = crs,
                                 hull_convex_ratio = hull_convex_ratio,
                                 sf_use_s2 = sf_use_s2)
  
  # return a data frame of all the statistics
  summary <- data.frame(k = k,
                        space_wcd = space_wcd,
                        time_wcr = time_wcr,
                        time_wce = time_wce,
                        intersects = intersects,
                        clustsize_f = clustsize_f)
  return(summary)
}

# print.blob ------------------------------------------------------------------
print.blob <- function(x, ...) {
  cat("STblob solution\n")
  cat("# clusters:\n")
  print(head(x$clust, 20))
  if (length(x$clust) > 20) cat(paste0(" <", length(x$clust) - 20," more points>\n"))
  cat("\n")
  cat("# summary:\n")
  print(x$summary)
  cat("\n")
  cat("# trace:\n")
  print(head(x$trace, 5))
  if (nrow(x$trace) > 5) cat(paste0(" <", nrow(x$trace) - 5," more rows>\n"))
  invisible(x)
}

# summary.blob ----------------------------------------------------------------
summary.blob <- function(x, ...) {
  print(x$summary)
  invisible(x)
}

# helpers ---------------------------------------------------------------------
## remove_outliers ------------------------------------------------------------
# testing
# find_outliers <- function(data, coords = c(1L, 2L),
#                             space_distmat = NULL, space_distmethod = NULL,
#                             knn = 3) {
#   # check space_distmat
#   space_distmat <- check_space_distmat(data = data,
#                                        coords = coords,
#                                        space_distmat = space_distmat,
#                                        space_distmethod = space_distmethod)
#   
#   # mean distances to the k nearest neighbours
#   mu <- colMeans(apply(space_distmat, 2, sort)[2:(knn+1), ])
#   q <- quantile(mu, c(0.25, 0.75)) # Q1 and Q3
#   # interquartile range
#   iqr <- q[2] - q[1]
#   # upper limit
#   ul <- q[2] + 1.5 * iqr 
#   outliers <- which(mu > ul)
#   
#   return(outliers)
# }

# This version is based on the distances to medoid and remove them adhoc
remove_outliers <- function(data,
                            coords = c(1L, 2L),
                            space_distmat = NULL,
                            space_distmethod = NULL,
                            weights = NULL,
                            iqrm = 1.5) {

  # clust
  clust <- data$clust
  # total number of clusters
  k <- length(unique(stats::na.omit(data$clust))) # NA is excluded

  # check if is.null(space_distmat)
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod)

  # remove outliers
  space_costs <- numeric()
  order <- numeric()
  for (i in 1:k) {
    clust_points <- which(clust == i)
    # skip if there is no point in the cluster
    if (length(clust_points) == 0) next
    space_costs <- append(space_costs,
                          compute_space_costs(space_distmat = space_distmat,
                                              clust_points = clust_points,
                                              weights = weights)[clust_points])
    order <- append(order, clust_points)
  }
  space_costs <- space_costs[order(order)]
  # https://www.geeksforgeeks.org/data-analysis/what-is-outlier-detection/
  # https://www.geeksforgeeks.org/machine-learning/interquartile-range-to-detect-outliers-in-data/
  q <- quantile(space_costs, c(0.25, 0.75)) # Q1 and Q3
  # interquartile range
  iqr <- q[2] - q[1]
  # upper limit
  # ub <- q[2] + 1.5 * iqr
  ub <- q[2] + iqrm * iqr
  outliers <- which(space_costs > ub)
  # outliers <- which(space_costs > quantile(space_costs, probs = quantile_prob))
  n_outliers <- length(outliers)
  data$clust[outliers] <- NA

  out <- list(data = data, n_outliers = n_outliers)
  return(out)
}

## check_intersects -----------------------------------------------------------
check_intersects <- function(data, clust = NULL, coords = c(1,2),
                             crs = 4326, hull_convex_ratio = 0.5,
                             sf_use_s2 = TRUE) {
  # st_union breaks when s2 is on, switch it off until on.exit()
  old <- suppressMessages(sf::sf_use_s2(sf_use_s2))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  
  data <- as.data.frame(data)
  
  if (!is.null(clust)) data$clust <- clust
  data_sf <- sf::st_as_sf(data, coords = coords, crs = crs)
  
  # obtain convex/concave hulls
  # If by_feature is TRUE each feature geometry is unioned individually.
  # This can for instance be used to resolve internal boundaries after
  # polygons were combined using st_combine.
  # https://r-spatial.github.io/sf/reference/geos_combine.html
  suppressMessages({
    hulls <- stats::aggregate(data_sf$geometry,
                              by = list(clust = data_sf$clust),
                              function(x){
                                x <- sf::st_combine(x)
                                x <- sf::st_union(x, by_feature = TRUE)
                                return(x)
                              }) 
    hulls <- sf::st_as_sf(hulls)
    hulls <- sf::st_concave_hull(hulls, ratio = hull_convex_ratio)
    
    # mark TRUE/FALSE for each comparison
    hulls <- sf::st_make_valid(hulls)
    intersects <- sf::st_intersects(hulls$geometry, sparse = F)
    diag(intersects) <- NA
  })
 
  out <- if (any(intersects == TRUE, na.rm = TRUE)) TRUE else FALSE
  return(out)
}

## check_input_bs -------------------------------------------------------------
check_input_bs <- function(data,
                           k,
                           r,
                           iter,
                           converge_ari,
                           coords,
                           age,
                           space_distmat,
                           filter_intersects,
                           hull_convex_ratio,
                           filter_clustsize,
                           random_init,
                           weights,
                           sf_use_s2) {
  stopifnot(
    "`data` must be a data frame or matrix" = is.data.frame(data) || is.matrix(data),
    "`data` must have at least 3 columns" = ncol(data) >= 3L,
    "`k` must be greater than 1" = k > 1L,
    "`r` must be between 0 and 1" = r >= 0 && r <= 1,
    "`iter` must be at least 3" = iter >= 3L,
    "`converge_ari` must be between 0 and 1" = converge_ari >= 0 && converge_ari <= 1,
    "`filter_intersects` must be logical" = is.logical(filter_intersects),
    "`hull_convex_ratio` must be between 0 and 1" = hull_convex_ratio >= 0 && hull_convex_ratio <= 1,
    "`filter_clustsize` must be logical" = is.logical(filter_clustsize),
    "`random_init` must be logical" = is.logical(random_init),
    "`sf_use_s2` must be logical" = is.logical(sf_use_s2)
  )
  
  check_col_ref(coords, data, "coords", expected_length = 2L)
  check_col_ref(age, data, "age", expected_length = 1L)

  n <- if (!is.null(space_distmat)) nrow(space_distmat) else nrow(data)
  
  if (!is.null(space_distmat)) {
    stopifnot("`space_distmat` must be a numeric matrix or a `dist` object" =
                is.numeric(space_distmat),
              "`space_distmat` must have the same number of rows as `data`" =
                nrow(data) == nrow(space_distmat))
  }
  if (!is.null(weights)) {
    stopifnot("`weights` must be numeric" =
                is.numeric(weights),
              "`weights` must have length equal to number of rows in `data`" =
                length(weights) == n)
  }
  
  invisible(NULL)
}

## select_bs ------------------------------------------------------------------
select_summary_bs <- function(x) {
  # select blob_search() summary columns
  x <- x[ , c("k", "k_o", "r", "space_wcd", "time_wcr", "time_wce", "iter", "ari",
              "intersects", "clustsize_f", "n_outliers")]
  rownames(x) <- NULL
  return(x)
}

select_trace_bs <- function(x) {
  # select blob_search() summary columns
  x <- x[ , c("k", "k_o", "r", "space_wcd", "time_wcr", "time_wce", "iter", "ari")]
  rownames(x) <- NULL
  return(x)
}

## new_blob -------------------------------------------------------------------
new_blob <- function(clust, data, summary, trace, params) {
  stopifnot(is.numeric(clust),
            is.data.frame(summary),
            is.data.frame(trace),
            is.data.frame(data),
            is.list(params))
  
  structure(
    list(clust = clust, summary = summary, trace = trace, data = data, params = params),
    class = "blob"
  )
}