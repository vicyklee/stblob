#' Core local-search algorithm
#' 
#' @description
#' `stblob_lsearch()` performs a bi-objective local-search algorithm to
#' assign clusters for a given number of clusters (`k`).
#' 
#' @param data a data frame or matrix with spatial coordinates, age and optionally
#' type of the data.
#' @param k number of clusters.
#' @param w_space relative spatial weight of range \eqn{\[0,1\]}.
#' @param optim_type_diversity a logical value to optimise data type diversity.
#' Default is `FALSE`.
#' @param w_time relative temporal weight of range \eqn{\[0,1\]}. Default is `NULL`.
#' @param w_type relative data type diversity weight of range \eqn{\[0,1\]}. Default is `NULL`.
#' @param iter number of iterations. Default is `10L`.
#' @param ls_tol tolerance for local-search optimisation. The difference to the
#' highest Adjusted Rand Index (ARI) (i.e. `1` -- identitical) for three
#' consecutive iterations to be considered reaching convergence. Default is `0`.
#' @param coords a vector of character strings of length `2L`.
#' Default is `NULL` with the 2nd and 3rd columns of `data` as longitude and latitude. 
#' @param age a character string. Default is `NULL` with the 4th column of `data`.
#' @param type a character string. Default is `NULL` with the 5th column of `data`.
#' @param space_distmat a spatial distance matrix or `dist` object.
#' Default is `NULL`.
#' @param space_distmethod spatial distance method used when `space_distmat`
#' is not specified. Either `"geodesic"` or `"euclidean"` is available (see
#' [compute_distmat()]). Default is `"geodesic"`.
#' @param filter_intersects a logical value to remove a solution with intersects
#' in space. Default is `TRUE`.
#' @param hull_convex_ratio a numeric value indicating the convexity of the
#' hulls passed onto [sf::st_concave_hull()] for checking intersects. `1`
#' returns convex and `0` maximally concave hulls. Default is `0.5`.
#' @param hull_crs coordinate reference system passed on to [sf::st_as_sf()] for
#' checking intersects. Default is `4326`.
#' @param filter_clustsize a logical value to remove a solution with clusters
#' below the expected size. Default is `TRUE`.
#' @param outlier_removal (testing) a logical value to remove outliers. Default is `FALSE`.
#' @param outlier_iqrm (testing) multiplier of IQR for detecting outliers. Default is `1.5`.
#' @param chebyshev_rho (testing) weight \eqn{\rho} in the weighted sum
#' component of augmented Chebyshev scalarisation. Default is `1e-4`.
#' @param outlier_iqrm (testing) multiplier of IQR for detecting outliers. Default is `1.5`.
#' @param random_init a logical value to choose random initial cluster points
#' instead of using a heuristic to choose points with approximately greatest 
#' separations. Default is `FALSE`.
#' @param sf_use_s2 a logical value to control spherical geometry.
#' See [sf::sf_use_s2()]. Default is `TRUE`.
#' 
#' @details
#' The core local-search algorithm searches for clusters that optimise
#' within-cluster spatial proximity, temporal coverage, and data type
#' information simultaneously.
#' Spatially, it uses k-medoid clustering approach (Kaufman & Rousseeuw, 1987),
#' implemented similarly to the standard k-means clustering algorithm 
#' (Lloyd, 1982). For time and data type diversity, greedy algorithms are
#' used to optimise within-cluster temporal range and evenness, as well as data
#' types contained in each cluster.
#' 
#' The cost functions are combined into a single function using augmented
#' weighted Chebyshev (Tchebycheff) (Steuer & Choo, 1983).
#' To ensure the costs are evaluated on comparable scales as specified by the
#' relative weights with the parameters `w_space`, `w_time` and `w_type`,
#' they are normalised with respect to the costs across \eqn{k} clusters at each 
#' step by min-max normalisation.
#' 
#' A partition is evaluated based on three objectives:
#'
#' * `space_distance`: average within-cluster spatial distances to the medoid
#' * `time_variance`:  average within-cluster temporal variance
#' * `time_evenness`: average within-cluster temporal variance of first differences
#' * `type_diversity`: average within-cluster Shannon entropy of data types
#' 
#' The objective values are returned in the `summary` by the terms `z_space`,
#' `z_time1`, `z_time2` and `z_type` respectively.
#' The maximisation objectives are multiplied by -1 as a minimisation problem
#' in multi-objective optimisation.
#' 
#' The search will terminate until either the ARI tolerance to the previous
#' iteration are at most as specified by `ls_tol` for three consecutive
#' iterations or the maximum iteration is reached as specified by `iter`.
#' 
#' An ideal partition should return clusters of roughly equal sizes. The
#' size threshold of a cluster is defined as \eqn{\frac{n}{2k}} where \eqn{n} is
#' the number of data points and \eqn{k} the number of clusters.
#' `filter_clustsize = T` will remove the solution with any flagged clusters.
#' `filter_intersects = T` constrains feasible solutions to
#' partitions with non-overlapping boundaries. The boundaries are constructed by
#' concave/convex hulls, of which the convexity is controlled by
#' `hull_convex_ratio` passed onto the `ratio` of [sf::st_concave_hull()].
#'
#' `outlier_removal` is being tested at the moment. The idea is to remove points
#' whose distance to their medoid is an outlier based on the
#' interquantile range (IQR) of distances to medoid across \eqn{k} clusters.
#'
#' @returns
#' an S3 object of class `stblob_sol` with the following components:
#'  * `clust`: a vector of cluster assignments.
#'  * `summary`: a data frame of summary statistics. They include the parameter
#'  value of (`k`), number of clusters of the output (`k_o`), three
#'  objective values (`z_space`, `z_time1`, `z_time2` and `z_type`),
#'  number of iterations (`iter`), ARI with the previous iteration (`ari`),
#'  state of local search convergence (`ls_convergence`)
#'  and number of outliers (`n_outliers`).
#'  * `trace`: a data frame of summary statistics for each iteration.
#'  * `data`: a data frame of the input data.
#'  * `status_code`:
#'    * `1`: degenerate case when there is only 1 returning cluster.
#'    * `2`: the solution has intersecting clusters
#'    * `3`: the solution has at least a cluster with lower than expected size
#'    * `4`: both cases `2` and `3`.
#'  * `params`: a list of input parameter values.
#'
#' @seealso [compute_distmat()], [sf::st_as_sf()], [sf::st_concave_hull()],
#' [sf::sf_use_s2()], [aricode::ARI()]
#' 
#' @references
#' Kaufman, L., & Rousseeuw, P. (1987). Clustering by means of medoids.
#' International Conference on Statistical Data Analysis Based on the L1-Norm
#' and Related Methods, 405–416.
#' 
#' Lloyd, S. (1982). Least squares quantization in PCM. IEEE Transactions on
#' Information Theory, 28(2), 129–137.
#' 
#' Steuer, R. E., & Choo, E.-U. (1983). An interactive weighted Tchebycheff
#' procedure for multiple objective programming. Mathematical Programming,
#' 26(3), 326–344.
#'
#' @export

stblob_lsearch <- function(data,
                           k,
                           w_space, 
                           optim_type_diversity,
                           w_time = NULL,
                           w_type = NULL,
                           iter = 10L,
                           coords = NULL, 
                           age = NULL, 
                           type = NULL, 
                           space_distmat = NULL,
                           space_distmethod = c("geodesic", "euclidean"),
                           ls_tol = 0,
                           chebyshev_rho = 1e-4, 
                           filter_intersects = TRUE,
                           hull_convex_ratio = 0.5,
                           hull_crs = 4326,
                           filter_clustsize = TRUE, 
                           outlier_removal = FALSE,
                           outlier_iqrm = 1.5,
                           random_init = FALSE,
                           sf_use_s2 = TRUE) {
  
  # return params as a list element
  params <- mget(ls(environment(), sorted = FALSE))
  
  # check data
  stopifnot("data must be a data.frame" = is.data.frame(data))
  data_input <- data
  
  # check columns
  coords <- if (!is.null(coords)) match(coords, names(data)) else names(data)[c(2,3)]
  age <- if (!is.null(age)) match(age, names(data)) else names(data)[4]
  type <- if (!is.null(type)) match(type, names(data)) else names(data)[5]
  
  # type may not be a column of optim_type_diversity is turned off
  if (!optim_type_diversity) type <- NULL
  
  stopifnot(
    "'coords' do not match any column names." = !any(is.na(coords)),
    "'age' does not match any column names." = !any(is.na(age)),
    "'type' does not match any column names." = !any(is.na(type)),
    "'coords' must be numeric." = is.numeric(data[[coords[1]]]) & is.numeric(data[[coords[2]]]),
    "'age' must be numeric." = is.numeric(data[[age]])
  )
  
  # calculate relative weights for NULL weights
  w_type <- w_type %||% 0
  w_time <- w_time %||% (1 - w_space - w_type)
  # floating point imprecision for 0
  if (abs(w_time) < .Machine$double.eps ^ 0.5) w_time <- 0
  
  # check w_space, w_time and w_type
  stopifnot(
    "'w_space', 'w_time' and 'w_type' must be between 0 and 1" =
      (w_space >= 0 && w_space <= 1) && (w_time >= 0 && w_time <= 1) && (w_type >= 0 && w_type <= 1),
    "'w_space', 'w_time' and 'w_type' must sum up to 1" =
      all.equal(sum(w_space, w_time, w_type), 1)
  )
  
  # check if space_distmat is supplied
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod)
  
  # in case crs is not specified otherwise
  if (space_distmethod == "euclidean") hull_crs <- NA
  
  check_opt_args(ls_tol = ls_tol,
                 chebyshev_rho = chebyshev_rho, 
                 filter_intersects = filter_intersects,
                 hull_convex_ratio = hull_convex_ratio,
                 hull_crs = hull_crs,
                 filter_clustsize = filter_clustsize, 
                 outlier_removal = outlier_removal,
                 outlier_iqrm = outlier_iqrm,
                 random_init = random_init,
                 sf_use_s2 = sf_use_s2)
  
  # init_lsearch() to initialise medoids
  clust <- lsearch_init(data = data,
                        k = k,
                        space_distmat = space_distmat,
                        random_init = random_init)
  
  # initialise trace table
  trace <- data.frame()

  # main loop
  for (t in 1:iter) {
    clust_prev <- clust
    # lsearch_assign
    clust <- lsearch_assign(data = data,
                            clust = clust,
                            k = k,
                            w_space = w_space,
                            w_time = w_time,
                            w_type = w_type,
                            space_distmat = space_distmat,
                            age = age,
                            type = type,
                            optim_type_diversity = optim_type_diversity,
                            chebyshev_rho = chebyshev_rho)
    
    # check convergence by ARI
    # assign 0 to NA for aricode::ARI
    clust[is.na(clust)] <- clust_prev[is.na(clust_prev)] <- 0
    ari <- if (t < 2) NA else aricode::ARI(clust, clust_prev)
    
    # evaluate the output
    summary <- eval_sol(data = data,
                        clust = clust,
                        coords = coords,
                        age = age,
                        type = type,
                        optim_type_diversity = optim_type_diversity,
                        hull_convex_ratio = hull_convex_ratio,
                        hull_crs = hull_crs,
                        space_distmat = space_distmat,
                        sf_use_s2 = sf_use_s2)
    
    # other summary columns
    summary$k <- as.integer(k)
    summary$w_space <- w_space
    summary$w_time <- w_time 
    summary$w_type <- w_type
    summary$ari <- ari
    summary$iter <- as.integer(t)
    summary$ls_convergence <- FALSE
    # append trace row
    trace_row <- summary
    trace <- rbind(trace, trace_row)
    
    # if converged between t and t-1 for three consecutive iter, break
    if (!is.null(ls_tol)) {
      if (t > 3) {
        if (all(tail(trace$ari, 3) >= 1 - ls_tol)) {
          summary$ls_convergence <- TRUE
          break
        } 
      }
    } 
  }
  
  # remove outliers
  if (outlier_removal) {
    o <- remove_outliers(data = data,
                         clust = clust,
                         space_distmat = space_distmat,
                         space_distmethod = space_distmethod,
                         iqrm = outlier_iqrm)
    clust <- o$clust
    n_outliers <- o$n_outliers
    # update summary after removing outliers
    summary_updated <- eval_sol(data = data,
                                clust = clust,
                                age = age,
                                type = type,
                                hull_convex_ratio = hull_convex_ratio,
                                hull_crs = hull_crs,
                                space_distmat = space_distmat,
                                sf_use_s2 = sf_use_s2)
    summary[names(summary_updated)] <- summary_updated
  } else {
    n_outliers <- 0
  }
  
  # add column n_outliers to summary
  summary$n_outliers <- n_outliers
  
  # filter solution with intersects and clustsize_f
  status <- status_intersects <- status_clustsizef <- 0
  
  if (length(unique(stats::na.omit(clust))) < 2) {
    message("Status code 1: Only 1 returning cluster.")
    status <- 1
  }
  
  if (filter_intersects) {
    status_intersects <- if (summary$intersects) 1
  }
  
  if (filter_clustsize) {
    status_clustsizef <- if (summary$clustsize_f > 0) 2
  }
  
  status_ic_sum <- sum(status_intersects, status_clustsizef)
  
  if (status_ic_sum > 0) {
    switch(
      status_ic_sum,
      { message("Status code 2: Intersecting clusters."); status <- 2 },
      { message("Status code 3: Cluster size below threshold."); status <- 3 },
      { message("Status code 4: Intersecting clusters & cluster size below threshold."); status <- 4 }
    )
  }
  
  if (status > 0) {
    clust <- summary <- trace <- NA
  } else {
    clust <- reorder_clust(clust)
    data <- data_input
    summary <- select_summary(summary) 
    trace <- select_trace(trace)
  }

  
  return(new_sol(clust = clust, summary = summary, trace = trace,
                 data = data, status = status, params = params))
}

# lsearch_init ----------------------------------------------------------------
lsearch_init <- function(data, space_distmat, k, random_init = FALSE) {
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

  return(clust)
}

# lsearch_assign -------------------------------------------------------------- 
lsearch_assign <- function(data,
                           clust,
                           k, 
                           w_space,
                           w_time,
                           w_type,
                           space_distmat,
                           optim_type_diversity,
                           age,
                           type,
                           chebyshev_rho) {
  
  N <- nrow(data)
  ages <- data[[age]]
  
  # initialise vectors for the loop
  space_costs <- time_costs <- n <- rep(NA, k)
  
  if (optim_type_diversity) {
    types <- data[[type]]
    type_costs <- rep(NA, k)
  }
  
  # precompute space_costs
  space_costmat <- vapply(1:k, function(j) {
    clust_points <- which(clust == j)
    compute_space_costs(space_distmat = space_distmat, clust_points = clust_points)
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
    
      # space
      # store space_cost of point i in clust j
      space_costs[j] <- space_costmat[i, j]
      
      # time
      # get the set of clust_points excluding point i
      clust_points_excl_i <- clust_points[clust_points != i]
      # temporal cost function
      clust_points_excl_i_type <- if (optim_type_diversity) {
        # include clust_points only of the same type
        clust_points_excl_i[types[clust_points_excl_i] == types[i]]
      } else { clust_points_excl_i }
      
      if (length(clust_points_excl_i_type) == 0) {
        # if only point i is in the cluster, assign max value
        # time_costs can only be assigned after the loop so NA for now
        time_costs[j] <- NA
      } else {
        # compute temporal cost for point i in clust j
        time_costs[j] <- min(abs(ages[i] - ages[clust_points_excl_i_type]))
      }
      
      # type diversity
      if (optim_type_diversity) {
        if (length(clust_points_excl_i) == 0) {
          type_costs[j] <- 0
        } else {
          # compute type diversity cost of point i in clust j
          # type sensitive
          clust_types_excl_i <- types[clust_points_excl_i]
          # number of data of the same type
          type_costs[j] <- sum(clust_types_excl_i == types[i])
        } 
      } else {
        # if information costs are irrelevant # w_info * info_costs = 0
        type_costs <- 0
      }
    }
    
    if (all(is.na(time_costs))) {
      # when type sensitive, it is possible there is no same type across all clusters
      time_costs[is.na(time_costs)] <- 0
    } else {
      # assign max value to the cluster without a single data point
      time_costs[is.na(time_costs)] <- max(time_costs, na.rm = TRUE)
    }
    
    # normalise the cost by min-max normalisation
    space_costs_norm <- minmax(space_costs)
    time_costs_norm <- minmax(time_costs, max = TRUE)
    type_costs_norm <- minmax(type_costs)
    
    # when max == min NA, denominator = 0, costs should be ideal
    utop <- -0.1 # utopia value
    space_costs_norm[is.na(space_costs_norm)] <-
      time_costs_norm[is.na(time_costs_norm)] <-
      type_costs_norm[is.na(type_costs_norm)] <- 0
    
    # https://www.austintripp.ca/blog/2025-05-12-chebyshev-scalarization/
    costs <- pmax(w_space * abs(space_costs_norm - utop),
                  w_time * abs(time_costs_norm - utop),
                  w_type * abs(type_costs_norm - utop)) +
      chebyshev_rho * (w_space * space_costs_norm +
                       w_time * time_costs_norm +
                       w_type * type_costs_norm)
  
    # cluster(s) with the minimum cost
    clust_i <- which(costs == min(costs, na.rm = TRUE))
    
    # if there is more than one candidate cluster
    if (length(clust_i) > 1) {
      # # pick one that has fewer points
      clust_n <- n[clust_i]
      clust_i <- clust_i[clust_n == min(clust_n)]
      if (length(clust_i) > 1) clust_i <- sample(clust_i, 1)
      clust[i] <- clust_i
    } else {
      # else just assign the cluster
      clust[i] <- clust_i
    }
  }
  
  return(clust)
}

# eval_sol --------------------------------------------------------------------
eval_sol <- function(data,
                     clust,
                     coords,
                     age,
                     type,
                     optim_type_diversity,
                     space_distmat,
                     hull_convex_ratio,
                     hull_crs,
                     sf_use_s2) {
  # total number of points
  N <- sum(!is.na(clust))
  # set of unique cluster indices
  K <- unique(stats::na.omit(clust))
  # total number of clusters in the output
  k_o <- length(K) # NA is excluded
  
  if (k_o < 2) {
    summary <- data.frame(k_o = k_o,
                          z_space = NA,
                          z_time1 = NA,
                          z_time2 = NA,
                          z_type = NA,
                          intersects = NA,
                          clustsize_f = NA)
    return(summary)
  }
  
  # initialise empty vectors for objective values and sizes of clusters
  # objective value
  z_space <- z_time1 <- z_time2 <- z_type <- numeric()
  # count flagged at the end
  clustsize_f <- logical()
  
  # extract the age column
  ages <- data[[age]]
  
  if (optim_type_diversity) {
    # extract the type column
    types <- data[[type]]
    # number of types in the data
    n_types <- length(unique(types))
  }
  
  # loop through k to obtain within cluster statistics
  for (j in K) {
    clust_points <- which(clust == j)
    n <- length(clust_points)
    clustsize_f <- n < N/k_o/2
    
    # space
    # sum of spatial distances to medoid for cluster j
    space_costs_j <- compute_space_costs(
      space_distmat = space_distmat,
      clust_points = clust_points)[clust_points]
    
    # space_wcd[j] <- sum(space_costs_j)
    z_space <- append(z_space, n * sum(space_costs_j)) # weighting n_j[j]
    
    # time
    # extract ages in cluster j
    ages_j <- ages[clust_points]
    
    # if type is concerned
    if (optim_type_diversity) {
      # extract types in cluster j
      types_j <- types[clust_points]
      
      # Shannon entropy for cluster j
      z_type <- append(z_type,  n * shannon_entropy(x = types_j, n = n_types))
      
      # temporal objectives
      for (q in unique(types)) {
        # subset the ages in clust j of type q
        clust_points_l <- clust_points[types_j == q]
        
        # skip if there is no point in the cluster
        if (length(clust_points_l) == 0) next
        
        ages_l <- ages[clust_points_l]
        n <- length(clust_points_l)
  
        # temporal objectives
        z_time1 <- append(z_time1, n * var(ages_l))
        z_time2 <- append(z_time2, n * stats::var(diff(sort(ages_l))))
      }
    } else {
      z_time1 <- append(z_time1, n * var(ages_j) )
      z_time2 <- append(z_time2, n * stats::var(diff(sort(ages_j))))
    }
  }
  
  # calculate summary statistics
  # https://stats.stackexchange.com/questions/122668/is-there-a-measure-of-evenness-of-spread
  # "n_qk / N" or "n_k / N" is the weight 
  
  # space_distance
  z_space <- sum(z_space, na.rm = TRUE) / N
  # time_variance
  z_time1 <- -1 * sum(z_time1, na.rm = TRUE) / N # multiply -1 for max objective 
  # time_evenness
  z_time2 <- -1 * sum(z_time2, na.rm = TRUE) / N # multiply -1 for max objective
  # type_diversity
  z_type <- -1 * sum(z_type) / N
  
  # evaluate if blobs are intersecting in space
  intersects <- check_intersects(data = data,
                                 clust = clust,
                                 coords = coords,
                                 hull_convex_ratio = hull_convex_ratio,
                                 hull_crs = hull_crs,
                                 sf_use_s2 = sf_use_s2)
  # count the number of clusters below threshold
  clustsize_f <- sum(clustsize_f)
  
  # return a data frame of all the statistics
  summary <- data.frame(k_o = k_o,
                        z_space = z_space,
                        z_time1 = z_time1,
                        z_time2 = z_time2,
                        z_type = z_type,
                        intersects = intersects,
                        clustsize_f = clustsize_f)
  return(summary)
}

# print.blobs -----------------------------------------------------------------
print.stblob_sol <- function(x, ...) {
  if(x$status == 0) {
    cat("STblob solution\n")
    cat("$clust\n# cluster assignments:\n")
    print(head(x$clust, 20))
    if (length(x$clust) > 20) cat(paste0(" <", length(x$clust) - 20," more points>\n"))
    cat("\n")
    cat("$summary\n# summary:\n")
    print(x$summary)
  } else {
    cat("No feasible solution was found.")
    cat("\n")
    cat(paste0("status code: ", x$status, "\n"))
  }
  invisible(x)
}

# summary.blobs ---------------------------------------------------------------
summary.stblob_sol <- function(sol, ...) {
  print(sol$summary)
  invisible(sol)
}

# helpers ---------------------------------------------------------------------
## remove_outliers ------------------------------------------------------------
# This version is based on the distances to medoid and remove them adhoc
remove_outliers <- function(data,
                            clust,
                            coords,
                            space_distmat = NULL,
                            space_distmethod = NULL,
                            iqrm = 1.5) {
  
  # total number of clusters
  K <- unique(stats::na.omit(clust)) # NA is excluded
  
  # check if is.null(space_distmat)
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod)
  
  # remove outliers
  space_costs <- numeric()
  order <- numeric()
  for (i in K) {
    clust_points <- which(clust == i)
    space_costs <- append(space_costs,
                          compute_space_costs(space_distmat = space_distmat,
                                              clust_points = clust_points)[clust_points])
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
  n_outliers <- length(outliers)
  if (n_outliers > 0) clust[outliers] <- NA
  
  out <- list(clust = clust, n_outliers = n_outliers)
  return(out)
}

## check_intersects -----------------------------------------------------------
check_intersects <- function(data,
                             clust,
                             coords,
                             hull_convex_ratio,
                             hull_crs,
                             sf_use_s2) { 
  # st_union breaks when s2 is on, switch it off until on.exit()
  old <- suppressMessages(sf::sf_use_s2(sf_use_s2))
  on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  
  data_sf <- sf::st_as_sf(data, coords = coords, crs = hull_crs)
  
  # obtain convex/concave hulls
  # If by_feature is TRUE each feature geometry is unioned individually.
  # This can for instance be used to resolve internal boundaries after
  # polygons were combined using st_combine.
  # https://r-spatial.github.io/sf/reference/geos_combine.html
  suppressMessages({
    hulls <- stats::aggregate(data_sf$geometry,
                              by = list(clust = clust),
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
  
  out <- if (any(intersects, na.rm = TRUE)) TRUE else FALSE
  return(out)
}

## select_bs ------------------------------------------------------------------
select_summary <- function(x) {
  # select stblob_lsearch() summary columns
  x <- x[ , c("k", "w_space", "w_time", "w_type",
              "k_o", "z_space", "z_time1", "z_time2", "z_type",
              "iter", "ari", "ls_convergence", "n_outliers")]
  rownames(x) <- NULL
  return(x)
}

select_trace <- function(x) {
  # select stblob_lsearch() summary columns
  x <- x[ , c("k", "w_space", "w_time", "w_type",
              "k_o",  "z_space", "z_time1", "z_time2", "z_type",
              "iter", "ari")]
  rownames(x) <- NULL
  return(x)
}

## new_blobs ------------------------------------------------------------------
new_sol <- function(clust, summary, trace, data, status, params) {
  stopifnot(is.numeric(clust) || (length(clust) == 1 && is.na(clust)),
            is.data.frame(summary) || (length(summary) == 1 && is.na(summary)),
            is.data.frame(trace) || (length(trace) == 1 && is.na(trace)),
            is.data.frame(data),
            is.numeric(status),
            is.list(params))
  
  structure(
    list(clust = clust, summary = summary, trace = trace, data = data,
         status = status, params = params),
    class = "stblob_sol"
  )
}