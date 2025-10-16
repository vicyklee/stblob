#------------------------------------------------------------------------------#
# Local search algorithm ####
#------------------------------------------------------------------------------#
#' Distance matrix computation
#' 
#' @description
#' This function computes and returns the distance matrix computed using the specified method
#' to compute the distances between the rows of a data matrix. Km is the unit when "geodesic" is specified.
#'
#' @param data a numeric matrix or data frame.
#' @param method a string of method. It must be one of "geodesic" or "euclidean". Default is "geodesic".
#'
#' @returns a numeric distance matrix. km is the unit when "geodesic" is specified.
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

#------------------------------------------------------------------------------#
#' Compute adjacency matrix using local scaling
#' 
#' @description
#' This function converts a distance matrix to a weighted adjacency matrix using local scaling
#' proposed by Zelnik-Manor and Perona (2004).
#'
#' @param distmat a distance matrix.
#' @param knn an integer of the k-th nearest neighbour. Default is 7.
#' 
#' @details
#' The method used to compute a weighted adjacency matrix was proposed by Zelnik-Manor and Perona (2004). 
#' Instead of using a fixed scaling parameter for all data points, a local scaling parameter is calculated for
#' every point based on their \eqn{k^th} nearest neighbour in the conversion of pairwise distance to pairwise adjacency (affinity/similarity).
#' See Zelnik-Manor and Perona (2004) for more information.
#' 
#' @returns a weighted adjacency matrix.
#' @references
#' Zelnik-Manor, L., & Perona, P. (2004). Self-Tuning Spectral Clustering. Neural Information Processing Systems, 17, 1601–1608.
#' 
#' @export

compute_adjacency <- function(distmat, knn = 7) {
  sorted_distmat <- t(apply(distmat, 1, sort)) # sort 
  knn_distmat <- sorted_distmat[, knn] # distance to the k-th nearest neighbour 
  V <- outer(knn_distmat, knn_distmat, "*") # a matrix of sigma_i * sigma_j
  W <- exp(-distmat^2 / V)
  diag(W) <- 0 # no self-loop
  return(W)
}

#------------------------------------------------------------------------------#
#' Compute graph Laplacian
#' 
#' @description
#' This function computes and returns a graph Laplacian.
#'
#' @param W a symmetric adjaceny matrix.
#' @param normalise a Boolean to normalise graph Laplacian. Default is TRUE.
#'
#' @returns a graph Laplacian matrix.
#' @references
#' Chung, F. R. (1997). Spectral graph theory (Vol. 92). American Mathematical Soc..
#' 
#' Smola, A. J., & Kondor, R. (2003, August). Kernels and regularization on graphs. In Learning theory and kernel machines: 16th annual conference on learning theory and 7th kernel workshop, COLT/kernel 2003, Washington, DC, USA, august 24-27, 2003. Proceedings (pp. 144-158). Berlin, Heidelberg: Springer Berlin Heidelberg.
#' 
#' @export

compute_laplacian <- function (W, normalise = TRUE) {
  d <- rowSums(W)
  if (normalise == FALSE) {
    L <- diag(d) - W
  } 
  if (normalise == TRUE) {
    D_inv_sqrt <- diag(1 / sqrt(d))
    D_inv_sqrt[!is.finite(D_inv_sqrt)] <- 0  # handle isolated nodes
    L <- diag(nrow(W)) - D_inv_sqrt %*% W %*% D_inv_sqrt
  }
  return(L)
}

#------------------------------------------------------------------------------#
#' Compute diffusion kernel
#' 
#' @description
#' This function computes and returns a kernel matrix using a diffusion kernel as defined in Kondor and Lafferty (2002).
#'
#' @param L a graph Laplacian.
#' @param beta a numeric parameter controlling the diffusion rate. 
#'
#' @returns a diffusion kernel matrix
#' @references 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' @export

diffusion_kernel <- function(L, beta) {
  eig <- eigen(L, symmetric = TRUE)
  U <- eig$vectors
  lambda <- eig$values
  
  exp_lambda <- exp(-beta * lambda)
  kmat <- U %*% diag(exp_lambda) %*% t(U)
  return(kmat)
}

#------------------------------------------------------------------------------#
#' Compute the effective rank from singular values of a matrix
#' 
#' @description
#' This function computes and returns the effective rank from singular values of a matrix following Roy and Vetterli (2007).
#'
#' @param sigma a numeric vector of singular values.
#' 
#' @details
#' Instead of treating \eqn{0\log{0} = 0} in Roy and Vetterli (2007), 100000000 is returned when any \eqn{\sigma_i = 0}
#' due to its application in the optimisation of \eqn{\beta} for diffusion kernel k-means clustering.
#' It avoids \eqn{\beta} that gives eigenvalues of 0 in the computation of the diffusion kernel. 
#'
#' @returns a numeric value of the effective rank.
#' @references 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' Roy, O., & Vetterli, M. (2007, September). The effective rank: A measure of effective dimensionality. In 2007 15th European signal processing conference (pp. 606-610). IEEE.

compute_erank_sigma <- function(sigma) {
  p_k <- sigma/sum(abs(sigma)) # sigma = 0 means redundancy
  if(any(p_k == 0)) return(1e8) 
  # p_k[p_k == 0] <- 1 # such that 0*log(0) = 0 as 1*log(1) = 0
  # the problem with treating 0*log(0) = 0 is that it doesn't optimise clustering given the graph.
  x_k <- p_k*log(p_k)
  H_k <- -sum(p_k*log(p_k))
  erank <- exp(H_k)
  return(erank)
}

#------------------------------------------------------------------------------#
#' Compute the effective rank of the diffusion kernel matrix
#' 
#' @description
#' This function computes and returns the effective rank of the diffusion kerenl matrix following Roy and Vetterli (2007).
#'
#' @inheritParams diffusion_kernel
#' 
#' @details
#' As the diffusion kernel is symmetric and positive definite,
#' its eigenvalues, defined as  \eqn{e^{\beta\lambda_i}},
#' are essentially singular values \eqn{\sigma_i}, where \eqn{\lamba_i} are eigenvalues of\eqn{L}.
#' 
#' Instead of treating \eqn{0\log{0} = 0}
#' when any singular values \eqn{\sigma_i = 0} in Roy and Vetterli (2007),
#' 1e8 is returned for its application in the optimisation of the diffusion rate (\eqn{\beta}) for kernel k-means clustering
#' with the constraint \eqn{e^{\beta\lambda_i}} > 0.
#'
#' @returns a numeric value of the effective rank.
#' @references 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' Roy, O., & Vetterli, M. (2007, September). The effective rank: A measure of effective dimensionality. In 2007 15th European signal processing conference (pp. 606-610). IEEE.
#' 
#' @export

compute_erank <- function(L, beta) {
  # PD matrix which kernel matrix is is symmetric and its eigenvalues are positive.
  # therefore eigenvalues are essentially the same as singular values.
  # we can just use exp(-beta*lambda) for optimising beta
  lambda <- eigen(L, symmetric = TRUE, only.values = T)$values
  sigma <- exp(-beta*lambda) 
  erank <- compute_erank_sigma(sigma)
  return(erank)
}

#------------------------------------------------------------------------------#
#' Loss function to minimise the difference in effective ranks
#' 
#' @description
#' This function calculates and returns the absolute difference between the targeted effective rank and
#' the effective rank of the diffusion kernel for a given `beta` (\eqn{\beta}) based on the effective rank defined by Roy and Vetterli (2007).
#' 
#' @inheritParams diffusion_kernel
#' @param k an integer of the targeted effective rank.
#' 
#' @details
#' As the diffusion kernel is symmetric and positive definite,
#' its eigenvalues, defined as  \eqn{e^{\beta\lambda_i}},
#' are essentially singular values \eqn{\sigma_i}, where \eqn{\lamba_i} are eigenvalues of\eqn{L}.
#' This is utilised to calculate the effective rank of the kernel for a given \eqn{\beta} following Roy and Vetterli (2007) and
#' this loss function to be minimised is defined as the absolute difference between that and the targeted effective rank.
#' 
#' Instead of treating \eqn{0\log{0} = 0}
#' when any singular values \eqn{\sigma_i = 0} in Roy and Vetterli (2007),
#' 1e8 is returned for its application in the optimisation of the diffusion rate (\eqn{\beta}) for kernel k-means clustering
#' with the constraint \eqn{e^{\beta\lambda_i}} > 0.
#' 
#' @returns a numeric value of the difference in effective ranks.
#' @references 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' Roy, O., & Vetterli, M. (2007, September). The effective rank: A measure of effective dimensionality. In 2007 15th European signal processing conference (pp. 606-610). IEEE.

min_erankdiff <- function(beta, L, k) {
  # PD matrix which kernel matrix is is symmetric and its eigenvalues are positive.
  # therefore eigenvalues are essentially the same as singular values.
  # we can just use exp(-beta*lambda) for optimising beta
  erank <- compute_erank(L, beta)
  diff <- abs(erank - k)
  return(diff)
}

#------------------------------------------------------------------------------#
#' Optimise beta for diffusion kernel
#' 
#' @description
#' This function optimises \eqn{\beta} for the diffusion kernel by minimising the absolute difference between
#' the number of clusters and the effective rank of the kernel matrix. 
#' 
#' @inheritParams diffusion_kernel
#' @param k an integer of the number of clusters.
#' @param par an numeric initial value for the parameter to be optimised over. The range is \eqn{[0.001, Inf]}. Default is 10.
#' 
#' @details
#' To perform diffusion kernel k-means clustering, the diffusion rate (\eqn{\beta}) has to be optimised to capture the global similarity of
#' a graph Lapclacian matrix \eqn{L} in the feature space.
#' 
#' This is done by minimising the absolute difference between the number of clusters and the effective rank of the kernel matrix where the effective rank,
#' which is a measure of effective dimensionality of a matrix, follows the definition in Roy and Vetterli (2007).
#' 
#' As the diffusion kernel is symmetric and positive definite,
#' its eigenvalues, defined as  \eqn{e^{\beta\lambda_i}},
#' are essentially singular values \eqn{\sigma_i}, where \eqn{\lamba_i} are eigenvalues of\eqn{L}.
#' This is utilised to calculate the effective rank of the kernel matrix for a given \eqn{\beta}
#' and the optimisation is performed by [stats::optim()] using "L-BFGS-B" method with a lower bound of 0.001.
#'  
#' @inherit stats::optim returns
#' @seealso [compute_erank()], [stats::optim()]
#' @references 
#' Roy, O., & Vetterli, M. (2007, September). The effective rank: A measure of effective dimensionality. In 2007 15th European signal processing conference (pp. 606-610). IEEE.
#' 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' @export

optim_beta <- function(L, k, par = 10) {
  out <- optim(par = par,
               fn = min_erankdiff,
               L = L,
               k = k,
               lower = 0.001,
               method = "L-BFGS-B",
               control = list(factr = 1e11) # less precise but does not affect the result (it actually converges at an optimised beta even though an error message on convergence is returned)
               )
  return(out)
}

#------------------------------------------------------------------------------#
#' Convert distance matrix to kernel matrix using diffusion kernel
#' 
#' @description
#' This function converts a distance matrix to a kernel matrix using a diffusion kernel.
#' 
#' @param distmat a distance matrix.
#' @param k an integer of the number of clusters.
#' @param w_knn an integer of the k-th nearest neighbour used in local scaling for computing the adjacency matrix, passed onto [compute_adjacency()]. Default is 7.
#' @param l_normalise a Boolean to normalise graph Laplacian, passed onto [compute_laplacian()]. Default is TRUE.
#' @param beta_par an numeric initial value for the parameter \eqn{\beta} to be optimised over for the diffusion kernel, passed onto [optim_beta()]. The range is \eqn{[0.001, Inf]}. Default is 10.
#' 
#' @details
#' This function takes a distance matrix as input and computes a weighted adjacency matrix using local scaling as proposed by Zelnik-Manor and Perona (2004).
#' A graph is constructed from the weighted adjacency matrix, represented by a graph Laplacian following Smola and Kondor (2003) and Chung (1997).
#' Then, it uses a diffusion kernel over the graph with optimised diffusion rate \eqn{\beta} to capture the global structure for kernel-based clustering (Kondor and Lafferty, 2002; Kondor and Vert, 2004).
#' \eqn{\beta} is optimised by minimising the absolute difference between the number of clusters and the effective rank of the kernel matrix where the effective rank,
#' which is a measure of effective dimensionality of a matrix, follows the definition in Roy and Vetterli (2007).
#' 
#' @returns a list of the following objects.
#' \itemize{
#'   \item \code{kmat}: a diffusion kernel matrix.
#'   \item \code{optim_out}: an output of [stats::optim()].
#' }
#' 
#' @seealso [compute_adjacency()], [compute_laplacian()], [diffusion_kernel()], [optim_beta()], [compute_erank()]
#' 
#' @references 
#' Chung, F. R. (1997). Spectral graph theory (Vol. 92). American Mathematical Soc..
#'
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' Roy, O., & Vetterli, M. (2007, September). The effective rank: A measure of effective dimensionality. In 2007 15th European signal processing conference (pp. 606-610). IEEE.
#' 
#' Smola, A. J., & Kondor, R. (2003, August). Kernels and regularization on graphs. In Learning theory and kernel machines: 16th annual conference on learning theory and 7th kernel workshop, COLT/kernel 2003, Washington, DC, USA, august 24-27, 2003. Proceedings (pp. 144-158). Berlin, Heidelberg: Springer Berlin Heidelberg.
#' 
#' @export

distmat_to_kmat <- function(distmat, k, w_knn = 7, l_normalise = TRUE, beta_par = 10) {
  # distance matrix to weighted adjaceny matrix
  W <- compute_adjacency(distmat = distmat, knn = w_knn)
  # compute normalised graph Laplacian
  L <- compute_laplacian(W = W, normalise = l_normalise)
  # optimise beta
  optim_out <- optim_beta(L = L, k = k, par = beta_par)
  beta <- optim_out$par
  # compute kernel matrix using diffusion kernel
  kmat <- diffusion_kernel(L, beta)
  return(list(kmat = kmat, optim_out = optim_out))
}

#------------------------------------------------------------------------------#
#' Compute kernel matrix using diffusion kernel from point data
#' 
#' @description
#' This function computes and returns a kernel matrix using a diffusion kernel taking the point data as input.
#' 
#' @inheritParams compute_distmat
#' @inheritParams distmat_to_kmat
#' 
#' @details 
#' This function is a wrapper of [compute_distmat()] and [distmat_to_kmat()]. Only geodesic or euclidean distance is available for computing the distances. See [distmat_to_kmat()] for more details.
#' 
#' @inherit distmat_to_kmat returns
#' 
#' @export

compute_kmat <- function(data, method = "geodesic", k, w_knn = 7, l_normalise = TRUE, beta_par = 10) {
  distmat <- compute_distmat(data = data, method = method)
  kmat_out <- distmat_to_kmat(distmat = distmat, k = k , w_knn = w_knn, l_normalise = l_normalise, beta_par = beta_par)
  return(kmat_out)
}

#------------------------------------------------------------------------------#
#' Spatial cost computation
#'
#' @description
#' This function computes and returns the spatial costs of all points for a given cluster.
#'
#' @param space_kmat a kernel matrix computed from the spatial distance matrix.
#' @param clust_points a numeric vector of point indices in the targeted cluster.
#' @param i an integer of a specific point index. Default is NULL.
#' 
#' @references
#' Dhillon, I. S., Guan, Y., & Kulis, B. (2004, August 22). Kernel k-means: spectral clustering and normalized cuts. Proceedings of the Tenth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. KDD04: ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, Seattle WA USA. https://doi.org/10.1145/1014052.1014118
#'
#' Schölkopf, B., Smola, A., & Müller, K.-R. (1998). Nonlinear component analysis as a kernel eigenvalue problem. Neural Computation, 10(5), 1299–1319.
#' 
#' @returns a numeric value or vector of spatial cost(s) for a given cluster.
#'
#' @export

compute_spacecost <- function(space_kmat, clust_points, i = NULL) {
  if (!is.null(i)) {
    k_ii <- space_kmat[i, i]
    k_ik <- space_kmat[i , clust_points, drop = FALSE]
  } else {
    k_ii <- diag(space_kmat)
    k_ik <- space_kmat[ , clust_points, drop = FALSE]
  }
  k_kk <- space_kmat[clust_points, clust_points, drop = FALSE]
  space_cost <- k_ii - 2*rowMeans(k_ik) + mean(k_kk)
  
  return(space_cost)
}

#------------------------------------------------------------------------------#
#' Intersect check 
#' 
#' @description
#' This function checks and returns a Boolean if any clusters are intersecting.
#'
#' @param data a numeric matrix or data frame.
#' @param clust a numeric vector of cluster assignment. Default is NULL when it is already present as a column in the data.
#' @param coords a vector of strings or numeric values indicating the columns of coordinates (longitude, latitide). Default is the first two columns.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()]. Default is 4326.
#' @param hull_convex_ratio a numeric value controlling the convexity of the hulls passed onto [sf::st_concave_hull()]. 1 returns convex and 0 maximally concave hulls. Default is 0.
#'
#' @returns TRUE or FALSE
#'
#' @export

intersects_bool <- function(data, clust = NULL, coords = c(1,2), crs = 4326, hull_convex_ratio = 0) {
  data <- as.data.frame(data)
  
  # if (is.null(crs)) { crs <- NA } # Euclidean
  if (!is.null(clust)) { data$clust <- clust }
  
  data_sf <- sf::st_as_sf(data, coords = coords, crs = crs)
  
  # obtain convex/concave hulls
  hulls <- stats::aggregate(data_sf$geometry, by = list(clust = data_sf$clust), sf::st_union) 
  hulls <- sf::st_as_sf(hulls) 
  hulls <- sf::st_concave_hull(hulls, ratio = hull_convex_ratio) # convex hulls
  
  # indicate TRUE/FALSE if there are any intersects
  intersects <- sf::st_intersects(hulls$geometry, sparse = F)
  diag(intersects) <- NA
  bool <- if (any(intersects == TRUE, na.rm = TRUE)) TRUE else FALSE
  
  return(bool)
}

#------------------------------------------------------------------------------#
#' Reorder clusters
#' 
#' @description
#' This function reorders and returns a vector of cluster assignment in ascending order.
#'
#' @param clust a numeric vector of cluster assignment.
#'
#' @returns a numeric vector of reordered cluster assignment.

reorder_clust <- function(clust) {
  # e.g. c(1,4,4,2) will become c(1,2,2,3)
  
  if(length(unique(clust)) == 1) return(rep(1, length(clust)))
  
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
    sets <- as.matrix(stats::na.omit(sets))
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
#'
#' @inheritParams intersects_bool
#' @inheritParams compute_spacecost
#' @param age a string or numeric value indicating the column of age. Default is the third column. 
#' @param space_distmat a spatial distance matrix used when `space_kmat` is not supplied. Default is NULL.
#' @param space_distmethod a string of method used when `space_kmat` and `space_distmat` are not specified. It must be one of "geodesic" or "euclidean". Default is NULL.
#' @param w_knn an integer of the k-th nearest neighbour used when `space_kmat` used in local scaling for computing the adjacency matrix, passed onto [compute_adjacency()]. Default is 7 when `space_kmat` is not supplied.
#' @param l_normalise a Boolean to normalise graph Laplacian, passed onto [compute_laplacian()]. Default is TRUE when `space_kmat` is not supplied.
#' @param beta_par an numeric initial value for the parameter \eqn{\beta} to be optimised over for the diffusion kernel, passed onto [optim_beta()]. The range is \eqn{[0.001, Inf]}. Default is 10 when `space_kmat` is not supplied.
#' 
#' @details
#' The critical size of a cluster is defined as \eqn{\frac{N}{2K}} where \eqn{N} is the number of data points and \eqn{k} is the number of clusters.
#'
#' @returns a list of the following objects.
#' \itemize{
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{clust_below_size}: a numeric vector of clusters below the critical size.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} when `space_kmat` is not supplied.
#' }
#'
#' @export

eval_blobs <- function(data,
                       coords = c(1,2),
                       age = 3,
                       crs = 4326,
                       hull_convex_ratio = 0,
                       space_kmat = NULL,
                       space_distmat = NULL,
                       space_distmethod = NULL,
                       w_knn = NULL,
                       l_normalise = NULL,
                       beta_par = NULL) {
  
  # total number of points
  N <- nrow(data)
  # total number of clusters
  k <- length(unique(stats::na.omit(data$clust))) # NA is excluded
  # initialise empty vectors 
  space_ss <- time_range <- time_evenness <- n <- numeric(k)
  
  # compute space_kmat
  if (is.null(space_kmat)) {
    if (is.null(w_knn)) w_knn <- 7
    if (is.null(l_normalise)) l_normalise <- TRUE
    if (is.null(beta_par)) beta_par <- 10
    
    if(is.null(space_distmat)) {
      if (is.null(space_distmethod)) {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
        message(paste0(space_distmethod," is used to compute space_distmat"))
      } else {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      }
      space_kmat_out <- compute_kmat(data = data[, coords],
                                     method = space_distmethod,
                                     k = k,
                                     w_knn = w_knn,
                                     l_normalise = l_normalise,
                                     beta_par = beta_par)
    } else {
      space_kmat_out <- distmat_to_kmat(distmat = space_distmat,
                                        k = k,
                                        w_knn = w_knn,
                                        l_normalise = l_normalise,
                                        beta_par = beta_par)
    }
    space_kmat <- space_kmat_out$kmat
    space_kmat_optim_out <- space_kmat_out$optim_out
  } else {
    space_kmat_optim_out <- NULL # As it is one of the returned items
  }
  
  # loop over k to obtain within cluster statistics
  for (j in 1:k) {
    clust_points <- which(data$clust == j)
    if (length(clust_points) == 0) next
    # data_k <- subset(data, data$clust == j)
    
    # spatial objective
    # kernel k means cost function
    space_ss[j] <- sum(compute_spacecost(space_kmat = space_kmat, clust_points = clust_points, i = NULL)[clust_points])
    # temporal objectives
    time_range[j] <- max(data_k[ , age], na.rm = T) - min(data_k[ ,age], na.rm = T)
    time_evenness[j] <- 1 / (1 + stats::var(diff(sort(data_k[ , age])))) # NA if there are fewer than 3 data points
    # other constraint related statistics
    n[j] <- length(clust_points)
  }
  
  # Calculate the summary statistics
  space_wcss <- sum(space_ss)
  time_range_mean <- mean(time_range, na.rm = TRUE)
  time_range_sd <- stats::sd(time_range, na.rm = TRUE)
  time_evenness_mean <- mean(time_evenness, na.rm = TRUE) # na.rm = TRUE to acess the overall performance as we normally do not care clusters with fewer than 3 points
  time_evenness_sd <- stats::sd(time_evenness, na.rm = TRUE) # na.rm = TRUE to acess the overall performance as we normally do not care clusters with fewer than 3 points
  size_diff <- max(n, na.rm = TRUE) - min(n, na.rm = TRUE)
  size_mean <- mean(n, na.rm = TRUE)
  size_sd <- stats::sd(n, na.rm = TRUE)
  
  # evaluate failed clust
  # n.points must be > N / k / 2
  clust_below_size <- which(n < N/k/2)
  n_fail <- sum(n[clust_below_size])
  
  # evaluate if blobs are intersecting in space
  intersects <- intersects_bool(data = data, coords = coords, crs = crs, hull_convex_ratio = hull_convex_ratio)
  
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
  
  return(list(summary = summary, clust_below_size = clust_below_size, space_kmat_optim_out = space_kmat_optim_out))
}

#------------------------------------------------------------------------------#
#' Assign the starting cluster members
#' 
#' @description
#' Assign the starting cluster members by approximating the maximum spatial separation between points or by random assignment.
#'
#' @param data a data matrix or data frame.
#' @param k an integer of the number of clusters.
#' @param random_start a Boolean to use random starting cluster members. Default is FALSE.
#' @inheritParams compute_spacecost
#'
#' @returns a data frame with assigned starting clusters as a column.
#' @export

start_blobs <- function(data, k, space_kmat, random_start = FALSE) {
  
  if (random_start == TRUE) {
    data <- as.data.frame(data)
    data$clust <- NA
    start <- sample(1:nrow(data), k)
    data$clust[start] <- 1:k
    data$order <- 1:nrow(data)
    return(data)
  }
  
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
  max_similarity <- numeric(N)
  for(n in 1:N){
    i <- sort(sample(1:nrow(data), size = k)) # sample k points from the data, sort them 
    mat[n, ] <- i # store sorted index in the matrix by row
    cb <- utils::combn(k, 2) # all combinations of k chooses 2 by column
    NC <- ncol(cb) # NC number of combinations
    similarities <- numeric(NC) # a vector of distances of NC long
    for(c in 1:NC) similarities[c] <- space_kmat[ i[cb[1, c]], i[cb[2, c]] ] # extract from distance matrix the distances for all combinations of points
    max_similarity[n] <- max(similarities) # find the maximum hence minimum separation
  }
  # pick the set with largest minimum distance between two points to ensure maximum separation between k points;
  # pick the first one if two are tied 
  start <- mat[which(max_similarity == min(max_similarity))[1], ]
  # Assign cluster memberships to the starting points
  data[start, "clust"] <- 1:k
  data <- as.data.frame(data)
  
  return(data)
}

#------------------------------------------------------------------------------#
#' Assign clusters
#' 
#' @description
#' This function assigns clusters for a given k and r.
#' 
#' @param data a data matrix or data frame.
#' @param k an integer of the number of clusters.
#' @param r an numeric value of spatial relative weight. It must be \eqn{[0,1]}.
#' @param age a string or numeric value indicating the column of age. Default is the third column.
#' @inheritParams compute_spacecost
#'
#' @details
#' A diffusion kernel is applied to compute the distance to the centroid in space.
#' 
#' @seealso [compute_spacecost()]
#' 
#' @references
#' Dhillon, I. S., Guan, Y., & Kulis, B. (2004, August 22). Kernel k-means: spectral clustering and normalized cuts. Proceedings of the Tenth ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. KDD04: ACM SIGKDD International Conference on Knowledge Discovery and Data Mining, Seattle WA USA. https://doi.org/10.1145/1014052.1014118
#' 
#' Kondor, R. I., & Lafferty, J. (2002, July). Diffusion kernels on graphs and other discrete structures. In Proceedings of the 19th international conference on machine learning (Vol. 2002, pp. 315-322).
#' 
#' Kondor, R., & Vert, J. P. (2004). Diffusion kernels. kernel methods in computational biology, 171-192.
#' 
#' Schölkopf, B., Smola, A., & Müller, K.-R. (1998). Nonlinear component analysis as a kernel eigenvalue problem. Neural Computation, 10(5), 1299–1319.
#'
#' @returns a data frame with assigned clusters as a column.
#' @export

find_blobs <- function(data, k, r, space_kmat, age = 3) {
  if (is.null(data$order)) stop("Did you forget to run start_blobs()?")
  
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
  
  # Reorder kmat to match the reordered data
  order <- data[ ,"order"]
  space_kmat <- space_kmat[order, order]
  
  # Extract clust to make the code cleaner
  clust <- data[,"clust"]
  
  # initialise stat matrices for the unassigned points
  start <- if (length(a_points) > nrow(data)) 1 else length(a_points) + 1
  N <- length(start:nrow(data))
  space_cost <- time_cost <- n <- numeric(k)
  
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
      space_cost[j] <- compute_spacecost(space_kmat = space_kmat,
                                         clust_points = clust_points,
                                         i = i)
      
      # compute the temporal cost
      clust_points_tmp <- if (i %in% clust_points) clust_points[clust_points != i] else clust_points
      
      if (length(clust_points_tmp) == 0) {
        time_cost[j] <- 0  # only i was in the cluster
      } else {
        time_cost[j] <- min(abs(data[i, age] - data[clust_points_tmp, 3]))
      }
      
    }
    
    # normalise the cost
    space_cost_norm <- (space_cost - min(space_cost, na.rm = T)) / (max(space_cost, na.rm = T) - min(space_cost, na.rm = T))
    time_cost_norm <- (time_cost - min(time_cost, na.rm = T)) / (max(time_cost, na.rm = T) - min(time_cost, na.rm = T))
    space_cost_norm[is.na(space_cost_norm)] <- 0
    time_cost_norm[is.na(time_cost_norm)] <- 0
    
    # weighted scalarising (into a single cost)
    cost <- space_cost_norm*r - time_cost_norm*(1-r)
    
    # assign cluster
    # check if there are tied clusters
    clust_tmp <- which(cost == min(cost, na.rm = T))
    
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
#' Core local-search algorithm
#' 
#' @description
#' This function performs the core bi-objective optimisation algorithm to assign clusters for a given k and r in an iterative fashion.
#' 
#' @inheritParams start_blobs
#' @inheritParams find_blobs
#' @inheritParams eval_blobs
#' @param iter an integer of the number of iterations. Default is 3.
#' @param converge_ari a numeric value of Adjusted Rand Index (ARI) that sets convergence threshold between two searches. It must be \eqn{[0,1]}. Default is NULL.
#' @param filter_intersects a Boolean to remove an assignment with intersects in space? Default is TRUE.
#' @param filter_clustsize a Boolean to assign NA to clustes below the critical size. Default is TRUE.
#' @param max_na a numeric value of the maximum proportion of NAs allowed. It must be \eqn{[0,1]}. Default is 0.05.
#'
#' @details
#' A diffusion kernel is applied to compute the distance to the centroid in space.]
#' See [distmat_to_kmat()] for more details for converting a distance matrix to a kernel matrix.
#'
#' Clusters are assigned in every iteration. It iterates until the set length or convergence. 
#'
#' When `converge_ari` is specified, convergence is defined and activated when ARI between the latest and the previous search is
#' above the specified threshold and at least three iterations are run.
#'
#' The critical size of a cluster is defined as \eqn{\frac{N}{2k}} where \eqn{N} is the number of data point and \eqn{k} is the number of clusters.
#'
#' @returns 
#' a numeric value of 1 if the search result is removed due to the presence of intersects,
#' 2 if removed due to the proportion of NAs exceeding `max_na`,
#' 3 if removed due to 1 resulting cluster,
#' or a list of the following objects if the search is meets the constraints,
#' \itemize{
#'   \item \code{data}: a data frame of the input data with assigned clusters as a column.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{clust_below_size}: a numeric vector of clusters below the critical size.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} in [distmat_to_kmat()] when `space_kmat` is not supplied.
#' }
#'
#' @seealso [distmat_to_kmat()] and [sf::st_as_sf()]

blob_search <- function(data,
                        k,
                        r,
                        iter = 3,
                        converge_ari = NULL,
                        coords = c(1,2),
                        age = 3,
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
                        beta_par = NULL) {
  
  #----------------------------------------------------------------------------#
  # spatial kernal matrix
  #----------------------------------------------------------------------------#
  # compute space_kmat
  if (is.null(space_kmat)) {
    if (is.null(w_knn)) w_knn <- 7
    if (is.null(l_normalise)) l_normalise <- TRUE
    if (is.null(beta_par)) beta_par <- 10
    
    if(is.null(space_distmat)) {
      if (is.null(space_distmethod)) {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
        message(paste0(space_distmethod," is used to compute space_distmat"))
      } else {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      }
      space_kmat_out <- compute_kmat(data = data[, coords],
                                     method = space_distmethod,
                                     k = k,
                                     w_knn = w_knn,
                                     l_normalise = l_normalise,
                                     beta_par = beta_par)
    } else {
      space_kmat_out <- distmat_to_kmat(distmat = space_distmat,
                                        k = k,
                                        w_knn = w_knn,
                                        l_normalise = l_normalise,
                                        beta_par = beta_par)
    }
    space_kmat <- space_kmat_out$kmat
    space_kmat_optim_out <- space_kmat_out$optim_out
  } else {
    space_kmat_optim_out <- NULL # As it is one of the returned items
  }
  
  #----------------------------------------------------------------------------#
  # search algorithm
  #----------------------------------------------------------------------------#
  # start_blobs() to pick centroids
  data <- start_blobs(data = data, k = k, space_kmat = space_kmat, random_start = random_start)
  
  # initialise counter counting find_blobs()
  t <- 0
  
  # intialise trace table
  trace <- data.frame()
  
  for (i in 1:iter) {
    
    data_old <- data
    # find_blobs()
    data <- find_blobs(data = data, k = k, r = r, age = age, space_kmat = space_kmat)
    # count find_blobs() executed
    t <- t + 1
    
    if (t > 0) {
      # check convergence
      ari <- mclust::adjustedRandIndex(data$clust, data_old$clust)
      
      # eval_blobs()
      eval_out <- eval_blobs(data,
                             coords = coords,
                             age = age,
                             crs = crs,
                             hull_convex_ratio = hull_convex_ratio,
                             space_kmat = space_kmat)
      clust_below_size <- eval_out$clust_below_size
      trace_newrow <- eval_out$summary
      trace_newrow$iter <- t
      trace_newrow$ari <- ari
      trace <- rbind(trace, trace_newrow)
      
      # if converged between t and t-1, break
      if (!is.null(converge_ari)) {
        if (all(ari >= converge_ari) == TRUE & t >= 3) break
      }
    }
  }
  
  trace$r <- r
  summary <- trace[nrow(trace), ]
  
  # order is not useful in the result, remove the column
  data$order <- NULL
  
  #----------------------------------------------------------------------------#
  # filters / constraints
  #----------------------------------------------------------------------------#
  N <- nrow(data)
  
  # 1. intersects
  if (filter_intersects == T) {
    intersects <- summary$intersects
    if (intersects == T) {
      # message("At least two clusters intersect. Return NULL.")  
      return(1)
    }
  }
  
  # 2. cluster size
  # initialise n_removed column. If no point is removed, the entries will all be 0
  summary$n_removed <- 0
  
  if (filter_clustsize == T) {
    if (length(clust_below_size) > 0) {
      # assign NA to removed clusters
      n_removed <- summary$n_fail
      data$clust[which(data$clust %in% clust_below_size)] <- NA
      # reassign the cluster number
      # make sure there is no gap in the sequence of k
      data$clust <- reorder_clust(data$clust)
      # update eval_out$summary
      eval_out_updated <- eval_blobs(data,
                                     coords = coords,
                                     age = age,
                                     crs = crs,
                                     hull_convex_ratio = hull_convex_ratio,
                                     space_kmat = space_kmat)
      updated_cols <- intersect(names(summary), names(eval_out_updated$summary))
      summary[ , updated_cols] <- eval_out_updated$summary[ , updated_cols]
      summary$n_removed <- n_removed
    }
  } else {
    # reassign the cluster number
    # make sure there is no gap in the sequence of k
    data$clust <- reorder_clust(data$clust)
  }
  
  # 3. filter too many NAs and k < 2
  # return NULL if too many points are removed
  if (summary$n_removed > N * max_na) return(2)
  if (summary$k < 2) return(3)
  
  # 4. remove some columns
  summary$n_fail <- NULL # confusing to have both n_removed and this in the output
  trace$n_fail <- NULL # confusing to have both n_removed and this in the output
  
  # 5. remove clust.below.size as filter has been applied so it is irrelevant to the downstream
  clust_below_size <- NULL
  
  #-------------------------------------------------------------------#
  # polishing the output
  #-------------------------------------------------------------------#
  summary <- summary[, c("k", "r",
                         "space_wcss",
                         "time_range_mean", "time_range_sd",
                         "time_evenness_mean", "time_evenness_sd",
                         "size_mean", "size_sd", "size_diff",
                         "intersects", "n_removed",
                         "iter", "ari")]

  trace <- trace[, c("k", "r",
                     "space_wcss",
                     "time_range_mean", "time_range_sd",
                     "time_evenness_mean", "time_evenness_sd",
                     "size_mean", "size_sd", "size_diff",
                     "intersects",
                     "iter", "ari")]
  
  rownames(summary) <- NULL
  rownames(trace) <- NULL
  
  return(list(data = data, summary = summary, trace = trace,
              clust_below_size = clust_below_size, space_kmat_optim_out = space_kmat_optim_out))
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
#' @param ari a numeric value of Adjusted Rand Index (ARI) that sets duplication threshold between two assignments. It must be \eqn{[0,1]}. Default is 1. See also [mclust::adjustedRandIndex()].
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
  
  # number of solutions
  n_sol <- nrow(clust)
  # get all combination of pairs
  pairs <- utils::combn(n_sol, 2)
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
#' Convert a list of blob objects to a pop object
#' 
#' @description
#' This function convert a list of blob objects to a pop object with combined output
#' from [blob_search()] and counts of filtered solutions by category.
#' 
#' @param blob_list a list of output from [blob_search()]
#' 
#' @returns a `pop` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of cluster assignments. Each row is a solution.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#'   \item \code{n_filtered}: a data frame of numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: a list of [stats::optim()] output from the optimisation of \eqn{\beta} in [distmat_to_kmat()].
#' }

convert_to_pop <- function(blob_list) {
  pop <- blob_list
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
  
  # extract clust
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
    # index the summary and trace
    summary$idx <- 1:nrow(summary)
    trace <- merge(trace, summary[, c("idx","run")], by = "run", all.x = TRUE)
  }
  
  n_filtered <- data.frame(intersects = filtered_intersects, size = filtered_size, k1 = filtered_k1, dup = filtered_dup)
  # n_filtered <- c(filtered_intersects, filtered_size, filtered_k1, filtered_dup)
  # names(n_filtered) <- c("intersects", "size", "k1", "dup")
  
  #-------------------------------------------------------------------#
  # polishing the output
  #-------------------------------------------------------------------#
  summary <- summary[, c("idx",
                         "k", "r", "run",
                         "space_wcss",
                         "time_range_mean", "time_range_sd",
                         "time_evenness_mean", "time_evenness_sd",
                         "size_mean", "size_sd", "size_diff",
                         "intersects", "n_removed",
                         "iter", "ari", "dup")]
  
  trace <- trace[, c("idx",
                     "k", "r", "run",
                     "space_wcss",
                     "time_range_mean", "time_range_sd",
                     "time_evenness_mean", "time_evenness_sd",
                     "size_mean", "size_sd", "size_diff",
                     "intersects",
                     "iter", "ari")]
  
  rownames(summary) <- NULL
  rownames(trace) <- NULL

  pop <- list(clust = clust, summary = summary, trace = trace, n_filtered = n_filtered)
  return(pop)
}

#------------------------------------------------------------------------------#
#' Populate solutions by weighted sum scalarisation
#' 
#' @description
#' This function populates solutions by weighted sum scalarisation of the bi-objective function in [blob_search()] for a given k. 
#' 
#' @inheritParams start_blobs
#' @inheritParams find_blobs
#' @inheritParams blob_search
#' @param k an integer value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the number of clusters.
#' @param r a numeric value or vector of length 2L. If a vector is supplied, they specify the lower and upper bounds of the relative spatial weight. They must be \eqn{[0,1]}. Default is c(0.5,1).
#' @param run an integer of the number of runs.
#' 
#' @details
#' A diffusion kernel is applied to compute the distance to the centroid in space.]
#' See [distmat_to_kmat()] for more details for converting a distance matrix to a kernel matrix.
#'
#' Clusters are assigned in every iteration. It iterates until the set length or convergence. 
#'
#' When `converge_ari` is specified, convergence is defined and activated when ARI between the latest and the previous search is
#' above the specified threshold and at least three iterations are run.
#'
#' The critical size of a cluster is defined as \eqn{\frac{N}{2k}} where \eqn{N} is the number of data point and \eqn{k} is the number of clusters.
#' 
#' Scalarisation is achieved by varying the relative spatial weight generated by Latin hypercube sampling using [lhs::randomLHS()].
#' 
#' To parallelise runs, [future.apply::future_lapply()] is implemented. See [future::future] and [future.apply::future.apply] for more information.
#' 
#' `k` is set to \eqn{[2, \lfloor\sqrt{N}\rfloor]} when NULL.
#' 
#' @returns a `pop` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of cluster assignments. Each row is a solution.
#'   \item \code{summary}: a data frame of summary statistics.
#'   \item \code{trace}: a data frame of summary statistics for tracing.
#'   \item \code{n_filtered}: a data frame of numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} in [distmat_to_kmat()] when `space_kmat` is not supplied.
#' }
#' 
#' @seealso [distmat_to_kmat()], [sf::st_as_sf()], [lhs::randomLHS()], [mclust::adjustedRandIndex()], [future::future], [future.apply::future.apply]
#' @export

blob_populate <- function(data,
                          k,
                          r = c(0.5,1),
                          iter = 3,
                          run,
                          converge_ari = NULL,
                          coords = c(1,2),
                          age = 3,
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
                          beta_par = NULL) {
  
  # compute space_kmat
  if (is.null(space_kmat)) {
    if (is.null(w_knn)) w_knn <- 7
    if (is.null(l_normalise)) l_normalise <- TRUE
    if (is.null(beta_par)) beta_par <- 10
    
    if(is.null(space_distmat)) {
      if (is.null(space_distmethod)) {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
        message(paste0(space_distmethod," is used to compute space_distmat"))
      } else {
        space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      }
      space_kmat_out <- compute_kmat(data = data[, coords],
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
  
  # sample r
  if (length(r) == 2) {
    # LHS sampling for more evenly distributed parameters
    lhs_samples <- lhs::randomLHS(run, 1)
    # scale to the range
    r_samples <- sort(as.vector(min(r) + lhs_samples * (max(r) - min(r))))
  } else {
    r_samples <- rep(r, run)
  }
  
  # assign k's lb and ub when NULL
  N <- nrow(data)
  if (is.null(k)) {
    k <- integer(2L)
    k[1] <- 2L # lower bound
    k[2] <- as.integer(floor(sqrt(N))) # upper bound, an heuristic to maximise information e.g. 100 points: 10 blobs, 10 points each 
  }
  
  # for single k
  if (length(k) == 1) {
    # run blob_search() in parallel
    blob_list <- future.apply::future_lapply(r_samples, function (r) {
      blob_search(data = data,
                  k = k,
                  r = r,
                  iter = iter,
                  converge_ari = converge_ari,
                  coords = coords,
                  age = age,
                  crs = crs,
                  hull_convex_ratio = hull_convex_ratio,
                  random_start = random_start,
                  filter_intersects = filter_intersects,
                  filter_clustsize = filter_clustsize,
                  max_na = max_na,
                  space_kmat = space_kmat) 
    }, future.seed = T)
    
    pop <- convert_to_pop(blob_list)
    pop$space_kmat_optim_out <- space_kmat_optim_out
    
  } else {
    # for a range of k
    if (length(k) == 2) {
      # create a grid of all combinations
      k_vec <- min(k):max(k)
      r_vec <- r_samples
      grid <- expand.grid(k_vec = k_vec, r_vec = r_vec)
      
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
                             crs = crs,
                             hull_convex_ratio = hull_convex_ratio,
                             random_start = random_start,
                             filter_intersects = filter_intersects,
                             filter_clustsize = filter_clustsize,
                             max_na = max_na,
                             space_kmat = space_kmat),
          k = k
        )
      }, grid$k_vec, grid$r_vec, future.seed = TRUE)
      
      # group runs of same k
      pop_list <- vector("list", k[2])
      for(i in 1:length(blob_list)) {
        m <- blob_list[[i]][["k"]]
        # append to pop_list
        pop_list[[m]] <- append(pop_list[[m]], list(blob_list[[i]][["blob"]]))
      }
      
      pop_list <- pop_list[-1] # as the loop starts at 2
      # convert to pop object for each k group
      pop_list <- lapply(pop_list, convert_to_pop) 
      
      # extract each element and add a column to indicate the initial k
      pop <- do.call(rbind, pop_list)
      
      summary_list <- lapply(seq_along(pop[ , "summary"]), function (i) cbind(pop[ , "summary"][[i]], k_o = i + 1))
      trace_list <- lapply(seq_along(pop[ , "trace"]), function (i) cbind(pop[ , "trace"][[i]], k_o = i + 1)) 
      clust_list <- pop[ ,"clust"]
      n_filtered_list <- pop[ , "n_filtered"]
      
      # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
      clust <- do.call(rbind, clust_list)
      summary <- do.call(rbind, summary_list)
      trace <- do.call(rbind, trace_list)
      n_filtered <- do.call(rbind, n_filtered_list)
      n_filtered <- cbind(k_o = 2:(nrow(n_filtered) + 1), n_filtered)
      
      # reindex the solutions
      summary$idx <- NULL
      summary$idx <- 1:nrow(summary)
      trace$idx <- NULL
      trace <- merge(trace, summary[, c("idx","k_o","run")], by = c("k_o","run"), all.x = TRUE)
      
      summary <- summary[, c("idx",
                             "k_o", "k", "r", "run",
                             "space_wcss",
                             "time_range_mean", "time_range_sd",
                             "time_evenness_mean", "time_evenness_sd",
                             "size_mean", "size_sd", "size_diff",
                             "intersects", "n_removed",
                             "iter", "ari", "dup")]

      trace <- trace[, c("idx",
                         "k_o", "k", "r", "run",
                         "space_wcss",
                         "time_range_mean", "time_range_sd",
                         "time_evenness_mean", "time_evenness_sd",
                         "size_mean", "size_sd", "size_diff",
                         "intersects",
                         "iter", "ari")]

      rownames(summary) <- NULL
      rownames(trace) <- NULL
      
      pop <- list(clust = clust, summary = summary, trace = trace, n_filtered = n_filtered, space_kmat_optim_out = space_kmat_optim_out)
    }
  }
  return(pop)
}

