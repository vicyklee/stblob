#' STblob (SpatioTemporal blob) optimisation algorithm
#' 
#' @description
#' The function performs multi-objective optimisation (MOO)
#' for a range of cluster numbers
#' and relative spatial weights
#' to approximate the Pareto front
#' of spatiotemporal blobs for regional time series analyses.
#' 
#' @inheritParams stblob_lsearch
#' @param k number of clusters. Either an integer or integer sequence.
#' @param w_space relative spatial weight of range \eqn{[0,1]}.
#' Either a numeric value or a vector of length 2.
#' The vector should specify the lower and upper bounds.
#' @param optim_type_diversity a logical value to optimise data type diversity.
#' Default is `FALSE`.
#' @param w_time relative spatial weight of range \eqn{[0,1]}. Default is `NULL`,
#' Supply a numeric value when fixed weights across runs are used.
#' @param w_time relative spatial weight of range \eqn{[0,1]}. Default is `NULL`,
#' Supply a numeric value when fixed weights across runs are used.
#' @param run number of runs per \eqn{k_i}. Default is `100L`.
#' @param batch number of batches. Default is `50L`.
#' @param parallel a logical value to parallelise the search. Default is `TRUE`.
#' @param workers number of parallel workers passed onto [future::plan()].
#' Default is `NULL` with [future::availableCores()] ` - 1`.
#' @param adapt_search a logical value to use adaptive search. `k` and
#' the lower bound of `w_space` are updated to optimise search space. Default is `TRUE`.
#' @param moo_tol tolerance for MOO convergence. The minimum proportional change
#' in hypervolume for three consecutive batches.
#' @inheritDotParams stblob_lsearch
#' 
#' @details
#' STblob (SpatioTemporal blob) optimisation algorithm performs multi-objective
#' optimisation (MOO) powered by a local search algorithm
#' (see [stblob_lsearch()]) in a MOO framework.
#' 
#' By sampling a range of parameter combinations of `w_space`, `w_time` and `w_type`,
#' it populates solutions that locally optimise trade-offs between spatial,
#' temporal and data type cost functions using augmented weighted
#' Chebyshev (Tchebycheff) scalarisation (Steuer & Choo, 1983)
#' across \eqn{k} clusters.
#' Sampling is done via Latin hypercube sampling ([lhs::randomLHS()]) for a
#' single free parameter, otherwise from a Dirichlet distribution
#' with \eqn{\alpha = (1,1,1)} ([MCMCpack::rdirichlet()]).
#' If `w_space` is specified by a range of values, input values of
#' `w_time` and `w_type` will be ignored.
#' 
#' To approximate the global optimal trade-offs,
#' non-dominated solutions are sorted across cumulative batches
#' to assess the advance of the Pareto front based on the following objectives:
#' 1. Average within-cluster spatial distances to the medoid
#' 2. Average within-cluster temporal variance
#' 3. Average within-cluster temporal variance of first differences
#' 4. Average within-cluster Shannon entropy of data types (`optim_type_diversity = TRUE`)
#' 
#' Objective values (\eqn{z}) are returned in `summary`
#' by the terms `z_space`, `z_time1`, `z_time2` and `z_type` respectively.
#' The objectives to maximise are multiplied by -1 as objective functions
#' are minimised.
#' 
#' Three commonly used quality indicators,
#' the inverted generational distance (IGD) (Coello Coello & Reyes Sierra, 2004),
#' IGD+ (Ishibuchi et al., 2015)
#' and hypervolume (Zitzler & Thiele, 1998)
#' (see also [moocore::igd()], [moocore::igd_plus()] and [moocore::hypervolume()]),
#' are used as diagnostics of the advance of Pareto front across successive batches.
#' For IGD and IGD+, the non-dominated set of all solutions is used as the
#' reference whereas a vector of \eqn{(1.1,1.1,1.1,...)} for \eqn{m} objectives
#' is used as the reference point for computing hypervolumes.
#' Note that objective values are normalised using the non-dominated set of all
#' solutions.
#' 
#' `adapt_search = TRUE` does two things to accelerate the search.
#' First, it removes converged \eqn{k_i} from `k` in the subsequent batches. 
#' Second, it updates the lower bound of `w_space` with the upper quartile (Q3)
#' of `w_space` across solutions with intersecting boundaries for each \eqn{k_i}.
#' 
#' `outlier_removal` is being tested at the moment. The idea is to remove points
#' whose distance to their medoid is considered an outlier based on the
#' interquantile range (IQR) of distances to medoid across \eqn{k} clusters.
#' 
#' @returns
#' an S3 object of class `stblob_moo` with the following components:
#'  * `clust_pf`: a matrix of cluster assignments on the Pareto front,
#'  ordered by `z_space`, `z_time1`, `z_time2` and `z_type`.
#'  Row: solution; column: data point.
#'  * `clust`: a matrix of all cluster assignments.
#'  Row: solution; column: data point.
#'  * `summary`: a data frame of summary statistics. They include the parameter
#'  value of (`k`), number of clusters of the output (`k_o`), four
#'  objective values (`z_space`, `z_time1`, `z_time2` and `z_type`),
#'  number of iterations (`iter`), ARI with the previous iteration (`ari`),
#'  state of local search convergence (`ls_convergence`)
#'  and number of outliers (`n_outliers`).
#'  `idx_pf` and `idx` index the solutions corresponding to the rows in
#'  `clust_pf` and `clust` respectively.
#'  `run` and `batch` indicate the run and batch from which the solution was generated.
#'  * `trace`: a data frame of summary statistics for each iteration.
#'  * `moo_convergence`: a logical value to indicate MOO convergence.
#'  * `moo_diagnostics`: a data frame of MOO quality
#'  indicators for every cumulative batch.
#'  * `data`: a data frame of the input data.
#'  * `space_distmat`: a spatial distance matrix used for clustering.
#'  * `filter_summary`: a data frame of filter summary.
#'  * `adapt_search_params`: a list of parameter values for `adapt_search`.
#'  * `params`: a list of parameter values.
#' 
#' @seealso [compute_distmat()], [geosphere::distGeo()], [stblob_lsearch()],
#' [moocore::is_nondominated()], [moocore::igd()],
#' [moocore::igd_plus()], [moocore::hypervolume()],
#' [sf::st_as_sf()], [sf::st_concave_hull()],
#' [sf::sf_use_s2()], [aricode::ARI()], [future::future],
#' [future.apply::future.apply], [progressr::progressr]
#' 
#' @references
#' Steuer, R. E., & Choo, E.-U. (1983). An interactive weighted Tchebycheff
#' procedure for multiple objective programming. Mathematical Programming,
#' 26(3), 326–344.
#' 
#' Coello Coello, C. A., & Reyes Sierra, M. (2004). A study of the
#' parallelization of a coevolutionary multi-objective evolutionary algorithm.
#' In MICAI 2004: Advances in Artificial Intelligence (pp. 688–697).
#' Springer Berlin Heidelberg.
#' 
#' Ishibuchi, H., Masuda, H., Tanigaki, Y., & Nojima, Y. (2015). Modified
#' distance calculation in generational distance and inverted generational
#' distance. In Lecture Notes in Computer Science (pp. 110–125).
#' Springer International Publishing.
#'
#' Zitzler, E., & Thiele, L. (1998). Multiobjective optimization using
#' evolutionary algorithms — A comparative case study. In Lecture Notes in
#' Computer Science (pp. 292–301). Springer Berlin Heidelberg.
#' 
#' Miettinen, K. (1998). Nonlinear Multiobjective Optimization (1998th ed.).
#' Springer.
#' 
#' @examples
#' # load demo data
#' data <- threeblobs_s
#' 
#' # optimise within-cluster spatial proximity and temporal coverage
#' set.seed(123)
#' stblob(data = data, k = 2:3, w_space = c(0.5,1))
#' 
#' # optimise within-cluster spatial proximity, temporal coverage and data type diversity
#' set.seed(123)
#' stblob(data = data, k = 2:3, w_space =  c(0.5,1), optim_type_diversity = TRUE)
#' 
#' # optimise with a precomputed spatial distance matrix
#' space_distmat <- threeblobs_s_sdistm
#' set.seed(123)
#' stblob(data = data, k = 2:3, w_space = c(0.5,1), space_distmat = space_distmat)
#' 
#' # specify 'data' columns
#' set.seed(123)
#' stblob(data = data, k = 2:3, w_space = c(0.5,1), optim_type_diversity = TRUE, coords = c("x","y"), age = "age", type = "type")
#'
#' @export

stblob <- function(data,
                   k,
                   w_space,
                   optim_type_diversity = FALSE,
                   iter = 100L,
                   run = 100L,
                   batch = 50L,
                   w_time = NULL, 
                   w_type = NULL,
                   
                   # additional input arguments like column names
                   coords = NULL,
                   age = NULL,
                   type = NULL,
                   space_distmat = NULL,
                   space_distmethod = c("geodesic", "euclidean"),
                   distGeo_a = 6378137,
                   distGeo_f = 1/298.257223563,
                   
                   # MOO convergence setting
                   adapt_search = TRUE,
                   moo_tol = 1e-4,
   
                   # future.apply
                   parallel = TRUE,
                   workers = NULL,

                   # constraints and additional local-search arguments
                   ...
                   ) {
  
  params <- mget(ls(environment(), sorted = FALSE))
  params$data <- params$space_distmat <- NULL
  
  # check data
  check_args(data = data,
             optim_type_diversity = optim_type_diversity,
             iter = iter,
             coords = coords,
             age = age,
             type = type)
  
  # assign column names to NULL
  coords <- coords %||% names(data)[c(2,3)]
  age <- age %||% names(data)[4]
  type <- if (optim_type_diversity) type %||% names(data)[5] else NULL
  
  # check if space_distmat is supplied, compute space_distmat otherwise
  space_distmat <- check_space_distmat(data = data,
                                       coords = coords,
                                       space_distmat = space_distmat,
                                       space_distmethod = space_distmethod,
                                       distGeo_a = distGeo_a,
                                       distGeo_f = distGeo_f)
  
  # in case crs is not specified otherwise
  if (match.arg(space_distmethod, c("geodesic", "euclidean")) == "euclidean") hull_crs <- NA
  
  # check_opt_args
  check_opt_args(...)
  
  # other checks
  stopifnot(
    "'k' must be at least 2 and smaller than the number of samples" = all(k > 1L) & all(k < nrow(data)),
    "'run' must be at least 1" = run > 0L,
    "'parallel' must be logical" = is.logical(parallel),
    "'adapt_search' must be logical" = is.logical(adapt_search),
    "'moo_tol' must be non-negative numeric" = is.numeric(moo_tol) && moo_tol >= 0
  )
  
  if (!adapt_search) {
    pop <- stblob_populate(data = data,
                           k = k,
                           w_space = w_space,
                           optim_type_diversity = optim_type_diversity,
                           iter = iter,
                           run = run,
                           batch = batch,
                           w_time = w_time,
                           w_type = w_type,
                           coords = coords,
                           age = age,
                           type = type,
                           space_distmat = space_distmat,
                           space_distmethod = space_distmethod,
                           distGeo_a = distGeo_a,
                           distGeo_f = distGeo_f,
                           parallel = parallel,
                           workers = workers,
                           ...)
    
    moo <- eval_moo(pop, moo_tol)
    
  } else {
    
    for (i in 1:batch) {
      pop <- stblob_populate(data = data,
                             k = k,
                             w_space = w_space,
                             optim_type_diversity = optim_type_diversity,
                             iter = iter,
                             run = run,
                             batch = 1L,
                             w_time = w_time,
                             w_type = w_type,
                             coords = coords,
                             age = age,
                             type = type,
                             space_distmat = space_distmat,
                             space_distmethod = space_distmethod,
                             distGeo_a = distGeo_a,
                             distGeo_f = distGeo_f,
                             parallel = parallel,
                             workers = workers,
                             .batch = i,
                             ...)
      
      if (i > 1) {
        pop <- combine_pop(pop_prev, pop)
      }
      
      moo <- eval_moo(pop, moo_tol)
      # break if converged
      if (!is.null(moo$moo_convergence) && moo$moo_convergence) break
      
      # record pop and dup from the last cum batch for next batch
      pop_prev <- pop
      k_prev <- k
      
      # update k
      k <- moo$adapt_search_params$k
      # update w_space
      w_space <- update_w_space(w_space = w_space,
                                w_space_q3 = pop$adapt_search_params$w_space_q3,
                                k_prev = k_prev,
                                k = k)
    }
  }
  
  if (!is.null(moo$moo_convergence) && moo$moo_convergence) {
    if (moo$moo_convergence) {
      message(
        paste0("MOO converged (moo_tol: ", moo_tol,")")
      )
    } else {
      warning(
        paste0("MOO not converged (moo_tol: ", moo_tol,")"), call. = FALSE
      )
    }
  }
  
  params_dots <- moo$params[setdiff(names(moo$params), names(params))]
  params <- c(params, params_dots)
  
  moo$params <- params
  return(moo)
}

# eval_moo --------------------------------------------------------------------
eval_moo <- function(pop, moo_tol = 1e-4) {
  # checks
  stopifnot("'pop' must be an S3 'stblob_pop' object" = inherits(pop, "stblob_pop"))
  
  # batch and optim_type_diversity values are essential
  k <- pop$adapt_search_params$k
  batch <- pop$params$batch
  optim_type_diversity <- pop$params$optim_type_diversity
  adapt_search_params <- pop$adapt_search_params
  pop$params$moo_tol <- moo_tol
  
  clust <- pop$clust
  summary <- pop$summary
  trace <- pop$trace
  
  
  # no solution
  if (is.null(clust)) {
    return(new_moo(clust_pf = NULL,
                   clust = pop$clust, # NULL
                   summary = summary, # NULL
                   trace = trace, # NULL
                   moo_convergence = NULL,
                   moo_diagnostics = NULL,
                   filter_summary = pop$filter_summary,
                   data = pop$data,
                   space_distmat = pop$space_distmat,
                   adapt_search_params = adapt_search_params,
                   params = pop$params))
  }
  
  # objective values and space
  z_cols <- c("z_space", "z_time1", "z_time2", "z_type")
  z_cols <- if (optim_type_diversity) z_cols else z_cols[1:3]
  # normalised objective space z
  z <- pop$summary[z_cols]
  z <- as.data.frame(apply(z, 2, minmax, simplify = FALSE))
  # stratify z by k for adapt_search
  z_k <- lapply(k, function(j) z[pop$summary$k == j, ])
  
  # initialise cumulative indices
  b <- integer()
  
  # initialise summary statistics
  igd <- igd_plus <- hv <- numeric()
  # matrix to store stratified hv 
  hv_k <- matrix(nrow = batch, ncol = length(k))
  
  # Pareto front across all solutions
  paretof_idx <- which(moocore::is_nondominated(z))
  paretof_z <- z[paretof_idx, ]
  
  # Pareto front across all solutions stratified by k
  paretof_idx_k <- lapply(z_k, function(x) {
    if (nrow(x) == 0) NA else which(moocore::is_nondominated(x))
  })
  
  paretof_z_k <- Map(function(x,y) {
    if (all(is.na(x))) NA else y[x,]
  }, paretof_idx_k, z_k)
  
  # rank idx
  paretof_idx <- paretof_idx[do.call(order, paretof_z[, z_cols])]
  # output clust_pf to simplify the steps for user to get selected clust
  clust_pf <- clust[paretof_idx, ]
  
  for (i in 1:batch) {
    # subset solutions till batch i
    b <- append(b, i)
    z_i <- z[pop$summary$batch %in% b, ]
    
    # skip empty cumulative batch
    if (nrow(z_i) == 0) next
    
    # find solutions on the Pareto front
    # paretof_idx_i <- which(moocore::is_nondominated(z_i))
    #  objective values of solutions on the Pareto front
    # paretof_z_i <- z_i[paretof_idx_i, ]
    paretof_z_i <- z_i[which(moocore::is_nondominated(z_i)), ]
    
    # compute IGD, IGD+ and Hypervolume
    igd[i] <- moocore::igd(x = paretof_z_i, reference = paretof_z)
    igd_plus[i] <- moocore::igd_plus(x = paretof_z_i, reference = paretof_z)
    hv[i] <- moocore::hypervolume(x = paretof_z_i, reference = rep(1.1, length(z_cols))) # 1.1 such that it would not become 0
    
    # per k_j
    z_ki <- lapply(k, function(j) z[pop$summary$batch %in% b & pop$summary$k == j, ])
    paretof_z_ki <- lapply(z_ki, function(x) {
      if (nrow(x) == 0) NA else x[which(moocore::is_nondominated(x)),]
    })
    
    hv_k[i,] <- vapply(paretof_z_ki, function(x) {
      if (all(is.na(x))) {NA} else {
        moocore::hypervolume(x = x, reference = rep(1.1, length(z_cols)))
      }
    }, numeric(1L))
  }
  
  # label idx of solutions on the Pareto front ordered by z_space, z_time1, z_time2, z_type
  summary$idx_pf <- match(summary$idx, paretof_idx)
  trace$idx_pf <- match(trace$idx, paretof_idx)
  
  # count paretof_n per batch
  paretof_n <- as.data.frame(table(subset(summary, !is.na(idx_pf))[["batch"]]))
  names(paretof_n) <- c("batch", "paretof_n")
  paretof_n <- merge(data.frame(batch = 1:batch), paretof_n, all.x = TRUE)
  paretof_n[is.na(paretof_n)] <- 0L
  
  # MOO indicators 
  moo_diagnostics <- data.frame(batch = 1:batch,
                                igd = igd,
                                igd_plus = igd_plus,
                                hv = hv,
                                paretof_n = paretof_n[["paretof_n"]])
  
  # reorder summary and trace columns
  summary <- summary[c("idx_pf","idx", setdiff(names(summary), c("idx","idx_pf")))]
  summary <- summary[order(summary$idx_pf), ]
  rownames(summary) <- NULL
  trace <- trace[c("idx_pf","idx", setdiff(names(trace), c("idx","idx_pf")))]
  trace <- trace[order(trace$idx_pf), ]
  rownames(trace) <- NULL
  
  # reindex idx in trace
  trace$idx <- match(trace$idx, summary$idx)
  # reorder clust
  clust <- clust[summary$idx, ]
  # reindex idx in summary
  summary$idx <- 1:nrow(summary)
  
  # check convergence
  # percentage change in hv in the last 3 batches
  moo_convergence <- FALSE
  
  # 3 consecutive front = 2 consecutive deltas
  if (i > 2) {
    hv_deltas <- diff(hv) / utils::head(hv, -1)
    hv_deltas <- utils::tail(hv_deltas, 2)
    moo_convergence <- all(hv_deltas < moo_tol)
    
    hv_deltas_k <- diff(hv_k) / utils::head(hv_k, -1)
    hv_deltas_k <- utils::tail(hv_deltas_k, 2)
    # sum == 2 means TRUE * 2 completed, remove from k in the next round by updating adapt_search_param$k
    k_keep <- apply(hv_deltas_k < moo_tol, 2, sum) < 2
    # continue to search for k without valid solutions (thus NA)
    k_keep[is.na(k_keep)] <- TRUE
    adapt_search_params$k <- k[k_keep]
    if (length(adapt_search_params$k) == 0) moo_convergence <- TRUE
  }
  
  moo <- new_moo(clust_pf = clust_pf,
                 clust = clust,
                 summary = summary,
                 trace = trace,
                 moo_convergence = moo_convergence,
                 moo_diagnostics = moo_diagnostics,
                 filter_summary = pop$filter_summary,
                 data = pop$data,
                 space_distmat = pop$space_distmat,
                 adapt_search_params = adapt_search_params,
                 params = pop$params)
  
  return(moo)
}

# helpers ---------------------------------------------------------------------
## combine_pop ----------------------------------------------------------------
combine_pop <- function(pop1, pop2, ...) {
  pop_list <- list(pop1, pop2, ...)
  # in case other arguments are passed here
  # pop_list <- pop_list[vapply(pop_list, function(x) inherits(x, "stblob_pop"), logical(1L))]
  
  # check pop1 and pop2
  stopifnot(
    "Input objects must be of S3 class 'stblob_pop'" =
      all(vapply(pop_list,
                 function(x) inherits(x, "stblob_pop"),
                 logical(1L)))
  )
  
  N <- length(pop_list)
  
  # check for same data and space_distmat
  for (i in 2:N) {
    stopifnot(
      "All pops must be generated with the same 'data'" =
        identical(pop1$data, pop_list[[i]]$data),
      "All pops must be generated with the same 'space_distmat'" =
        identical(pop1$space_distmat, pop_list[[i]]$space_distmat)
    )
  }
  
  # use the last pop$adapt_search_params
  adapt_search_params <- pop_list[[N]]$adapt_search_params
  
  # keep the params k, batch and those have to be identical / in ... (for convenience)
  params_names <- c("k","batch",
                    "optim_type_diversity",
                    "ls_tol","chebyshev_rho",
                    "filter_intersects", "hull_convex_ratio","hull_crs",
                    "filter_clustsize", 
                    "outlier_removal", "outlier_iqrm",
                    "random_init", "sf_use_s2")
  
  params <- lapply(params_names, function(x) {
    lapply(pop_list, function(y) y$params[[x]])
  })
  
  names(params) <- params_names
  
  for (e in params_names) {
    v <- unlist(params[[e]])
    # save k across pops
    if (e == "k") params[[e]] <- sort(unique(v))
    # save the total numbers of batches
    if (e == "batch") params[[e]] <- sum(v)
    # save the param values that have to be identical
    if (!e %in% c("k","batch")) {
      if (length(unique(v)) > 1) stop(sprintf("All pops must be generated with the same'%s' value", e))
      params[[e]] <- v[1]
    }
  }
  
  data <- pop_list[[1]]$data
  space_distmat <- pop_list[[1]]$space_distmat
  
  # compute cumulative batch and solution indices to update summary and trace
  batches <- cum_batches <- counts <- cum_counts <- integer()
  for (i in 1:(N-1)) {
    batches <- append(batches, pop_list[[i]]$params$batch) 
    cum_batches[i+1] <- sum(batches)
    counts <- append(counts, if (!is.data.frame(pop_list[[i]]$summary)) 0L else max(pop_list[[i]]$summary$idx))
    cum_counts[i+1] <- sum(counts)
  }
  
  # the first is always NA 
  cum_batches[is.na(cum_batches)] <- cum_counts[is.na(cum_counts)] <- 0L
  
  # combine clust, summary, trace, filter_summary and adapt_search_params
  pop <- list()
  for (l in c("clust", "summary", "trace", "filter_summary")) {
    # extract list of data frame
    e <- lapply(pop_list, function(x) x[[l]])
    
    # update batch and solution indices
    if (l %in% c("summary", "trace")) {
      e <- Map(function(x,y) { x$batch <- x$batch + y; x }, e, cum_batches)
      e <- Map(function(x,y) { x$idx <- x$idx + y; x }, e, cum_counts)
      # in case of NA
      e <- lapply(e, function(x) { x <- if(!is.data.frame(x)) NULL else x; x })
    }
    
    if (l %in% c("clust", "summary", "trace")) {
      e <- do.call(rbind, e);
      # assign NA if NULL
      e <- e %||% NA
    }
    
    # update filter_summary
    if (l == "filter_summary") {
      # get all k and make sure all status_counts df has same dimension
      K <- unlist(lapply(e, function(x) x$k)); K <- min(K):max(K)
      e <- lapply(e, function(x) { merge(data.frame(k = K), x, all.x = TRUE) })
      # fill 0 for pass, k1, intersects, clustsize, intersects_clustsize, dup
      e <- lapply(e, function(x) { x[2:7][is.na(x[2:7])] <- 0L; x })
      # e: pass, k1, intersects, clustsize, intersects_clustsize, dup
      e <- lapply(e, function(x) x[2:7])
      
      # put them back as a data.frame
      e <- cbind(
        k = K, # k
        Reduce("+", e) # e
      )
    }
    pop[[l]] <- e
  }
  
  clust <- pop$clust
  summary <- pop$summary
  trace <- pop$trace
  filter_summary <- pop$filter_summary
  
  if(!is.null(clust) && nrow(clust) > 1) {
    dup <- duplicated(clust)
    # duplicated indices
    dup_idx <- which(dup)
    
    if (length(dup_idx) > 0) {
      # parameter k to split dup count
      dup_idx_k <- split(dup_idx, summary$k[dup_idx])
      
      # record the freq
      dup_counts <- vapply(dup_idx_k, function(x) length(x), integer(1L))
      dup_counts <- as.data.frame(dup_counts)
      dup_counts <- cbind(as.integer(rownames(dup_counts)), dup_counts)
      names(dup_counts) <- c("k", "dup")
      
      # update status_counts
      dup_counts <- merge(data.frame(k = filter_summary$k), dup_counts, all.x = TRUE)
      dup_counts[is.na(dup_counts)] <- 0L
      # add dup counts
      filter_summary$dup <- filter_summary$dup + dup_counts$dup
      # minus pass counts
      filter_summary$pass <- filter_summary$pass - dup_counts$dup
      
      # remove the runs from the output
      clust <- clust[-dup_idx, , drop = FALSE]
      summary <- subset(summary, !idx %in% dup_idx)
      trace <- subset(trace, !idx %in% dup_idx)
      # reindex the output
      rownames(summary) <- NULL
      summary$idx <- match(summary$idx, unique(summary$idx))
      rownames(trace) <- NULL
      trace$idx <- match(trace$idx, unique(trace$idx))
    }
  }
  
  return(new_pop(clust = clust,
                 summary = summary,
                 trace = trace,
                 filter_summary = filter_summary,
                 data = data,
                 space_distmat = space_distmat,
                 adapt_search_params = adapt_search_params,
                 params = params))
}

## update_w_space -------------------------------------------------------------
update_w_space <- function(w_space, w_space_q3, k_prev, k) {
  
  if (!is.list(w_space)) w_space <- list(w_space)
  lb_prev <- vapply(w_space, min, numeric(1L)) # of k_prev 
  ub <- vapply(w_space, max, numeric(1L))
  
  lb <- w_space_q3$w_space_q3[w_space_q3$k %in% k]
  # take previous value if NA
  if(length(lb_prev) == 1) {
    lb[is.na(lb)] <- lb_prev
    ub <- rep(ub, length(lb))
  } else {
    lb_prev <- lb_prev[k_prev %in% k]
    lb[is.na(lb)] <- lb_prev[is.na(lb)]
    ub <- ub[k_prev %in% k]
  }
  
  # create new w_space parameter output
  w_space <- Map(append, as.list(lb), as.list(ub))
  w_space <- lapply(w_space, unique)
  if (length(unique(w_space)) == 1) w_space <- w_space[[1]]
  return(w_space)
}

## new_stblob -----------------------------------------------------------------
new_moo <- function(clust_pf,
                    clust,
                    summary,
                    trace,
                    moo_convergence,
                    moo_diagnostics,
                    filter_summary,
                    data,
                    space_distmat,
                    adapt_search_params,
                    params) {
  
  stopifnot(is.matrix(clust_pf) || is.numeric(clust_pf) || is.null(clust_pf) ,
            is.matrix(clust) || is.numeric(clust) || is.null(clust),
            is.data.frame(summary) || is.null(summary),
            is.data.frame(trace) || is.null(trace),
            is.logical(moo_convergence) || is.null(moo_diagnostics),
            is.data.frame(moo_diagnostics) || is.null(moo_diagnostics),
            is.data.frame(filter_summary),
            is.data.frame(data),
            is.matrix(space_distmat),
            is.list(adapt_search_params),
            is.list(params))
  
  structure(
    list(clust_pf = clust_pf,
         clust = clust,
         summary = summary,
         trace = trace,
         moo_convergence = moo_convergence,
         moo_diagnostics = moo_diagnostics,
         filter_summary = filter_summary,
         data = data,
         space_distmat = space_distmat,
         adapt_search_params = adapt_search_params,
         params = params),
    
    class = c("stblob_moo")
  )
}
