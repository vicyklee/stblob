#' STblob (SpatioTemporal blob) optimisation algorithm
#' 
#' @description
#' `stblob()` performs multi-objective optimisation (MOO)
#' for a range of cluster numbers (`k`) 
#' and relative spatial weights (`w_space`)
#' to approximate the Pareto front
#' of spatial partitions for regional time series analyses.
#' 
#' @inheritParams stblob_lsearch
#' @inheritParams stblob_populate
#' @param adapt_search a logical value to use adaptive search. `k` and
#' the lower bound of `w_space` are updated every batch to reduce infeasible
#' searches. Default is `TRUE`.
#' @param moo_tol tolerance for multi-objective optimisation. The minimum change
#' in hypervolume for three consecutive batches to be considered reaching
#' convergence.
#' @param dup_tol tolerance for duplicates when `adapt_search = TRUE`. Searches
#' for \eqn{k_i} with a proportion of duplicates more than `dup_tol` will be
#' considered complete.
#' 
#' @details
#' STblob (SpatioTemporal blob) optimisation algorithm performs multi-objective
#' optimisation (MOO) using a local-search algorithm
#' (see [stblob_lsearch()]) in a MOO framework.
#' 
#' By sampling a range of parameter combinations of `w_space`, `w_time` and `w_type`,
#' it populates solutions that locally optimise trade-offs between spatial,
#' temporal and data type cost functions using augmented weighted
#' Chebyshev (Tchebycheff) scalarisation (Steuer & Choo, 1983)
#' across \eqn{k} clusters.
#' Sampling is done via Latin hypercube sampling ([lhs::randomLHS()]) for a
#' single free parameter, otherwise via sampling from a Dirichlet distribution
#' with \eqn{\alpha = (1,1,1)} ([MCMCack::rdirichlet()]).
#' 
#' To approximate the global optimal trade-offs,
#' non-dominated solutions are sorted across cumulative batches
#' to assess the advance of the Pareto front based on the following objectives:
#' 
#' * `space_distance`: average within-cluster spatial distances to the medoid
#' * `time_variance`:  average within-cluster temporal variance
#' * `time_evenness`: average within-cluster temporal variance of first differences
#' * `type_diversity` (`optim_type_diversity = TRUE`): average within-cluster Shannon entropy of data types
#' 
#' Objective values (\eqn{z}) are returned in the `summary`
#' by the terms `z_space`, `z_time1`, `z_time2` and `z_type` respectively.
#' The maximisation objectives are multiplied by -1
#' as a minimisation problem in MOO.
#' 
#' Three commonly used quality indicators,
#' the inverted generational distance (IGD) (Coello Coello & Reyes Sierra, 2004),
#' IGD+ (Ishibuchi et al., 2015)
#' and hypervolume (Zitzler & Thiele, 1998)
#' (see also [moocore::igd()], [moocore::igd_plus()] and [moocore::hypervolume()]),
#' are used as diagnostics of the advance of Pareto front across successive batches.
#' For IGD and IGD+, the non-dominated set of all solutions is used as the
#' reference whereas a vector of \eqn{(1.1,1.1,1.1,...)} for \eqn{m} objectives
#' is used as the reference point for the computation of hypervolume.
#' Note that objective values are normalised using the non-dominated set of all
#' solutions.
#' 
#' `adapt_search = TRUE` does two things to reduce duplicated or infeasible
#' searches.
#' First, it stops searching for \eqn{k_i}
#' that exceeds the proportion of duplicates as specified by `dup_tol`.
#' Second, it updates the lower bound of `w_space` with the upper quartile (Q3)
#' of `w_space` across solutions with intersecting boundaries for each \eqn{k}.
#' 
#' `outlier_removal` is being tested at the moment. The idea is to remove points
#' whose distance to their medoid is considered an outlier based on the
#' interquantile range (IQR) of distances to medoid across \eqn{k} clusters.
#' 
#' @returns
#' an S3 object of class `stblob_moo` with the following components:
#'  * `paretof_idx`: indices of selected solutions on the Pareto front, ordered by
#'  `z_space`, `z_time1`, `z_time2` and `z_type`.
#'  The indices correspond to the rows of `clust`.
#'  * `clust`: a matrix of cluster assignments;
#'  row: solution; column: data point.
#'  * `summary`: a data frame of summary statistics. They include the parameter
#'  value of (`k`), number of clusters of the output (`k_o`), four
#'  objective values (`z_space`, `z_time1`, `z_time2` and `z_type`),
#'  number of iterations (`iter`), ARI with the previous iteration (`ari`),
#'  state of local search convergence (`ls_convergence`)
#'  and number of outliers (`n_outliers`).
#'  `idx`, `run` and `batch` indicate the row in `clust`, run and batch from
#'  which the solution was generated. `paretof` indicates whether the solution
#'  is on the Pareto front.
#'  * `trace`: a data frame of summary statistics for each iteration.
#'  * `moo_convergence`: a logical value to indicate MOO convergence.
#'  * `moo_diagnostics`: a data frame of MOO quality
#'  indicators for every cumulative batch.
#'  * `data`: a data frame of the input data.
#'  * `filter_summary`: a data frame of filter summary.
#'  * `params`: a list of parameter values.
#' 
#' @seealso [compute_distmat()], [stblob_lsearch()], [stblob_populate()],
#' [eval_moo()],
#' [moocore::is_nondominated()], [moocore::igd()],
#' [moocore::igd_plus()], [moocore::igd_hypervolume()],
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
#' @export

stblob <- function(data,
                   k,
                   w_space,
                   optim_type_diversity = FALSE,
                   iter = 10L,
                   run = 10L,
                   batch = 10L,
                   w_time = NULL, 
                   w_type = NULL,
                   
                   # additional input arguments like column names
                   coords = NULL,
                   age = NULL,
                   type = NULL,
                   space_distmat = NULL,
                   space_distmethod = c("geodesic", "euclidean"),
                   
                   # MOO convergence setting
                   adapt_search = TRUE, ###
                   moo_tol = 1e-4, ###
                   dup_tol = 0.9, ### 
   
                   # future.apply
                   parallel = TRUE,
                   workers = NULL,

                   # constraints and additional local-search arguments
                   ...
                   ) {
  
  params <- c(mget(ls(environment(), sorted = FALSE)),
              match.call(expand.dots = FALSE)$...)
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
                                       space_distmethod = space_distmethod)
  
  # in case crs is not specified otherwise
  if (space_distmethod == "euclidean") hull_crs <- NA
  
  # check_opt_args
  check_opt_args(...)
  
  # check populate related args
  stopifnot(
    "'run' must be at least 2" = run > 1L,
    "'batch' must be at least 2" =  batch > 1L,
    "'parallel' must be logical" = is.logical(parallel)
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
                           parallel = parallel,
                           workers = workers,
                           ...)
    
    # eval_moo() ###
    moo <- eval_moo(pop, optim_type_diversity, moo_tol)
    
  } else {
    # progress bar for the loop
    # it will override the one in blob_populate()
    p <- progressr::progressor(along = 1:batch) 
    
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
                             parallel = parallel,
                             workers = workers,
                             ...)
      
      # progress
      p()
      
      if (i > 1) {
        pop <- combine_pop(pop_prev, pop)
        # count the dup from the second pop
        dup <- pop$filter_summary$dup
        dup <- dup - dup_prev
        K <-  pop$filter_summary$k
        
        # check which k has dup > tol, remove k to save search effort
        k_done <- K[dup/run > dup_tol]
        # update k parameter
        if (length(k_done) > 0) k <- k[!k %in% k_done]
        
        w_space <- update_w_space(w_space = w_space,
                                  filter_summary = pop$filter_summary,
                                  k_prev = k_prev,
                                  k = k)
        
        # eval_moo
        moo <- eval_moo(pop, optim_type_diversity, moo_tol)
        # converged then break
        if (!is.na(moo$moo_convergence) && moo$moo_convergence) break
        # k will be emptied, no more searches but not techniquely the same as convergence defined here
        if (length(k) == 0) {
          message(paste0("MOO terminated with too many duplicates (dup_tol:", dup_tol, ")"))
          break
        }
      }
      
      # record pop and dup from the last cum batch for next batch
      pop_prev <- pop
      dup_prev <- pop_prev$filter_summary$dup
      k_prev <- k
    }
  }
  
  if (!is.na(moo$moo_convergence) && moo$moo_convergence) {
    message(
      paste0("MOO converged (moo_tol: ", moo_tol,")")
    )
  } else {
    if (!is.na(moo$moo_convergence) && !moo$moo_convergence)
    warning(
      paste0("MOO not converged (moo_tol: ", moo_tol,")"), call. = FALSE
    )
  }
    
  moo$params <- params
  return(moo)
}

# eval_moo --------------------------------------------------------------------
#' Evaluate multi-objective optimisation (MOO)
#' 
#' @description
#' `eval_moo()` evaluates MOO.
#' 
#' @inheritParams stblob
#' @param pop an S3 `stblob_pop` object.
#' 
#' @inheritSection stblob returns
#' 
#' @export

eval_moo <- function(pop, optim_type_diversity, moo_tol) {
  batch <- sum(unlist(pop$params$batch))
  
  # checks
  stopifnot("'batch' must be at least 2" = batch > 1L,
            "'pop' must be an S3 'stblob_pop' object" = inherits(pop, "stblob_pop"))
  
  # no solution
  if (all(is.na(pop$clust)) && length(pop$clust) == 1) {
    return(new_moo(paretof_idx = NA,
                   clust = pop$clust, # NA
                   summary = pop$summary, # NA
                   trace = pop$trace, # NA
                   moo_convergence = NA,
                   moo_diagnostics = NA,
                   filter_summary = pop$filter_summary,
                   data = pop$data, 
                   params = pop$params))
  }
  
  # 1 solution from 2 batches
  if (nrow(pop$clust) == 1) {
    # insert paretof column as NA
    summary$paretof <- trace$paretof <- NA
    summary <- summary[c("idx","paretof", setdiff(names(summary), c("idx","paretof")))]
    trace <- trace[c("idx","paretof", setdiff(names(trace), c("idx","paretof")))]
    
    message("Only one solution was found. More runs are recommended.")
    return(new_moo(paretof_idx = 1,
                   clust = pop$clust, 
                   summary = summary, 
                   trace = trace, 
                   moo_convergence = NA,
                   moo_diagnostics = NA,
                   filter_summary = pop$filter_summary,
                   data = pop$data, 
                   params = pop$params))
  }
  
  summary <- pop$summary
  trace <- pop$trace
  
    # only for output info 
    # obj <- c("space_distance", "time_variance", "time_evenness", "type_diversity")
    # obj <- if (optim_type_diversity) obj else obj[1:3]
  
  # objective values and space
  z_cols <- c("z_space", "z_time1", "z_time2", "z_type")
  z_cols <- if (optim_type_diversity) z_cols else z_cols[1:3]
  # normalised objective space z
  z <- pop$summary[z_cols]; z <- as.data.frame(apply(z, 2, minmax))
  
  # initialise cumulative indices
  b <- integer()
    # initialise Pareto objective space from the previous cumulative batch
    # paretof_idx_l <- list()
  # initialise summary statistics
  # igd <- igd_plus <- hv <- cum_n <- numeric()
  igd <- igd_plus <- hv <- numeric()
  summary$paretof <- trace$paretof <- FALSE
  
  # Pareto front across all solutions
  # s <- subset(obj_space, select = -batch)
  paretof_idx <- which(moocore::is_nondominated(z))
  paretof_z <- z[paretof_idx, ]
  # rank idx
  paretof_idx <- paretof_idx[do.call(order, paretof_z[, z_cols])]

  for (i in 1:batch) {
    # subset solutions till batch i
    b <- append(b, i)
    z_i <- z[pop$summary$batch %in% b, ]
    
    # skip empty cumulative batch
    if (nrow(z_i) == 0) next
    
    # find solutions on the Pareto front
    paretof_idx_i <- which(moocore::is_nondominated(z_i))
    
    ##################################
    # record the idx for every batch
    # paretof_idx_l[[i]] <- paretof_idx_i
    ##################################
    
    paretof_z_i <- z_i[paretof_idx_i, ]
    
    # compute IGD, IGD+ and Hypervolume
    igd[i] <- moocore::igd(x = paretof_z_i, reference = paretof_z)
    igd_plus[i] <- moocore::igd_plus(x = paretof_z_i, reference = paretof_z)
    hv[i] <- moocore::hypervolume(x = paretof_z_i, reference = rep(1.1, length(z_cols))) # 1.1 such that it would not become 0
  }
  
  # label summary and trace
  summary$paretof[which(summary$idx %in% paretof_idx)] <-
    trace$paretof[which(trace$idx %in% paretof_idx)] <- TRUE

  # count paretof_n per batch
  paretof_n <- as.data.frame(table(subset(summary, paretof == TRUE)[["batch"]]))
  names(paretof_n) <- c("batch", "paretof_n")
  paretof_n <- merge(data.frame(batch = 1:batch), paretof_n, all.x = TRUE)
  paretof_n[is.na(paretof_n)] <- 0

  # MOO indicators 
  moo_diagnostics <- data.frame(batch = 1:batch,
                                igd = igd,
                                igd_plus = igd_plus,
                                hv = hv,
                                paretof_n = paretof_n[["paretof_n"]])

  # reorder summary and trace columns
  summary <- summary[c("idx","paretof", setdiff(names(summary), c("idx","paretof")))]
  trace <- trace[c("idx","paretof", setdiff(names(trace), c("idx","paretof")))]
  
  # check convergence
  # percentage change in hv in the last 3 batches
  moo_convergence <- FALSE

  if (i > 3) {
    hv_deltas <- diff(hv) / head(hv, -1)
    hv_deltas <- tail(hv_deltas, 3)
    moo_convergence <- all(hv_deltas < moo_tol)
  }

  moo <- new_moo(paretof_idx = paretof_idx,
                 clust = pop$clust,
                 summary = summary,
                 trace = trace,
                 moo_convergence = moo_convergence,
                 moo_diagnostics = moo_diagnostics,
                 filter_summary = pop$filter_summary,
                 data = pop$data,
                 # obj = obj,
                 # paretof_idx_l = paretof_idx_l,
                 params = pop$params)
  
  return(moo)
}

# filter_paretof ---------------------------------------------------------------
#' Filter solutions on the Pareto front
#' 
#' @description
#' `filter_paretof()` filters for diverse solutions on the Pareto front.
#'
#' @param stblob an S3 `stblob_moo` object.
#' @param paretof_ari a numeric value of the Adjusted Rand Index (ARI)
#' that sets the similarity threshold between two Pareto optimal solutions.
#' Default is `1`.
#' @param z_dist a numeric value of the minimum normalised Euclidean distance
#' between solutions in the objective space. Default is `0`.
#' 
#' @details
#' Solutions are ranked by their objective values in the order of `z_space`,
#' `z_time1`, `z_time2` and `z_type`.
#' They are selected based on their ARI in the decision space,
#' and/or the normalised Euclidean distance in the objective space
#' to the previous "better" solution.
#' 
#' @returns a vector of selected solution indices. 
#' 
#' @export

filter_paretof <- function(moo, paretof_ari = 0.8, z_dist = 0) {
  stopifnot(
    "'paretof_ari' must be less than 1 and above 0" =
      (is.numeric(paretof_ari) && paretof_ari >= 0 && paretof_ari <= 1),
    "'z_dist' must be between 0 and 1" =
      (is.numeric(z_dist) && z_dist >= 0 && z_dist <= 1),
    "'moo' must be an S3 `stblob_moo` object" =
      inherits(moo, "stblob_moo"),
    "At least 2 solutions are required for filtering" =
      length(moo$paretof_idx) > 1L
  )
  
  summary <- moo$summary
  clust <- moo$clust
  optim_type_diversity <- unlist(moo$params$optim_type_diversity)[1]
  m <- if (optim_type_diversity) 4 else 3
  z_cols <- c("z_space", "z_time1", "z_time2", "z_type"); z_cols <- z_cols[1:m]
  paretof_idx <- moo$paretof_idx
  
  # paretof objective space and clust 
  z <- summary[z_cols]; z <- apply(z, 2, minmax)
  paretof_z <- z[paretof_idx,]
  paretof_clust <- clust[paretof_idx, ]
  paretof_clust[is.na(paretof_clust)] <- 0
  
  # initialise vectors for the loop
  ref <- 1
  l <- length(paretof_idx)
  sel <- logical(l)
  # the highest rank is always selected
  sel[1] <- TRUE

  for (i in 2:l) {
    if (ref == l) break
    # decision space
    ari <- aricode::ARI(paretof_clust[ref, ], paretof_clust[i, ])
    # objective space
    # / sqrt(1*m) to get normalised Euclidean distance
    d <- as.numeric(stats::dist(paretof_z[c(ref,i), ])) / sqrt(1*m)
    
    # if distinct, keep i and compare i to i+1
    if (ari < paretof_ari && d >= z_dist) {
      sel[i] <- TRUE
      ref <- i
    }
  }
  
  paretof_idx <- paretof_idx[sel]
  return(paretof_idx)
}

# print.stblob_moo ------------------------------------------------------------
print.stblob_moo <- function(x, ...) {
  if (!is.null(x$clust)) {
    cat("STblob MOO output of", nrow(x$clust), "solutions\n")
    
    cat("$paretof_idx\n#Indices of Pareto optimal solutions\n")
    print(head(x$paretof_idx, 10))
    if (length(x$paretof_idx) > 10) cat(paste0(" <", length(x$paretof_idx) - 10," more solutions>\n"))
    cat("\n")
    
    cat("$clust\n# cluster assignments:\n")
    print(head(x$clust, c(5, 10)))
    if (ncol(x$clust) > 10) cat(paste0(" <", ncol(x$clust) - 10," more columns (points)>\n"))
    if (nrow(x$clust) > 5) cat(paste0(" <", nrow(x$clust) - 5," more rows (solutions)>\n"))
    cat("\n")
    
    cat("$summary\n# summary:\n")
    print(head(x$summary, 5))
    if (nrow(x$summary) > 5) cat(paste0(" <", nrow(x$summary) - 5," more rows>\n"))
    cat("\n")
    
    cat("$moo_convergence\n# MOO convergence status:\n")
    cat(paste0(x$moo_convergence, " (moo_tol: ", x$params$moo_tol, ")\n"))
    cat("\n")
    
    cat("$moo_diagnostics\n# MOO diagnostics:\n")
    print(head(x$moo_diagnostics, c(5, ncol(x$moo_diagnostics))))
    if (nrow(x$moo_diagnostics) > 5) cat(paste0(" <", nrow(x$moo_diagnostics) - 5," more batches>\n"))
    cat("\n")
    
    cat("$filter_summary\n# number of passed and filtered solutions: \n")
    print(x$filter_summary[c("k",
                             "pass",
                             "k1",
                             "intersects",
                             "clustsize",
                             "intersects_clustsize",
                             "dup")])
    cat("\n")
    cat("# upper quartiles (Q3) of 'w_space' of solutions with intersecting boundaries: \n")
    print(x$filter_summary[c("k",
                             "w_space_q3")])
    cat("\n")
    
    
  } else {
    cat("No feasible solution was found.\n")
    cat("\n")
    
    cat("# number of filtered solutions:\n")
    print(x$status_counts)
  }
  
  invisible(x)
}

# summary.blob ----------------------------------------------------------------
summary.stblob_moo <- function(x, ...) {
  print(x$summary)
  invisible(x)
}

# helpers ---------------------------------------------------------------------
## update_w_space -------------------------------------------------------------
update_w_space <- function(w_space, filter_summary, k_prev, k) {
  ws_q3 <- filter_summary[c("k","w_space_q3")]
  
  if (!is.list(w_space)) w_space <- list(w_space)
  lb_prev <- vapply(w_space, min, numeric(1L)) # of k_prev 
  ub <- vapply(w_space, max, numeric(1L))
  
  lb <- ws_q3$w_space_q3[ws_q3$k %in% k]
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
  return(w_space)
}

## new_stblob -----------------------------------------------------------------
new_moo <- function(paretof_idx,
                    clust,
                    summary,
                    trace,
                    moo_convergence,
                    moo_diagnostics,
                    filter_summary,
                    data,
                    # obj,
                    # paretof_idx_l,
                    params) {
  
  stopifnot(is.numeric(paretof_idx) || (length(paretof_idx) == 1 && all(is.na(paretof_idx))),
            is.matrix(clust) || (length(clust) == 1 && all(is.na(clust))) ,
            is.data.frame(summary) || (length(summary) == 1 && all(is.na(summary))),
            is.data.frame(trace) || (length(trace) == 1 && all(is.na(trace))),
            is.logical(moo_convergence),
            is.data.frame(moo_diagnostics) || (length(moo_diagnostics) == 1 && all(is.na(moo_diagnostics))),
            is.data.frame(filter_summary),
            is.data.frame(data),
            # is.character(obj),
            # is.null(paretof_idx_l) || is.list(paretof_idx_l),
            is.null(params) || is.list(params))
  
  structure(
    list(paretof_idx = paretof_idx,
         clust = clust,
         summary = summary,
         trace = trace,
         moo_convergence = moo_convergence,
         moo_diagnostics = moo_diagnostics,
         filter_summary = filter_summary,
         data = data,
         # obj = obj,
         # paretof_idx_l = paretof_idx_l,
         params = params),
    class = c("stblob_moo", "stblob_pop")
  )
}
