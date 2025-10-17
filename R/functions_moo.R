#------------------------------------------------------------------------------#
# Multi-objective optimisation ####
#------------------------------------------------------------------------------#
#' Populate solutions by weighted sum scalarisation for multiple batches
#' 
#' @description
#' This function populates solutions by weighted sum scalarisation of the bi-objective function in [blob_search()] for a given k. 
#' 
#' @inheritParams start_blobs
#' @inheritParams find_blobs
#' @inheritParams blob_search
#' @inheritParams blob_populate
#' @param k an integer value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the number of clusters.
#' @param r a numeric value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the relative spatial weight.
#' They must be \eqn{[0,1]}. Default is c(0.5,1).
#' @param batch an integer of the number of batches. Default is 5.
#'
#' @inherit blob_populate details
#' @returns a list of `pop` objects; see also [blob_populate()].
#' @inherit blob_populate seealso
#' @export

blob_populate_batch <- function(data,
                                k,
                                r = c(0.5,1),
                                iter = 10,
                                run = 100,
                                batch = 5,
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
    #-------------------------------------------------------------------#
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
  #-------------------------------------------------------------------#
  # sample r
  if (length(r) == 2) {
    # LHS sampling for more evenly distributed parameters
    lhs_samples <- lhs::randomLHS(run,1)
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
    k[1] <- 2L # lower bound
    k[2] <- as.integer(floor(sqrt(N))) # upper bound, an heuristic to maximise information e.g. 100 points: 10 blobs, 10 points each 
  }
  #-------------------------------------------------------------------#
  r_vec <- r_samples
  batch_vec <- 1:batch
  #-------------------------------------------------------------------#
  if(length(k) == 1) {
    grid <- expand.grid(r_vec = r_vec, batch_vec = batch_vec)
    #-------------------------------------------------------------------#
    blob_list <- future.apply::future_Map(function(r, batch) {
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
        batch = batch
      )
    }, grid$r_vec, grid$batch_vec, future.seed = TRUE)
    #-------------------------------------------------------------------#
    # group runs of same batch
    batch_list <- vector("list", batch)
    for (i in 1:length(blob_list)) {
      m <- blob_list[[i]][["batch"]]
      # append to pop_list
      batch_list[[m]] <- append(batch_list[[m]],
                                list(blob_list[[i]][["blob"]])
      )
    }
    #-------------------------------------------------------------------#
    for (i in 1:batch) {
      blob_list <- batch_list[[i]]
      pop <- convert_to_pop(blob_list)
      pop$summary$batch <- i
      pop$trace$batch <- i
      #-------------------------------------------------------------------#
      pop$summary <- pop$summary[, c("idx",
                                     "batch",
                                     "k", "r", "run",
                                     "space_wcss",
                                     "time_range_mean", "time_range_sd",
                                     "time_evenness_mean", "time_evenness_sd",
                                     "size_mean", "size_sd", "size_diff",
                                     "intersects", "n_removed",
                                     "iter", "ari", "dup")]
      pop$trace <- pop$trace[, c("idx",
                                 "batch",
                                 "k", "r", "run",
                                 "space_wcss",
                                 "time_range_mean", "time_range_sd",
                                 "time_evenness_mean", "time_evenness_sd",
                                 "size_mean", "size_sd", "size_diff",
                                 "intersects",
                                 "iter", "ari")]
      #-------------------------------------------------------------------#
      pop$space_kmat_optim_out = space_kmat_optim_out
      #-------------------------------------------------------------------#
      class(pop) <- "pop"
      batch_list[[i]] <- pop
    }
  } else {
    if (length(k) == 2) {
      k_vec <- k[1]:k[2]
      grid <- expand.grid(k_vec = k_vec, r_vec = r_vec, batch_vec = batch_vec)
      #-------------------------------------------------------------------#
      blob_list <- future.apply::future_Map(function(k, r, batch) {
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
          k = k,
          batch = batch
        )
      }, grid$k_vec, grid$r_vec, grid$batch_vec, future.seed = TRUE)
      #-------------------------------------------------------------------#
      # group runs of same batch
      batch_list <- vector("list", batch)
      for (i in 1:length(blob_list)) {
        m <- blob_list[[i]][["batch"]]
        # append to pop_list
        batch_list[[m]] <- append(batch_list[[m]],
                                  list(
                                    list(blob = blob_list[[i]][["blob"]],
                                         k = blob_list[[i]][["k"]])
                                  )
        )
      }
      #-------------------------------------------------------------------#
      # group runs of same k
      for (i in 1:batch) {
        pop_list <- vector("list", k[2])
        for (j in 1:length(batch_list[[i]])) {
          n <- batch_list[[i]][[j]][["k"]]
          # append to pop_list
          pop_list[[n]] <- append(pop_list[[n]], list(batch_list[[i]][[j]][["blob"]]))
        }
        pop_list <- pop_list[-1] # as the loop starts at 2
        pop_list <- lapply(pop_list, convert_to_pop)
        #-------------------------------------------------------------------#
        # extract each element and add a column to indicate the initial k
        pop <- do.call(rbind, pop_list)
        summary_list <- lapply(seq_along(pop[ , "summary"]), function (i) cbind(pop[ , "summary"][[i]], k_o = i + 1))
        trace_list <- lapply(seq_along(pop[ , "trace"]), function (i) cbind(pop[ , "trace"][[i]], k_o = i + 1)) 
        clust_list <- pop[ ,"clust"]
        n_filtered_list <- pop[ , "n_filtered"]
        #-------------------------------------------------------------------#
        # redundant to store the whole blob data frame at this step so only clust is kept as output from blobs
        clust <- do.call(rbind, clust_list)
        summary <- do.call(rbind, summary_list)
        trace <- do.call(rbind, trace_list)
        n_filtered <- do.call(rbind, n_filtered_list)
        n_filtered <- cbind(k_o = 2:(nrow(n_filtered) + 1), n_filtered)
        #-------------------------------------------------------------------#
        # reindex the solutions
        summary$idx <- NULL
        summary$idx <- 1:nrow(summary)
        trace$idx <- NULL
        trace <- merge(trace, summary[, c("idx","k_o","run")], by = c("k_o","run"), all.x = TRUE)
        summary$batch <- i
        trace$batch <- i
        summary <- summary[, c("idx",
                               "batch",
                               "k_o", "k", "r", "run",
                               "space_wcss",
                               "time_range_mean", "time_range_sd",
                               "time_evenness_mean", "time_evenness_sd",
                               "size_mean", "size_sd", "size_diff",
                               "intersects", "n_removed",
                               "iter", "ari", "dup")]
        trace <- trace[, c("idx",
                           "batch",
                           "k_o", "k", "r", "run",
                           "space_wcss",
                           "time_range_mean", "time_range_sd",
                           "time_evenness_mean", "time_evenness_sd",
                           "size_mean", "size_sd", "size_diff",
                           "intersects",
                           "iter", "ari")]
        
        rownames(summary) <- NULL
        rownames(trace) <- NULL
        #-------------------------------------------------------------------#
        pop <- list(clust = clust,
                    summary = summary,
                    trace = trace,
                    n_filtered = n_filtered,
                    space_kmat_optim_out = space_kmat_optim_out)
        class(pop) <- "pop"
        batch_list[[i]] <- pop
      }
    }
  }
  return(batch_list)
}

#------------------------------------------------------------------------------#
#' Combine pop objects
#' 
#' @description
#' This function combine two pop objects.
#'
#' @param pop_a,pop_b a `pop`, `pop_pareto` or `pop_moo` object;
#' see also [blob_populate()], [find_pareto()] and [eval_moo()].
#'
#' @returns a `pop` object.
#' @export

combine_pop <- function(pop_a, pop_b) {
  # check if any pop is NULL
  if (is.null(pop_a) & is.null(pop_b)) return(NULL)
  # if pop_a went through find_pareto()
  if (!is.null(pop_a)) {
    if(!is.null(pop_a$pareto_idx)) {
      pop_a$summary$pareto <- NULL
      pop_a$summary$pareto_similar <- NULL
      pop_a$trace$pareto <- NULL
      pop_a$trace$pareto_similar <- NULL
      pop_a$pareto_idx <- NULL
      pop_a$obj <- NULL
      message("Pareto optimality has to be re-evaluated.")
    }
  }
  #-------------------------------------------------------------------#
  # if pop_b went through find_pareto()
  if (!is.null(pop_b)) {
    if(!is.null(pop_b$pareto_idx)) {
      pop_b$summary$pareto <- NULL
      pop_b$summary$pareto_similar <- NULL
      pop_b$trace$pareto <- NULL
      pop_b$trace$pareto_similar <- NULL
      pop_b$pareto_idx <- NULL
      pop_b$obj <- NULL
      message("Pareto optimality has to be re-evaluated.")
    }
  }
  #-------------------------------------------------------------------#
  # if one of the pop is NULL, return the !NULL one
  if (is.null(pop_a)){
    if (is.null(pop_b$summary$batch)) {
      pop_b$summary$batch <- 2
      pop_b$trace$batch <- 2
      return(pop_b)
    }
  } else {
    if (is.null(pop_b)) {
      if (is.null(pop_a$summary$batch)) {
        pop_a$summary$batch <- 1
        pop_a$trace$batch <- 1
      }
      return(pop_a)
    }
  }
  #-------------------------------------------------------------------#
  # if both pop is !NULL
  # index the batch
  # check if batch is a column in the output
  if (is.null(pop_a$summary$batch)) {
    pop_a$summary$batch <- 1
    pop_a$trace$batch <- 1
  }
  #-------------------------------------------------------------------#
  if (is.null(pop_b$summary$batch)) {
    pop_b$summary$batch <- max(pop_a$summary$batch) + 1
    pop_b$trace$batch <- max(pop_a$trace$batch) + 1
  }
  #-------------------------------------------------------------------#
  # combine pop
  pop <- Map(rbind, pop_a, pop_b)
  #-------------------------------------------------------------------#
  # reindex summary and trace
  pop$summary$idx <- 1:nrow(pop$summary)
  pop$trace$idx <- NULL # remove column for the new
  if("k_o" %in% names(pop$trace)) {
    pop$trace <- merge(pop$trace, unique(pop$summary[ , c("idx","batch","k_o","run")]), by = c("batch","k_o","run"), all.x = TRUE)
    pop$trace <- pop$trace[, c(length(pop$trace), 1:length(pop$trace)-1)]
    #-------------------------------------------------------------------#
    # sum the filtered counts by k_o
    pop$n_filtered <- cbind(
      k_o = sort(unique(pop$n_filtered[, 1])),
      rowsum(x = pop$n_filtered[ , -1], group = pop$n_filtered[ , 1])
      )
    rownames(pop$n_filtered) <- NULL
  } else {
    pop$trace <- merge(pop$trace, unique(pop$summary[ , c("idx","batch","run")]), by = c("batch","run"), all.x = TRUE)
    pop$trace <- pop$trace[, c(length(pop$trace), 1:length(pop$trace)-1)]
    #-------------------------------------------------------------------#
    # sum the filtered counts
    pop$n_filtered <- rowsum(x = pop$n_filtered, group = rep(1, nrow(pop$n_filtered))) # colSums can does the same but just want to keep it as a data frame without additional functions.
  }
  #-------------------------------------------------------------------#
  # find and filter exact duplicates
  if (!is.null(pop$clust)) {
    if (nrow(pop$clust) > 1) {
      dup <- find_dup(pop$clust, ari = 1)
      #-------------------------------------------------------------------#
      if (length(dup$idx) > 0) {
        # record the duplicate freq
        pop$summary$dup[as.numeric(names(dup$freq))] <- pop$summary$dup[as.numeric(names(dup$freq))] + as.vector(dup$freq)
        # filter the duplicates
        pop$clust <- pop$clust[-dup$idx, , drop = F]
        pop$summary <- pop$summary[-dup$idx, ]
        pop$trace <- subset(pop$trace, !idx %in% dup$idx)
        #-------------------------------------------------------------------#
        # append the dup count
        if ("k_o" %in% names(pop$n_filtered)) {
          # extract from summary
          dup_count <- pop$summary[ , c("k_o","dup")]
          dup_count <- cbind(
            k_o = sort(unique(dup_count$k_o)),
            rowsum(x = dup_count[ , "dup", drop = FALSE], group = dup_count$k_o)
          )
          pop$n_filtered$dup <- dup_count$dup
        } else {
          pop$n_filtered$dup <- pop$n_filtered$dup + length(dup$idx)
        }
        #-------------------------------------------------------------------#
        # reindex the solution
        pop$summary$idx <- 1:nrow(pop$summary)
        pop$trace$idx <- NULL # remove column for the new
        if("k_o" %in% names(pop$trace)) {
          pop$trace <- merge(pop$trace, unique(pop$summary[ , c("idx","batch","k_o","run")]), by = c("batch", "k_o","run"), all.x = TRUE)
          pop$trace <- pop$trace[, c(length(pop$trace), 1:length(pop$trace)-1)]
        } else {
          pop$trace <- merge(pop$trace, unique(pop$summary[ , c("idx","batch","run")]), by = c("batch","run"), all.x = TRUE)
          pop$trace <- pop$trace[, c(length(pop$trace), 1:length(pop$trace)-1)]
        }
      }
    }
  }
  #-------------------------------------------------------------------#
  rownames(pop$summary) <- NULL
  rownames(pop$trace) <- NULL
  #-------------------------------------------------------------------#
  class(pop) <- "pop"
  return(pop)
}

#------------------------------------------------------------------------------#
#' Parse objective parameter
#' 
#' @description
#' This function parses `obj` string vector.
#' 
#' @param obj a string vector of objectives. "-" prefix to indicate objectives to be maximised.
#' @details
#' A leading "-" (minus) is used to indicate an objective to be maximised (e.g. "-obj").
#' 
#' @returns a list of the following objects.
#' \itemize{
#'   \item \code{obj}: a string vector of parsed objectives.
#'   \item \code{maximise_obj_idx}: an numeric vector of indices of the objective to be maximised.
#' }

parse_obj <- function(obj) {
  maximise_obj_idx <- grep(pattern = "^-", obj)
  obj <- sub("^-", "", obj)
  out <- list(obj = obj, maximise_obj_idx = maximise_obj_idx)
  return(out)
}

#------------------------------------------------------------------------------#
#' Find Pareto optimal solutions
#' 
#' @description
#' This function finds a set of Pareto optimal solutions from a given set of solutions using [moocore::is_nondominated()] and updates `pop`.
#' 
#' @param pop a `pop` object; see also [blob_populate()].
#' @param obj a string vector of objectives. "-" prefix to indicate objectives to be maximised.
#' 
#' @returns a `pop_pareto` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of cluster assignments. Each row is a Pareto optimal solution.
#'   \item \code{summary}: a data frame of summary statistics of all feasible solutions.
#'   \item \code{trace}: a data frame of summary statistics for tracing of all feasible solutions.
#'   \item \code{n_filtered}: a data frame of numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} in [distmat_to_kmat()] when `space_kmat` is not supplied.
#'   \item \code{pareto_idx}: a numeric vector of indices of Pareto optimal cluster assignment.
#'   \item \code{obj}: a string vector of objectives.
#' }
#' 
#' @seealso [moocore::is_nondominated()]
#' @export

find_pareto <- function(pop, obj = NULL) {
  if (is.null(obj)) obj <- pop$obj
  if (is.null(obj)) stop("Missing obj!")
  #-------------------------------------------------------------------#
  # parse obj
  obj_input <- obj
  parse_obj_out <- parse_obj(obj)
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  obj <- parse_obj_out$obj
  # subset objspace from summary
  objspace <- pop$summary[ , obj]
  #-------------------------------------------------------------------#
  # multiply obj to be maximised by -1
  if (length(maximise_obj_idx) > 0) {
    for(j in maximise_obj_idx) {
      objspace[[ obj[j] ]] <- -objspace[[ obj[j] ]]
    }
  }
  #-------------------------------------------------------------------#
  # extract the idx of Pareto optimal solutions
  pareto_idx <- which(moocore::is_nondominated(objspace) == TRUE)
  #-------------------------------------------------------------------#
  # update pop
  # create a column in summary to indicate Pareto front solutions
  pop$summary$pareto <- NULL
  pop$summary$pareto <- 0
  pop$summary$pareto[pareto_idx] <- 1
  pop$trace$pareto <- NULL
  pop$trace <- merge(pop$trace, pop$summary[,c("idx","pareto")], by = "idx", all.x = TRUE)
  #-------------------------------------------------------------------#
  # output orignal idx for id in the summary table
  pop$pareto_idx <- pareto_idx
  pop$obj <- obj_input
  #-------------------------------------------------------------------#
  pop_pareto <- pop
  class(pop_pareto) <- "pop_pareto"
  return(pop_pareto)
}

#------------------------------------------------------------------------------#
#' Evaluate multi-objective optimisation performance
#' 
#' @description
#' This function evaluates multi-objective optimisation (MOO) performance.
#' 
#' @param batch_list a list of `pop` objects from different batches; see also [blob_populate_batch()].
#' @inheritParams find_pareto
#'
#' @details
#' The quality indicators include IGD, IGD+ and hypervolume commonly used in MOO,
#' computed using [moocore::igd()], [moocore::igd_plus()], and [moocore::hypervolume()].
#' 
#' The reference for IGD and IGD+ use the Pareto optimal solutions of all batches combined and
#' the reference point for hypervolume uses \eqn{(x_1, x_2, ..., x_n)},
#' where n is the number of objectives and \eqn{x_i = 1} for all \eqn{i}.
#' 
#' The lower and upper bounds for normalisation are obtained from the Pareto front considering all batches combined.
#' 
#' See [moocore::moocore] for more information.
#'
#' @seealso [combine_pop()], [find_pareto()], [moocore::moocore]
#' @returns a `pop_pareto` object includes a list of the following objects.
#' \itemize{
#'   \item \code{clust}: a numeric matrix of cluster assignments. Each row is a Pareto optimal solution.
#'   \item \code{summary}: a data frame of summary statistics of all feasible solutions.
#'   \item \code{trace}: a data frame of summary statistics for tracing of all feasible solutions.
#'   \item \code{n_filtered}: a data frame of numbers of filtered solutions.
#'   \item \code{space_kmat_optim_out}: an output of [stats::optim()] from the optimisation of \eqn{\beta} in [distmat_to_kmat()] when `space_kmat` is not supplied.
#'   \item \code{pareto_idx}: a numeric vector of indices of Pareto optimal cluster assignment.
#'   \item \code{obj}: a string vector of objectives.
#'   \item \code{moo_quality}: a data frame of MOO quality indicators.
#' }
#'
#' @export

eval_moo <- function(batch_list, obj) {
  # parse obj
  obj_input <- obj
  parse_obj_out <- parse_obj(obj)
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  obj <- parse_obj_out$obj
  n_obj <- length(obj)
  #-------------------------------------------------------------------#
  # get a list of Pareto optimal objspaces for the accumulated batches
  N <- length(batch_list)
  pareto_objspace_list <- list()
  for (i in 1:N) {
    # problem!
    if(i > 1) { pop <- suppressMessages(combine_pop(pop, batch_list[[i]])) } else {
      pop <- batch_list[[i]]
    }
    pop_pareto <- find_pareto(pop, obj = obj_input)
    pareto_summary <- subset(pop_pareto$summary, pareto == 1)
    pareto_objspace <- pareto_summary[, obj]
    # inverse the max obj
    if (length(maximise_obj_idx) > 0) {
      for(j in maximise_obj_idx) {
        pareto_objspace[[ obj[j] ]] <- -pareto_objspace[[ obj[j] ]]
      }
    }
    pareto_objspace_list[[i]] <- pareto_objspace
  }
  #-------------------------------------------------------------------#
  # obtain upper and lower bounds from the total Pareto front
  ub <- vapply(1:n_obj, function(i) max(pareto_objspace[[obj[i]]]), numeric(1))
  lb <- vapply(1:n_obj, function(i) min(pareto_objspace[[obj[i]]]), numeric(1))
  #-------------------------------------------------------------------#
  # normalise the objectives
  pareto_objspace_list <- lapply(pareto_objspace_list, function(x) moocore::normalise(x, to_range = c(0,1), lower = lb, upper = ub))
  #-------------------------------------------------------------------#
  # compute IGD, IGD plus and HV
  igd <- vapply(1:(N - 1),
                function(i)
                  moocore::igd(x = pareto_objspace_list[[i]], reference = pareto_objspace_list[[i+1]]),
                numeric(1))
  
  igd_plus <- vapply(1:(N - 1),
                     function(i)
                       moocore::igd_plus(x = pareto_objspace_list[[i]], reference = pareto_objspace_list[[i+1]]),
                     numeric(1))
  
  hv <- vapply(1:N,
               function(i)
                 moocore::hypervolume(x = pareto_objspace_list[[i]], reference = rep(1, n_obj)),
               numeric(1))
  
  igd <- append(NA, igd)
  igd_plus <- append(NA, igd_plus)
  #-------------------------------------------------------------------#
  moo_quality <- data.frame(igd = igd, igd_plus = igd_plus, hv = hv)
  pop_pareto$moo_quality <- moo_quality
  pop_moo <- pop_pareto
  class(pop_moo) <- "pop_moo"
  return(pop_moo)
}

#------------------------------------------------------------------------------#
#' Find similar Pareto optimal solutions 
#' 
#' @description
#' This function finds similar Pareto optimal solutions and updates a `pop_pareto` object.
#'
#' @inheritParams find_pareto
#' @inheritParams find_dup
#' @param pop a `pop_pareto` or `pop_moo` object; see also [find_pareto()] and [eval_moo()].
#' 
#' @details
#' The designation of a similar (redundant) solution is dependent on the order of objectives.
#' Objectives should be ranked in descending order based on their priority.
#' 
#' @returns a `pop_pareto` or `pop_moo` object depending on the input.
#' @seealso [find_dup()]
#'
#' @export

find_pareto_similar <- function(pop, ari) {
  obj <- pop$obj
  # when this is to be reapplied
  pop$summary$pareto_similar <- NULL
  
  #-------------------------------------------------------------------#
  # parse parameter obj to get the column indices and parsed obj names
  parse_obj_out <- parse_obj(obj)
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  # use the parsed obj names in the following steps
  obj <- parse_obj_out$obj
  # select the columns of the objective space
  pareto_objspace <- subset(pop$summary, pareto == 1)[ , obj]
  pareto_idx <- subset(pop$summary, pareto == 1)$idx
  pareto_clust <- pop$clust[pareto_idx, ]
  
  # multiply maximising obj by -1
  if (length(maximise_obj_idx) > 0) {
    for(j in maximise_obj_idx) {
      pareto_objspace[[ obj[j] ]] <- -pareto_objspace[[ obj[j] ]]
    }
  }
  #-------------------------------------------------------------------#
  # find_dup() is order dependent (i.e. it keeps the first row between the dup), therefore when it comes to similarity we want to prioritise based on the objectives
  # order the solutions by the objectives
  pareto_objspace$idx <- 1:nrow(pareto_objspace)
  ordered_idx <- eval(
    parse(
      text = paste0("pareto_objspace[order(", paste0(paste0("pareto_objspace$", obj), collapse = ","), "), ]$idx")
    ) 
  )
  ordered_clust <- pareto_clust[ordered_idx, ]
  #-------------------------------------------------------------------#
  similar <- find_dup(clust = ordered_clust, ari = ari)
  # map the idx in the find_dup() to the original idx
  # idx
  similar$idx <- pareto_idx[ordered_idx[similar$idx]]
  # freq
  names(similar$freq) <- pareto_idx[ordered_idx[as.numeric(names(similar$freq))]]
  # pairs_dup
  similar$pairs_dup <- apply(similar$pairs_dup, c(1,2), function(x) pareto_idx[ordered_idx[x]])
  #-------------------------------------------------------------------#
  # flag similar if any
  pop$summary$pareto_similar <- NULL
  pop$summary$pareto_similar <- 0
  pop$trace$pareto_similar <- NULL
  pop$trace$pareto_similar <- 0
  if (length(similar$idx) > 0) {
    pop$summary$pareto_similar[similar$idx] <- 1
    pop$trace$pareto_similar[pop$trace$idx %in% similar$idx] <- 1
    pop$pareto_idx <- subset(pop$summary, pareto == 1 & pareto_similar == 0)$idx
  }
  return(pop)
}

#------------------------------------------------------------------------------#
#' Multi-objective optimisation (MOO) for spatiotemporal clustering 
#'
#' @description
#' This function populates solutions by weighted sum scalarisation of the bi-objective function in [blob_search()],
#' optimises multiple objectives under constraints and returns a set of Pareto optimal solutions together with MOO quality indicators.
#' 
#' @inheritParams start_blobs
#' @inheritParams find_blobs
#' @inheritParams blob_search
#' @inheritParams blob_populate
#' @inheritParams blob_populate_batch
#' @inheritParams find_pareto
#' @param k an integer value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the number of clusters.
#' @param r a numeric value or vector of length 2. If a vector is supplied, they specify the lower and upper bounds of the relative spatial weight.
#' They must be \eqn{[0,1]}. Default is c(0.5,1).
#' @param pareto_similar_ari a numeric value of Adjusted Rand Index (ARI) that
#' sets similarity threshold between two Pareto optimal solutions.
#' It must be \eqn{[0,1]}. Default is NULL.
#'
#' @details
#' This function is a wrapper of [blob_populate_batch()], [find_pareto()], [eval_moo()] and [find_pareto_similar()].
#' 
#' @inherit eval_moo return
#'
#' @seealso [distmat_to_kmat()], [blob_populate_batch()],
#' [find_pareto()], [eval_moo()], [find_pareto_similar()],
#' [sf::st_as_sf()], [lhs::randomLHS()], [mclust::adjustedRandIndex()],
#' [future::future], [future.apply::future.apply],
#' [moocore::is_nondominated()]
#'
#' @export

blob_moo <- function (data,
                      k,
                      r = c(0.5,1),
                      iter = 10,
                      run = 100,
                      batch = 5,
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
                      beta_par = NULL,
                      pareto_similar_ari = NULL,
                      obj = c("space_wcss","-time_range_mean","-time_evenness_mean","n_removed")) {
                      
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
  #-------------------------------------------------------------------#
  # run batches of blob_populate() using blob_populate_batch()
  batch_list <- blob_populate_batch(data = data,
                                    k = k,
                                    r = r,
                                    iter = iter,
                                    run = run,
                                    batch = batch,
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
  #-------------------------------------------------------------------#
  # evaluate MOO
  pop_moo <- eval_moo(batch_list = batch_list, obj = obj)
  #-------------------------------------------------------------------#
  # if there are feasible solutions, there must be optimal solutions in the set.
  if (is.null(pop_moo$pareto_idx)) {
    message("No assignment was found.")
    return(NULL)
  }
  #-------------------------------------------------------------------#
  # flag similar solutions on the Pareto front
  if (!is.null(pareto_similar_ari)) {
    pop_moo <- find_pareto_similar(pop = pop_moo, ari = pareto_similar_ari)
  }
  #-------------------------------------------------------------------#
  return(pop_moo)
}
