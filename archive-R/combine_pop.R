#' Combine populations of searches
#' 
#' @description
#' `combine_pop()` combines `stblob_pop` objects.
#' 
#' @param pop1,pop2 an S3 `stblob_pop` object.
#' @param ... further S3 `stblob_pop` objects to combine.
#' 
#' @inherit stblob_populate return
#'
#' @export

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
  
  # keep the params k, batch and those have to be identical
  params_names <- c("k","batch",
                    "optim_type_diversity","filter_intersects",
                    "hull_convex_ratio","hull_crs","filter_clustsize",
                    "outlier_removal","outlier_iqrm","sf_use_s2")
  
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
  batches <- cum_batches <- counts <- cum_counts <- numeric()
  for (i in 1:(N-1)) {
    batches <- append(batches, pop$params$batch[[i]]) 
    cum_batches[i+1] <- sum(batches)
    counts <- append(counts, if (!is.data.frame(pop_list[[i]]$summary)) 0 else max(pop_list[[i]]$summary$idx))
    cum_counts[i+1] <- sum(counts)
  }
  
  # the first is always NA 
  cum_batches[is.na(cum_batches)] <- cum_counts[is.na(cum_counts)] <- 0

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
      e <- lapply(e, function(x) { merge(data.frame(k = K), x, all.x = T) })
      # fill 0 for pass, k1, intersects, clustsize, intersects_clustsize, dup
      e <- lapply(e, function(x) { x[2:7][is.na(x[2:7])] <- 0; x })
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
      dup_counts <- vapply(dup_idx_k, function(x) length(x), numeric(1L))
      dup_counts <- as.data.frame(dup_counts)
      dup_counts <- cbind(as.numeric(rownames(dup_counts)), dup_counts)
      names(dup_counts) <- c("k", "dup")
      
      # update status_counts
      dup_counts <- merge(data.frame(k = filter_summary$k), dup_counts, all.x = TRUE)
      dup_counts[is.na(dup_counts)] <- 0
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
