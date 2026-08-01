#' Combine pop objects
#' 
#' @description
#' `combine_pop()` combines `pop` objects.
#' 
#' @param ... at least two arguments of an S3 `pop` object.
#' 
#' @inheritSection blob_populate returns
#' 
#' @export

combine_pop <- function(...) {
  pop_list <- list(...)
  N <- length(pop_list)

  # checks
  stopifnot("Input must be S3 `pop` objects" = all(sapply(pop_list, function(x) inherits(x, "pop"))),
            "There must be at least two arguments" = length(pop_list) > 1)
  
  # combine the params
  params_names <- c("k", "r", "iter", "run", "batch", "converge_ari",
                    "filter_intersects", "hull_convex_ratio", "filter_clustsize",
                    "outlier_removal", "outlier_iqrm")
  params <- lapply(params_names, function(x) unlist(lapply(pop_list, function(y) y$params[[x]])))
  names(params) <- params_names
  # they should have the same weights
  params["weights"] <- list(pop_list[[1]]$params$weights) 
  
  # data to output
  data <- pop_list[[1]]$data
  
  # compute cumulative batch and solution indices to update summary and trace
  batches <- cum_batches <- counts <- cum_counts <- numeric()
  for (i in 1:(N-1)) {
    batches <- append(batches, params$batch[i])
    cum_batches[i+1] <- sum(batches)
    counts <- append(counts, max(pop_list[[i]]$summary$idx))
    cum_counts[i+1] <- sum(counts)
  }
  cum_batches[is.na(cum_batches)] <- cum_counts[is.na(cum_counts)] <- 0

  # combine clust, summary and trace
  pop <- list()
  # pop_idx <- lapply(1:length(pop_list), function(i) data.frame(pop = i))
  for (l in c("clust", "summary", "trace", "filtered_counts")) {
    # extract the df
    e <- lapply(pop_list, function(x) x[[l]])
    
    # update batch and solution indices
    if (l %in% c("summary", "trace")) {
      e <- Map(function(x,y) {x$batch <- x$batch + y; x}, e, cum_batches)
      e <- Map(function(x,y) {x$idx <- x$idx + y; x}, e, cum_counts)
    }
    
    if (l %in% c("clust", "summary", "trace")) e <- do.call(rbind, e)
    
    # update filtered counts
    if (l == "filtered_counts") {
      status <- pop_list[[1]]$filtered_counts$status
      n <- colSums(do.call(rbind, lapply(e, function(x) x$n)))
      e <- data.frame(status = status, n = n)
    }
    pop[[l]] <- e
  }
  
  clust <- pop$clust
  summary <- pop$summary
  trace <- pop$trace
  filtered_counts <- pop$filtered_counts
  
  # remove stblob columns from stblob objects
  summary$sel <- summary$pareto <- trace$sel <- trace$pareto <- NULL
  
  if(nrow(clust) > 1) {
    dup <- find_dup(clust)
    # duplicated indices
    dup_idx <- dup$idx
    
    if (length(dup_idx) > 0) {
      # record the freq
      a_idx <- dup$pairs_dup[1, ] # being duplicated
      b_idx <- dup$pairs_dup[2, ] # duplicate
      b_n <- summary$dup[b_idx] # previous dup freq from its own pop
      # append dup to include dup freq from previous its own pop
      tot_a_n <- dup$freq$n + sapply(split(b_n, a_idx),sum)
      
      summary$dup[a_idx] <- summary$dup[a_idx] + tot_a_n
      filtered_counts$n[5] <- filtered_counts$n[5] + sum(dup$freq$n)
      # remove the runs from the output
      clust <- clust[-dup_idx, ]
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
                 filtered_counts = filtered_counts,
                 data = data,
                 params = params))
}
