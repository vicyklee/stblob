# compute_spacecost -----------------------------------------------------------
compute_space_costs <- function(space_distmat, clust_points, weights = NULL) {
  weights <- weights %||% rep(1,nrow(space_distmat))
  clust_weights <- weights[clust_points]
  
  # find the medoid # weighted
  clust_distmat <- space_distmat[clust_points, clust_points, drop = FALSE]
  distsum <- clust_distmat %*% clust_weights
  medoid <- clust_points[which.min(distsum)[1]]
  
  # compute costs to the medoid # unweighted
  space_costs <- (space_distmat[ , medoid])
  return(space_costs) 
}

# check_space_distmat ---------------------------------------------------------
check_space_distmat <- function(data, coords = c(1L, 2L), space_distmat, space_distmethod) {
  if(is.null(space_distmat)) {
    if (is.null(space_distmethod)) {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      message(paste0(space_distmethod," is used to compute space_distmat"))
    } else {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
    }
    space_distmat <- compute_distmat(data = data[coords], method = space_distmethod)
  } else {
    space_distmat <- as.matrix(space_distmat)
  }
  return(space_distmat)
}

# check_col_ref ---------------------------------------------------------------
check_col_ref <- function(x, data, arg_name, expected_length) {
  if (length(x) != expected_length) {
    stop(sprintf("`%s` must have length %d", arg_name, expected_length))
  }
  
  if (is.character(x)) {
    if (!all(x %in% names(data))) {
      stop(sprintf("`%s` must reference valid column names in `data`", arg_name))
    }
  } else if (is.numeric(x)) {
    if (!all(x >= 1L & x <= ncol(data))) {
      stop(sprintf("`%s` must reference valid column indices in `data`", arg_name))
    }
  } else {
    stop(sprintf("`%s` must be either character (column names) or numeric (column indices)", arg_name))
  }
  
  invisible(NULL)
}

# reorder_clust ---------------------------------------------------------------
reorder_clust <- function(clust) {
  # e.g. c(1,4,4,2) will become c(1,2,2,3)
  if(length(unique(clust)) == 1) return(rep(1, length(clust)))
  
  sets <- list()
  for (i in 1:length(clust)) {
    sets[[i]] <- which(clust == clust[i])
    if (length(sets[[i]]) == 0) sets[[i]] <- NA
  }
  
  max_set_length <- max(sapply(sets, length))
  sets <- lapply(sets, function(x) if (length(x) < max_set_length) {
    c(x, rep(0, max_set_length - length(x)))} else {x})
  sets <- do.call(rbind, sets)
  sets <- unique(sets)
  # remove empty clusters
  sets <- as.matrix(sets[rowSums(sets) != 0, ])
  
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

# find_dup --------------------------------------------------------------------
find_dup <- function (clust, ari = 1) {
  # as NA is ignored by mclust::adjustedRandIndex()
  clust[is.na(clust)] <- 0
  # number of solutions
  n_sol <- nrow(clust)
  # get all combination of pairs
  pairs <- utils::combn(n_sol, 2)
  
  # calculate the pairwise ARI for all pairs of solutions
  pairwise_ari <- future.apply::future_vapply(
    1:ncol(pairs),
    function(x) mclust::adjustedRandIndex(clust[pairs[1,x], ], clust[pairs[2,x], ]),
    numeric(1),
    future.seed = TRUE
  )
  
  # index the pairs with ari >= ari
  pairs_dup_idx <- which(pairwise_ari >= ari)
  # subset the columns of duplicated pairs
  pairs_dup <- pairs[, pairs_dup_idx, drop = FALSE]
  
  # if a ~ b and b ~ c, then not necessarily a ~ c 
  if (ncol(pairs_dup) > 1) {
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
    
    # remove the dependent columns 
    if (length(pairs_dup_depend_idx) > 0) {
      pairs_dup <- pairs_dup[ , -pairs_dup_depend_idx, drop = FALSE]
    }
  }
  
  # note it is sensitive to order
  # second row is considered the duplicate of the first row
  idx <- unique(pairs_dup[2, ])
  
  # frequency of duplicates
  freq <- as.data.frame(table(pairs_dup[1,]), stringsAsFactors = FALSE)
  # as.data.frame turn it either character or factor
  names(freq) <- c("idx", "n")
  freq$idx <- as.numeric(freq$idx)
  
  dup_out <- list(idx = idx,
                  pairwise_ari = pairwise_ari, # all pairs
                  freq = freq,
                  pairs_dup = pairs_dup)
  return(dup_out)
}

# globalVariables -------------------------------------------------------------
utils::globalVariables(
  names = c(
    ".data", "batch", "clust", "geometry", "idx" , "iter", "k",
    "pareto", "pareto_similar", "r", "run", "stat", "value","obj_value"
  ), package = "stblob"
)