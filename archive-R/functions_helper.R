#------------------------------------------------------------------------------#
# Helper functions ####
#------------------------------------------------------------------------------#
# Reorder clust vector in ascending order
reorder_clust <- function(clust) {
  # e.g. c(1,4,4,2) will become c(1,2,2,3)
  if(length(unique(clust)) == 1) return(rep(1, length(clust)))
  #-------------------------------------------------------------------#
  sets <- list()
  for (i in 1:length(clust)) {
    sets[[i]] <- which(clust == clust[i])
    if (length(sets[[i]]) == 0) sets[[i]] <- NA
  }
  #-------------------------------------------------------------------#
  max_set_length <- max(sapply(sets, length))
  sets <- lapply(sets, function(x) if (length(x) < max_set_length) {c(x, rep(0, max_set_length - length(x)))} else {x})
  sets <- do.call(rbind, sets)
  sets <- unique(sets)
  sets <- as.matrix(sets[rowSums(sets) != 0, ]) # remove empty clusters
  #-------------------------------------------------------------------#
  # put NA to the back
  # count the number of rows with NA
  n_row <- nrow(sets)
  n_row_na <- length(unique(which(is.na(sets), arr.ind = T)[,"row"]))
  if (n_row_na > 0) {
    sets <- as.matrix(stats::na.omit(sets))
    sets <- rbind(sets,NA)
  }
  #-------------------------------------------------------------------#
  # reorder cluster labels
  for (i in 1:nrow(sets)) {
    for (j in 1:ncol(sets)) {
      clust[sets[i,j]] <- i
    }
  }
  #-------------------------------------------------------------------#
  return(clust)
}
#------------------------------------------------------------------------------#
# Check space_distmat
check_space_distmat <- function(space_distmat, space_distmat) {
  if(is.null(space_distmat)) {
    if (is.null(space_distmethod)) {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      message(paste0(space_distmethod," is used to compute space_distmat"))
    } else {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
    }
    space_distmat <- compute_distmat(data = data[, coords], method = space_distmethod)
  }
  return(space_distmat)
}
#------------------------------------------------------------------------------#
# Select columns
