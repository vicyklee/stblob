#' Filter solutions on the Pareto front
#' 
#' @description
#' `filter_paretof()` filters for diverse solutions on the Pareto front.
#'
#' @param moo an S3 `stblob_moo` object.
#' @param ari a numeric value of the Adjusted Rand Index (ARI)
#' that sets the similarity threshold between two solutions on the Pareto front.
#' Default is `0.8`.
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

filter_paretof <- function(moo, ari = 0.8, z_dist = 0) {
  stopifnot(
    "'ari' must be equal to or smaller than 1" =
      (is.numeric(ari) && ari <= 1),
    "'z_dist' must be between 0 and 1" =
      (is.numeric(z_dist) && z_dist >= 0 && z_dist <= 1),
    "'moo' must be an S3 `stblob_moo` object" =
      inherits(moo, "stblob_moo")
  )
  
  clust_pf <- moo$clust_pf
  
  if (is.null(clust_pf)) {
    message("No solution to filter")
    return(NULL)
  }
  
  if (!is.null(clust_pf) && !is.matrix(clust_pf)) {
    message("Only one soluton on the Pareto front: returning 'idx_pf' of 1")
    return(1)
  }

  summary <- moo$summary
  clust_pf <- moo$clust_pf
  idx_pf <- moo$summary$idx_pf[!is.na(moo$summary$idx_pf)]
  optim_type_diversity <- unlist(moo$params$optim_type_diversity)[1]
  m <- if (optim_type_diversity) 4L else 3L
  z_cols <- c("z_space", "z_time1", "z_time2", "z_type"); z_cols <- z_cols[1:m]
  
  # paretof objective space and clust_pf
  z <- summary[z_cols]; z <- apply(z, 2, minmax)
  paretof_z <- z[idx_pf, ]
  clust_pf[is.na(clust_pf)] <- 0L
  
  # initialise vectors for the loop
  ref <- 1
  l <- length(idx_pf)
  sel <- logical(l)
  # the highest rank is always selected
  sel[1] <- TRUE
  
  for (i in 2:l) {
    if (ref == l) break
    # decision space
    a <- aricode::ARI(clust_pf[ref, ], clust_pf[i, ])
    # objective space
    # / sqrt(1*m) to get normalised Euclidean distance
    d <- as.numeric(stats::dist(paretof_z[c(ref,i), ])) / sqrt(1*m)
    
    # if distinct, keep i and compare i to i+1
    if (a < ari && d >= z_dist) {
      sel[i] <- TRUE
      ref <- i
    }
  }
  
  idx_pf <- idx_pf[sel]
  return(idx_pf)
}
