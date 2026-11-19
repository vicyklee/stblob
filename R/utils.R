# compute_space_costs ---------------------------------------------------------
compute_space_costs <- function(space_distmat, clust_points) {
  weights <- rep(1,nrow(space_distmat)) # deprecated weights option
  clust_weights <- weights[clust_points]
  
  # find the medoid # weighted
  clust_distmat <- space_distmat[clust_points, clust_points, drop = FALSE]
  distsum <- clust_distmat %*% clust_weights
  medoid <- clust_points[which.min(distsum)[1]]
  
  # compute costs to the medoid # unweighted
  space_costs <- (space_distmat[ , medoid])
  return(space_costs) 
}

# check_data_cols -------------------------------------------------------------
check_args <- function(data,
                       optim_type_diversity,
                       iter,
                       coords,
                       age,
                       type) {
  
  stopifnot(
    "'data' must be a data frame" =
      is.data.frame(data),
    "'optim_type_diversity' must be logical of length 1" =
      is.logical(optim_type_diversity) && length(optim_type_diversity) == 1,
    "'iter' must be positive numeric of length 1" =
      is.numeric(iter) && length(iter) == 1 && iter >= 1
  )
  
  # check the length of the input
  args <- list(coords, age, type)
  arg_defaults <- list(c(2,3), 4, 5)
  arg_names <- c("coords", "age", "type")
  arg_lens <- c(2L,1L,1L)
  
  IDX <- numeric()
  for (i in 1:3) {
    # skip type check if !optim_type_diversity
    if (!optim_type_diversity) {
      if (arg_names[i] == "type") next
    }
    
    # check arg length
    if (!(length(args[[i]]) == arg_lens[i] || is.null(args[[i]]))) {
      stop(sprintf("'%s' must be of length %i", arg_names[i], arg_lens[i]))
    }
    
    # check if args match data columns
    idx <- if (!is.null(args[[i]])) match(args[[i]], names(data)) else arg_defaults[[i]]
    if (any(is.na(idx))) {
      stop(sprintf("'%s' does not match any column names", arg_names[[i]]))
    }
    IDX <- append(IDX, idx)
    
    # check arg types for coords and age
    if (arg_names[i] == "coords") stopifnot("'coords' column must be numeric" = is.numeric(data[[ idx[1] ]]) && is.numeric(data[[ idx[2] ]]))
    if (arg_names[i] == "age") stopifnot("'age' column must be numeric" = is.numeric(data[[idx]]))
  }
  
  # check if any column idx used twice
  stopifnot("Same 'data' column used twice" = all(table(IDX) == 1))
  
  invisible(NULL)
}

# check_space_distmat ---------------------------------------------------------
check_space_distmat <- function(data,
                                coords,
                                space_distmat,
                                space_distmethod,
                                distGeo_a,
                                distGeo_f) {
  if(is.null(space_distmat)) {
    # check space_distmethod
    if (is.null(space_distmethod)) {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
      message(paste0(space_distmethod," is used to compute space_distmat"))
    } else {
      space_distmethod <- match.arg(space_distmethod, choices = c("geodesic", "euclidean"))
    }
    
    space_distmat <- compute_distmat(data = data[coords],
                                     method = space_distmethod,
                                     distGeo_a = distGeo_a,
                                     distGeo_f = distGeo_f)
    
  } else {
    # check space_distmat
    stopifnot(
      "'space_distmat' must be a numeric matrix or an S3 'dist' object" =
        (is.numeric(space_distmat) && is.matrix(space_distmat)) || inherits(space_distmat, "dist"),
      "The dimension of 'space_distmat' does not correspond to the number of rows of 'data'" =
        nrow(data) == nrow(space_distmat) && nrow(data) == nrow(space_distmat)
    )
    
    space_distmat <- as.matrix(space_distmat)
    
  }
  return(space_distmat)
}

# shannon_entropy -------------------------------------------------------------
shannon_entropy <- function(x, n = length(unique(x)), base = 2, normalise = TRUE) {
  p <- as.numeric(table(x)/length(x))
  I <- log(p, base)
  I[p==0] <- 0
  H <- -sum(p*I)
  if (normalise == TRUE) H <- H/log(n, base)
  return(H)
}

# check_opt_args --------------------------------------------------------------
check_opt_args <- function(ls_tol = 0,
                           chebyshev_rho = 1e-4, 
                           filter_intersects = TRUE,
                           hull_convex_ratio = 0.5,
                           hull_crs = 4326,
                           filter_clustsize = TRUE, 
                           outlier_removal = FALSE,
                           outlier_iqrm = 1.5,
                           random_init = FALSE,
                           sf_use_s2 = TRUE,
                           ...) {
  
  stopifnot(
    "'ls_tol' must be between 0 and 1" =
      is.numeric(ls_tol) && ls_tol >= 0 && ls_tol <= 1,
    "'chebyshev_rho' must be non-negative numeric" =
      is.numeric(chebyshev_rho) && chebyshev_rho >= 0,
    "'filter_intersects' must be logical" =
      is.logical(filter_intersects),
    "'hull_convex_ratio' must be between 0 and 1" =
      is.numeric(hull_convex_ratio) && hull_convex_ratio >= 0 && hull_convex_ratio <= 1,
    "'hull_crs' cannot retrieve coordinate reference system" =
      is.na(hull_crs) | !is.na(sf::st_crs(hull_crs)$wkt),
    "'filter_clustsize' must be logical" =
      is.logical(filter_clustsize),
    "'outlier_removal' must be logical" =
      is.logical(outlier_removal),
    "'outlier_iqrm' must be non-negative numeric" =
      is.numeric(outlier_iqrm) && outlier_iqrm >= 0,
    "'random_init' must be logical" =
      is.logical(random_init),
    "'sf_use_s2' must be logical" = is.logical(sf_use_s2)
  )
  
  invisible(NULL)
}

# reorder_clust ---------------------------------------------------------------
reorder_clust <- function(clust) {
  clust <- as.factor(clust)
  levels(clust) <- 1:length(levels(clust))
  clust <- as.integer(clust)
  return(clust)
}

# minmax ----------------------------------------------------------------------
minmax <- function(x, lb = NULL, ub = NULL, max = FALSE, na.rm = TRUE) {
  # lower bound
  lb <- min(x, lb, na.rm = na.rm) %||% min(x, na.rm = na.rm)
  # upper bound
  ub <- max(x, ub, na.rm = na.rm) %||% max(x, na.rm = na.rm)
  
  # normalise
  if (!max) {
    x <- (x - lb) / (ub - lb)
    x <- if (!all(is.na(x))) x else rep(0, length(x))
  } else {
    x <- (x - ub) / (lb - ub)
    x <- if (!all(is.na(x))) x else rep(1, length(x))
  } 
  return(x)
}

# globalVariables -------------------------------------------------------------
utils::globalVariables(names = c("idx", "idx_pf"))
# utils::globalVariables(
#   names = c(
#     ".data", "batch", "clust", "geometry", "idx" , "iter", "k",
#     "pareto", "pareto_similar", "r", "run", "stat", "value","obj_value"
#   ), package = "stblob"
# )