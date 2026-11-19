#' Evaluate multi-objective optimisation (MOO)
#' 
#' @description
#' `eval_moo()` evaluates MOO.
#' 
#' @inheritParams stblob
#' @param pop an S3 `stblob_pop` object.
#' 
#' @inherit stblob return
#' 
#' @export

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
                   summary_pf = NULL,
                   clust = pop$clust, # NULL
                   summary = summary, # NULL
                   trace = trace, # NULL
                   moo_convergence = NULL,
                   moo_diagnostics = NULL,
                   filter_summary = pop$filter_summary,
                   data = pop$data,
                   space_distmat = space_distmat,
                   adapt_search_params = adapt_search_params,
                   params = pop$params))
  }
  
  # objective values and space
  z_cols <- c("z_space", "z_time1", "z_time2", "z_type")
  z_cols <- if (optim_type_diversity) z_cols else z_cols[1:3]
  # normalised objective space z
  z <- pop$summary[z_cols]; z <- as.data.frame(apply(z, 2, minmax))
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
  paretof_n[is.na(paretof_n)] <- 0
  
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
  
  # 3 consecutive front = 2 consecutive changes
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
                 clust = pop$clust,
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
