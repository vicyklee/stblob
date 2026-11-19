# print.stblob_sol ------------------------------------------------------------
#' @export
print.stblob_sol <- function(x, ...) {
  if(x$status == 0) {
    cat("---------------------------------------------------\n")
    cat("STblob solution\n")
    cat("---------------------------------------------------\n")
    
    cat("$clust\n# Cluster assignments:\n")
    print(utils::head(x$clust, 20))
    if (length(x$clust) > 20) cat(paste0(" <", length(x$clust) - 20," more points>\n"))
    cat("\n")
    cat("$summary\n# Summary:\n")
    print(x$summary)
  } else {
    cat("---------------------------------------------------\n")
    cat("STblob solution\n")
    cat("---------------------------------------------------\n")
    cat("* No feasible solution was found *\n")
    cat("\n")
    cat(paste0("Status code: ", x$status, "\n"))
  }
  invisible(x)
}

# summary.stblob_sol-----------------------------------------------------------
#' @export
summary.stblob_sol <- function(object, ...) {
  print(object$summary)
  invisible(object)
}

# print.stblob_moo ------------------------------------------------------------
#' @export
print.stblob_moo <- function(x, ...) {
  if (!is.null(x$clust)) {
    if (is.matrix(x$clust_pf)) {
      cat("---------------------------------------------------\n")
      cat("STblob MOO output of", nrow(x$clust_pf), "solutions\n")
    } else {
      cat("---------------------------------------------------\n")
      cat("STblob MOO output of 1 solution\n")
    }
    cat("---------------------------------------------------\n")
    
    cat("$clust_pf\n# Cluster assignments on the Pareto front:\n")
    if (is.matrix(x$clust_pf)) {
      print(utils::head(x$clust_pf, c(5, 20)))
      if (nrow(x$clust_pf) > 5) cat(paste0(" <", nrow(x$clust_pf) - 5," more rows (solutions)>\n"))
      if (ncol(x$clust_pf) > 20) cat(paste0(" <", ncol(x$clust_pf) - 20," more columns (points)>\n"))
    } else {
      print(utils::head(x$clust, 20))
      if (length(x$clust) > 20) cat(paste0(" <", length(x$clust) - 20," more points>\n"))
    }
    cat("\n")
    
    cat("$summary\n# Summary:\n")
    print(utils::head(x$summary, 5))
    if (nrow(x$summary) > 5) cat(paste0(" <", nrow(x$summary) - 5," more rows>\n"))
    cat("\n")
    
    cat("$moo_convergence\n# MOO convergence status:\n")
    cat(paste0(x$moo_convergence, " (moo_tol: ", x$params$moo_tol, ")\n"))
    cat("\n")
    
    cat("$moo_diagnostics\n# MOO diagnostics:\n")
    print(utils::head(x$moo_diagnostics, c(5, ncol(x$moo_diagnostics))))
    if (nrow(x$moo_diagnostics) > 5) cat(paste0(" <", nrow(x$moo_diagnostics) - 5," more batches>\n"))
    cat("\n")
    
    cat("$filter_summary\n# Number of passed and filtered solutions: \n")
    print(x$filter_summary[c("k",
                             "pass",
                             "k1",
                             "intersects",
                             "clustsize",
                             "intersects_clustsize",
                             "dup")])
    
  } else {
    cat("---------------------------------------------------\n")
    cat("STblob MOO output of 0 solutions\n")
    cat("---------------------------------------------------\n")
    cat("* No feasible solution was found *\n")
    cat("\n")
    
    cat("$filter_summary\n# Number of passed and filtered solutions:\n")
    print(x$filter_summary[c("k",
                             "pass",
                             "k1",
                             "intersects",
                             "clustsize",
                             "intersects_clustsize",
                             "dup")])
  }
  
  invisible(x)
}

# summary.stblob_moo ----------------------------------------------------------
#' @export
summary.stblob_moo <- function(object, ...) {
  print(object$summary)
  invisible(object)
}
