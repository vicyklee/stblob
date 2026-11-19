#------------------------------------------------------------------------------#
# Plotting data ####
#------------------------------------------------------------------------------#
#' Save a plot as PDF
#' 
#' @description
#' This function saves a plot as PDF. As [ggplot2::ggsave()] can be buggy with text spacing.
#' @param p plot.
#' @inheritParams grDevices::pdf
#' @inheritDotParams grDevices::pdf
#' @seealso [grDevices::pdf()]
#' @export

save_pdf <- function(p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = 'CP1250')
  print(p)
  invisible(grDevices::dev.off())
}

#------------------------------------------------------------------------------#
#' Plot objective space
#' 
#' @description
#' This function plots the objective space.
#' @param pop a list of objects returned from [stblob()].
#' @param obj a string vector of length 2 that indicates the objectives.
#' @param colour a string of the attribute to be coloured. It must be either "pareto", "r", "batch" or "k".
#' @param normalise a Boolen. Should values of objectives be normalised? Default is TRUE.
#' @param size size for [ggplot2::geom_point()].
#' @param alpha See <[`aes-colour-fill-alpha`][ggplot2::aes_colour_fill_alpha()]>.
#' @return a `ggplot` object.
#' @export

plot_objspace <- function (pop,
                           obj,
                           colour = c("r", "batch", "k", "pareto"),
                           normalise = TRUE,
                           size= 1,
                           alpha = 0.8) {
  
  # check for input
  colour <- match.arg(colour)
  if(!colour %in%  names(pop$summary)) {
    stop("colour option is not available for your input.")
  }
  
  parse_obj_out <- parse_obj(obj)
  signed_obj <- obj
  obj <- parse_obj_out$obj
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  
  objspace <- pop$summary %>%
    dplyr::select(c(dplyr::all_of(obj), dplyr::all_of(colour)))
  
  if (colour == "batch") {
    objspace <- objspace %>% dplyr::mutate(batch = as.factor(batch))
  }
  if (colour == "k") {
    objspace <- objspace %>% dplyr::mutate(k_o = as.factor(k))
  }
  if (colour == "pareto") {
    objspace <- objspace %>% dplyr::mutate(pareto = factor(pareto, levels=c("1","0")))
  }
  
  if(length(parse_obj_out$maximise_obj_idx) > 0) {
    objspace <- objspace %>%
      dplyr::mutate_at(obj[maximise_obj_idx], function(x) -1*x) %>%
      dplyr::rename_with(~paste0("-",.x), dplyr::all_of(obj[maximise_obj_idx]))
  }
  
  if (normalise == TRUE) {
    # obtain upper and lower bounds from the Pareto front
    ub <- vapply(1:2, function(i) max(objspace[[signed_obj[i]]]), numeric(1))
    lb <- vapply(1:2, function(i) min(objspace[[signed_obj[i]]]), numeric(1)) 
    # normalise the objectives
    objspace[, c(1,2)] <- moocore::normalise(objspace[, c(1,2)], to_range = c(0,1), lower = lb, upper = ub)
  }
  
  p <- ggplot2::ggplot(objspace) +
    ggplot2::geom_point(ggplot2::aes(x = .data[[signed_obj[1]]],
                                     y = .data[[signed_obj[2]]],
                                     colour = .data[[colour]]),
                        size = size, alpha = alpha)
  
  if (colour == "r") {
    p <- p + ggplot2::scale_color_viridis_c(option = "G", direction = -1, begin = 0.1, end = 0.8)
  }
  if (colour == "batch") {
    p <- p + ggplot2::scale_color_viridis_d(option = "H", direction = -1)
  }
  if (colour == "k") {
    p <- p + ggplot2::scale_color_viridis_d(option = "A", direction = -1, begin = 0.1, end = 0.8)
  }
  if (colour == "pareto") {
    p <- p + ggplot2::scale_colour_manual(values = c("#95c36e", "grey"))
  }
  
  if (normalise == TRUE) {
    p <- p + ggplot2::coord_fixed()
  }
  return(p)
}

#------------------------------------------------------------------------------#
#' Plot parallel coordinate plot
#' 
#' @description
#' This function plots the objective space.
#' @param pop a list of objects returned from [stblob()].
#' @param colour a string of the attribute to be coloured. It must be either "pareto", "r", "batch" or "k".
#' @param normalise a Boolen. Should values of objectives be normalised? Default is TRUE.
#' @param pareto a Boolen. Should only the Pareto solutions be displayed? Default is TRUE.
#' @param pareto_similar a Boolen. Should only the similar Pareto solutions be displayed? Default is FALSE.
#' @param linewidth linewidth for [ggplot2::geom_line()].
#' @param alpha See <[`aes-colour-fill-alpha`][ggplot2::aes_colour_fill_alpha()]>.
#' @return a `ggplot` object.
#' @export

plot_pcp <- function (pop,
                      colour = c("r","k","batch","pareto"),
                      normalise = TRUE,
                      pareto = TRUE,
                      pareto_similar = FALSE,
                      linewidth = 1,
                      alpha = 0.8) {
  
  # check for input
  colour <- match.arg(colour)
  if(!colour %in%  names(pop$summary)) {
    stop("colour option is not available for your input.")
  }
  
  obj <- pop$obj
  parse_obj_out <- parse_obj(obj)
  signed_obj <- obj
  obj <- parse_obj_out$obj
  maximise_obj_idx <- parse_obj_out$maximise_obj_idx
  
  objspace <- pop$summary %>%
    dplyr::select(c(dplyr::all_of(obj), dplyr::all_of(colour)), pareto, pareto_similar, idx)
  
  if (colour == "batch") {
    objspace <- objspace %>% dplyr::mutate(batch = as.factor(batch))
  }
  if (colour == "k") {
    objspace <- objspace %>% dplyr::mutate(k = as.factor(k))
  }
  if (colour == "pareto") {
    objspace <- objspace %>% dplyr::mutate(pareto = factor(pareto, levels=c("1","0")))
  }
  
  if (length(parse_obj_out$maximise_obj_idx) > 0) {
    objspace <- objspace %>%
      dplyr::mutate_at(obj[maximise_obj_idx], function(x) -1*x) %>%
      dplyr::rename_with(~paste0("-",.x), dplyr::all_of(obj[maximise_obj_idx]))
  }
  
  if (normalise == TRUE) {
    # obtain upper and lower bounds from the Pareto front
    ub <- vapply(1:length(obj), function(i) max(objspace[[signed_obj[i]]]), numeric(1))
    lb <- vapply(1:length(obj), function(i) min(objspace[[signed_obj[i]]]), numeric(1)) 
    # normalise the objectives
    objspace[, signed_obj] <- moocore::normalise(objspace[, signed_obj], to_range = c(0,1), lower = lb, upper = ub)
  }
  
  if (pareto == TRUE) {
    objspace <- objspace %>%
      dplyr::filter(pareto == 1)
  }
  
  if (pareto_similar == FALSE) {
    objspace <- objspace %>%
      dplyr::filter(pareto_similar == 0)
  }
  
  objspace_long <- tidyr::pivot_longer(objspace,
                                       cols = dplyr::all_of(signed_obj),
                                       names_to = "obj",
                                       values_to = "obj_value") %>%
    dplyr::mutate(obj = factor(obj, levels = signed_obj))
  
  p <- ggplot2::ggplot(objspace_long) +
    ggplot2::geom_line(ggplot2::aes(x = obj,
                                    y = obj_value,
                                    colour = .data[[colour]],
                                    group = idx),
                       linewidth = linewidth, alpha = alpha) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.1,0.1)))
  
  if (colour == "r") {
    p <- p + ggplot2::scale_color_viridis_c(option = "G", direction = -1, begin = 0.1, end = 0.8)
  }
  if (colour == "batch") {
    p <- p + ggplot2::scale_color_viridis_d(option = "H", direction = -1)
  }
  if (colour == "k") {
    p <- p + ggplot2::scale_color_viridis_d(option = "A", direction = -1, begin = 0.1, end = 0.8)
  }
  if (colour == "pareto") {
    p <- p + ggplot2::scale_colour_manual(values = c("#95c36e", "grey"))
  }
  
  return(p)
}

#------------------------------------------------------------------------------#
#' Pivot trace summary for plotting
#' 
#' @description
#' This function pivots trace summary for plotting.
#' @param data a data frame of trace summary.
#' @return a pivoted data frame of trace summary.
#' @seealso [tidyr::pivot_longer()]

pivot_trace <- function (data) {
  stat <- c("space_wcd", "time_wcr", "time_wce", "ari")
  
  optional_cols <- c("batch", "k")
  optional_cols <- optional_cols[optional_cols %in% names(data)]
  fct_cols <- c("run", optional_cols)

  data <- data %>%
    dplyr::select(c(dplyr::all_of(stat), "iter", "r", "run", dplyr::all_of(fct_cols))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(stat), names_to = "stat", values_to = "value") %>%
    dplyr::mutate(stat = factor(stat, levels = unique(stat)),
                  dplyr::across(dplyr::all_of(fct_cols), as.factor))

  return(data)
}

#------------------------------------------------------------------------------#
#' Plot trace
#' 
#' @description
#' This function plots the trace.
#' @inheritParams plot_objspace
#' @param colour a string of the attribute to be coloured. It must be either "r", "batch" or "k".
#' @return a ggplot.
#' @seealso [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @export

plot_trace <- function(pop, colour = c("r", "batch", "k"), alpha = 0.8) {
  # check for input
  colour <- match.arg(colour)
  
  optional_cols <- c("batch", "k")
  optional_cols <- optional_cols[optional_cols %in% names(pop$trace)]
  group_cols <- c("stat", "run", optional_cols)
  
  trace <- pop$trace %>% pivot_trace()
  
  p <- ggplot2::ggplot(trace) +
    ggplot2::geom_line(ggplot2::aes(x = iter, y = value,
                                    group = interaction(!!!rlang::syms(group_cols), sep = "_"),
                                    colour = .data[[colour]]), alpha = alpha) +
    ggplot2::facet_wrap(~stat, scales = "free", nrow = 1) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("iteration")
  
  if (colour == "r") {
    p <- p + ggplot2::scale_color_viridis_c(option = "G", direction = -1, begin = 0.1, end = 0.8)
  }
  if (colour == "batch") {
    p <- p + ggplot2::scale_color_viridis_d(option = "H", direction = -1)
  }
  if (colour == "k") {
    p <- p + ggplot2::scale_color_viridis_d(option = "A", direction = -1, begin = 0.1, end = 0.8)
  }
  
  return(p)
}

#------------------------------------------------------------------------------#
#' Plot data in space
#' 
#' @description
#' This function plots data in space.
#' @param data a data frame
#' @param clust a numeric vector of cluster assignment.
#' @param space a string of the space data plotted in. It must be either "earth" or "euclidean".
#' @param coords a vector of strings or numeric values indicating the columns of coordinates (longitude, latitide). Default is the first two columns.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()] and [sf::st_transform()]. Default is NULL.
#' @param hull a Boolean. Should convex hulls be drawn? Default is FALSE.
#' @param size size for [ggplot2::geom_point()].
#' @param alpha See <[`aes-colour-fill-alpha`][ggplot2::aes_colour_fill_alpha()]>.
#' @param weights a numeric vector of weights for each data point, used to indicate them by data point colours. Default is NULL.
#' @inheritParams intersects_bool
#' @return a `ggplot` object.
#' @seealso [sf::st_as_sf()], [sf::st_transform()]
#' @export

plot_space <- function(data,
                       clust = NA,
                       space = "earth",
                       coords = c(1,2),
                       crs = NULL,
                       size = 1,
                       alpha = 0.8,
                       hull = FALSE,
                       hull_convex_ratio = 0.8,
                       weights = NULL) {
  
  space <- match.arg(space, choices = c("earth", "euclidean"))
  if (space == "earth" & is.null(crs)) crs <- 4326
  if (space == "euclidean" & is.null(crs)) crs <- NA
  
  data$clust <- as.factor(clust)
  pts <- sf::st_as_sf(data, coords = coords, crs = crs)
  if (!is.na(crs)) {
    pts <- sf::st_transform(pts, crs = crs)
    # Get bounding box of your data
    bb <- sf::st_bbox(pts)
    # Add a buffer (in degrees, assuming EPSG:4326)
    if (crs == 4036 | crs == 4326) {
      buffer <- 5
    }
    if (crs == 3035) {
      buffer <- 500000 # 500 km buffer for 3035
    }
    xlim <- c(bb["xmin"] - buffer, bb["xmax"] + buffer)
    ylim <- c(bb["ymin"] - buffer, bb["ymax"] + buffer)
  }
  
  begin <- 0.1
  end <- 0.8
  
  if (is.null(weights)) {
    colour <- "clust"
  } else {
    pts$weight <- weights
    colour <- "weight"
  }
  
  if (hull == FALSE) {
    if (space == "euclidean") {
      if (!all(is.na(clust))) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = pts, size = size, alpha = alpha, ggplot2::aes(colour = clust))
      } else {
        p <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = pts, size = size, alpha = alpha, colour = "blue")
      }
      
    } else if (space == "earth") {
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      world <- sf::st_transform(world, crs = crs)
      if (!all(is.na(clust))) {
        p <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = world, alpha = 0.2) +
          ggplot2::geom_sf(data = pts, size = size, alpha = alpha, ggplot2::aes(colour = clust)) +
          ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
      } else {
        p <- ggplot2::ggplot() +
          ggplot2::geom_sf(data = world, alpha = 0.2) +
          ggplot2::geom_sf(data = pts, size = size, alpha = alpha, colour = "blue") +
          ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE)
      }
    }
  } else {
    if (all(is.na(clust))) {
      stop("clust is required to draw hulls.")
    }
    suppressWarnings({
      hulls <- pts %>%
        dplyr::filter(!is.na(clust)) %>%
        dplyr::group_by(clust) %>%
        dplyr::summarise(geometry = sf::st_combine(geometry)) %>% 
        sf::st_concave_hull(ratio = hull_convex_ratio) %>%
        sf::st_make_valid() %>%
        sf::st_collection_extract("POLYGON", warn = FALSE)
    })
    
    if (space == "euclidean") {
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = hulls, ggplot2::aes(fill = clust), colour = NA, lwd = 0.3, alpha = 0.3, show.legend = F) +
        ggplot2::geom_sf(data = pts, size = size, alpha = 0.8, ggplot2::aes(colour = .data[[colour]])) 
      
    } else if (space == "earth") {
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      world <- sf::st_transform(world, crs = crs)
      
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = world, alpha = 0.2) +
        ggplot2::geom_sf(data = hulls, ggplot2::aes(fill = clust), colour = NA, lwd = 0.3, alpha = 0.3, show.legend = F) +
        ggplot2::geom_sf(data = pts, size = size, alpha = 0.8, ggplot2::aes(colour = clust)) +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = crs)
    }
    p <- p +
      ggplot2::scale_fill_viridis_d(option = "A", direction = -1, begin = begin, end = end)
  }
  
  if (is.null(weights)) {
    if(!all(is.na(clust))) {
      p <- p + ggplot2::scale_colour_viridis_d(option = "A", direction = -1, begin = begin, end = end, na.value = "grey50")
    }
  } else {
    p <- p + ggplot2::scale_colour_viridis_c(option = "A", direction = -1, begin = begin, end = end)
  }
  return(p)
}

#------------------------------------------------------------------------------#
#' Plot data in time
#' 
#' @description
#' This function plots data in time
#' @param data a data frame
#' @param clust a numeric vector of cluster assignment.
#' @param age a string or numeric value indicating the column of age. Default is the third column. 
#' @export

plot_time <- function (data, clust, age = 3) {
  data$clust <- as.factor(clust)
  if (is.numeric(age)){
    age_col <- names(data)[age]
  } else {
    age_col <- age
  }
  
  begin <- 0.1
  end <- 0.8
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[age_col]],
                                          y = clust,
                                          colour = clust,
                                          fill = clust)) +
    ggridges::geom_density_ridges(alpha = 0.2, scale = 0.5, 
                                  jittered_points = T,
                                  position = ggridges::position_points_jitter(width = 0, height = 0),
                                  point_shape = '|', point_size = 3, point_alpha = 0.8) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_colour_viridis_d(option = "A", direction = -1, begin = begin, end = end, na.value = "grey50") +
    ggplot2::scale_fill_viridis_d(option = "A", direction = -1, begin = begin, end = end) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = c(0,0), add = c(0.2,0.6)))
  return(p)
}

#------------------------------------------------------------------------------#
#' Plot multi-objective optimisation quality
#' 
#' @description
#' This function plots multi-objective optimisation (MOO) quality.
#' @param pop_moo a `pop_moo` object; see also [stblob()].
#' @param indicator a string of MOO quality indicator. It must be either "igd", "igd_plus" or "hv".
#' @returns a `ggplot` object.
#' @export

plot_mooquality <- function(pop_moo, indicator = c("igd", "igd_plus", "hv")) {
  indicator <- match.arg(indicator)
  moo_quality <- cbind(batch = as.integer(1:nrow(pop_moo$moo_quality)), pop_moo$moo_quality)
  ggplot2::ggplot(moo_quality) +
    ggplot2::geom_line(ggplot2::aes(x = batch, y = .data[[indicator]]), colour = "#5a97c1")
    # ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
}
