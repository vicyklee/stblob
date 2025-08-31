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

savePDF <- function(p, file, width, height, ...) {
  grDevices::pdf(file = file, width = width, height = height, ...)
  grDevices::pdf.options(encoding = 'CP1250')
  print(p)
  invisible(grDevices::dev.off())
}

#' Plot objective space
#' 
#' @description
#' This function plots the objective space.
#' @param pop a list of objects returned from [blob_moo()].
#' @param obj a string vector of length 2 that indicate the objectives.
#' @param colour a string of the attribute to be coloured. It must be "pareto", "r", "batch" or "k_o".
#' @param normalise a logical operator. Should values of objectives be normalised? Default is T.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @param alpha See <[`aes-colour-fill-alpha`][ggplot2::aes_colour_fill_alpha()]>.
#' @return a ggplot.
#' @seealso [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @export

plot_objspace <- function (pop, obj, colour = c("r", "batch", "k_o", "pareto"), normalise = T, palette = NULL, alpha = 0.8) {
  colour <- match.arg(colour)

  objspace <- pop$summary %>%
    dplyr::select(c(dplyr::all_of(obj), "pareto", "r", "batch", "k_o")) %>%
    dplyr::mutate(pareto = as.factor(pareto),
                  batch = as.factor(batch),
                  k_o = as.factor(k_o))
  pareto_objspace <- subset(objspace, pareto == 1)
  
  if (normalise == T) {
    # obtain upper and lower bounds from the Pareto front
    ub <- vapply(1:2, function(i) max(pareto_objspace[[obj[i]]]), numeric(1))
    lb <- vapply(1:2, function(i) min(pareto_objspace[[obj[i]]]), numeric(1)) 
    # normalise the objectives
    objspace[, c(1,2)] <- moocore::normalise(objspace[, c(1,2)], to_range = c(0,1), lower = lb, upper = ub)
  }
  
  p <- ggplot2::ggplot(objspace) +
    ggplot2::geom_point(ggplot2::aes(x = .data[[obj[1]]], y = .data[[obj[2]]], colour = .data[[colour]]), alpha = alpha)
  
  if (!is.null(palette)) {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = palette) 
    }
    if (colour == "batch" | colour == "k_o" | colour == "pareto") {
      p <- p + ggplot2::scale_colour_manual(values = palette,
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
  } else {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = MetBrewer::met.brewer("Hokusai3"))
    }
    if (colour == "batch") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Egypt", n = length(unique(objspace$batch))),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
    if (colour == "k_o") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Archambault", n = length(unique(objspace$k_o))),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
    if (colour == "pareto") {
      p <- p + ggplot2::scale_colour_manual(values = c("grey","#95c36e"),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
  }
  
  if (normalise == T) {
    p <- p + ggplot2::coord_fixed()
  }
  return(p)
}

#' Plot data in space
#' 
#' @description
#' This function plots data in space.
#' @param data a data frame
#' @param clust a numeric vector of cluster assignemnt. Default is NA.
#' @param space a string of the space data plotted in. It must be either "earth" or "euclidean".
#' @param hull a logical operator. Should convex hulls be drawn? Default is F.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()] and [sf::st_transform()]. Default is NULL
#' @param lab a string vector of length 2. The first one is passed on to [ggplot2::xlab()] and the second [ggplot2::ylab()]. Default is NULL.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @return a ggplot.
#' @seealso [sf::st_as_sf()], [sf::st_transform()], [ggplot2::labs()], [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @export

plot_space <- function (data, clust = NA, space, hull = F, crs = NULL, lab = NULL, palette = NULL) {
  space <- match.arg(space, choices = c("earth", "euclidean"))
  if (space == "earth" & is.null(crs)) crs <- 4036
  if (space == "euclidean" & is.null(crs)) crs <- NA
  
  clust_levels <- levels(as.factor(clust))
  if (is.null(palette)) {
    clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels)) 
  }
  
  data$clust <- as.factor(clust)
  pts <- sf::st_as_sf(data, coords = c(1,2), crs = crs)
  if (!is.na(crs)) {
    pts <- sf::st_transform(pts, crs = crs)
    # Get bounding box of your data
    bb <- sf::st_bbox(pts)
    # Add a buffer (in degrees, assuming EPSG:4326)
    if (crs == 4036) {
      buffer <- 5
    }
    if (crs == 3035) {
      buffer <- 500000 # 500 km buffer for 3035
    }
    xlim <- c(bb["xmin"] - buffer, bb["xmax"] + buffer)
    ylim <- c(bb["ymin"] - buffer, bb["ymax"] + buffer)
  }
  
  if (hull == F) {
    if (space == "euclidean") {
      if (is.null(lab)) lab <- c("x", "y")
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
        ggplot2::scale_colour_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::xlab(lab[1]) +
        ggplot2::ylab(lab[2])
    } else if (space == "earth") {
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      world <- sf::st_transform(world, crs = crs)
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = world, alpha = 0.2) +
        ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::scale_colour_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::xlab(lab[1]) +
        ggplot2::ylab(lab[2])
    }
  } else if (hull == T) {
    hulls <- pts %>%
      dplyr::filter(!is.na(clust)) %>%
      dplyr::group_by(clust) %>%
      dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
      sf::st_convex_hull()
    if (space == "euclidean") {
      if (is.null(lab)) lab <- c("x", "y")
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = hulls, ggplot2::aes(colour = clust, fill = clust), lwd = 0.3, alpha = 0.1, show.legend = F) +
        ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
        ggplot2::scale_colour_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::scale_fill_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::xlab(lab[1]) +
        ggplot2::ylab(lab[2])
    } else if (space == "earth") {
      world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
      world <- sf::st_transform(world, crs = crs)
      p <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = world, alpha = 0.2) +
        ggplot2::geom_sf(data = hulls, ggplot2::aes(colour = clust, fill = clust), lwd = 0.3, alpha = 0.1, show.legend = F) +
        ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 3035) +
        ggplot2::scale_colour_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::scale_fill_manual(values = stats::setNames(clust_cols, clust_levels)) +
        ggplot2::xlab(lab[1]) +
        ggplot2::ylab(lab[2])
      
    }
  }
  return(p)
}

#' Plot data in time
#' 
#' @description
#' This function plots data in time
#' @param data a data frame
#' @param clust a numeric vector of cluster assignemnt.
#' @param lab a string of label passed on to [ggplot2::xlab()]. Default is NULL.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @return a ggplot.
#' @seealso [ggplot2::labs()], [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @export

plot_time <- function (data, clust, lab = NULL, palette = NULL) {
  clust_levels <- levels(as.factor(clust))
  if (is.null(palette)) {
    clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels)) 
  }
  
  data$clust <- as.factor(clust)
  
  if (is.null(lab)) lab <- "time"
  
  p <- ggplot2::ggplot(data, ggplot2::aes(y = clust, x = data[,3], colour = clust, fill = clust)) +
    ggridges::geom_density_ridges(alpha = 0.2, scale = 0.5,
                                  jittered_points = T,
                                  position = ggridges::position_points_jitter(width = 0, height = 0),
                                  point_shape = '|', point_size = 3, point_alpha = 0.8) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()) +
    ggplot2::scale_colour_manual(values = stats::setNames(clust_cols, clust_levels)) +
    ggplot2::scale_fill_manual(values = stats::setNames(clust_cols, clust_levels)) +
    ggplot2::xlab(lab)
  
  return(p)
}

#' Plot trace summary
#' 
#' @description
#' This function plots 2D objective space.
#' @inheritParams plot_objspace
#' @param colour a string of the attribute to be coloured. It must be "r", "batch" or "k_o".
#' @return a ggplot.
#' @seealso [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @export

plot_trace <- function(pop, alpha = 0.8, colour = c("r", "batch", "k_o"), palette = NULL) {
  colour <- match.arg(colour)

  trace <- pop$trace %>% pivot_trace()
  p <- ggplot2::ggplot(trace) +
    ggplot2::geom_line(ggplot2::aes(x = iter, y = value, group = interaction(stat, run, batch, r, k_o),
                                    colour = .data[[colour]]), alpha = alpha) +
    ggplot2::facet_wrap(~stat, scales = "free", nrow = 2) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())  +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = max(unique(trace$iter))))
  
  if (!is.null(palette)) {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = palette) 
    }
    if (colour == "batch" | colour == "k_o") {
      p <- p + ggplot2::scale_colour_manual(values = palette)
    }
  } else {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = MetBrewer::met.brewer("Hokusai3")) 
    }
    if (colour == "batch") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Egypt", n = length(unique(pop$trace$batch))))
    }
    if (colour == "k_o") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Archambault", n = length(unique(pop$trace$k_o))))
    }
  }
  return(p)
}

#' Pivot trace summary for plotting
#' 
#' @description
#' This function pivots trace summary for plotting.
#' @param df a data frame of trace summary.
#' @return a pivoted data frame of trace summary.
#' @seealso [tidyr::pivot_longer()]

pivot_trace <- function (df) {
  stat <- c("space_wcss", "time_range_mean", "time_range_sd", "time_evenness_mean", "time_evenness_sd", "ari")
  df <- df %>%
    dplyr::select(c(dplyr::all_of(stat), "iter", "r", "run", "batch", "k_o")) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(stat), names_to = "stat", values_to = "value") %>%
    dplyr::mutate(stat = factor(stat, levels = unique(stat)),
                  run = as.factor(run),
                  batch = as.factor(batch),
                  k_o = as.factor(k_o)) 
  return(df)
}

#' Plot multi-objective optimisation quality
#' 
#' @description
#' This function plots multi-objective optimisation (MOO) quality.
#' @inheritParams plot_objspace
#' @param indicator a string of MOO quality indicator. It must be either "igd", "igd_plus" or "hv".
#' @return a ggplot.
#' @export

plot_moo_quality <- function(pop, indicator = c("igd", "igd_plus", "hv")) {
  indicator <- match.arg(indicator)
  moo_quality <- cbind(batch = as.integer(1:nrow(pop$moo_quality)), pop$moo_quality)
  ggplot2::ggplot(moo_quality) +
    ggplot2::geom_line(ggplot2::aes(x = batch, y = .data[[indicator]]), colour = "#5a97c1") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
}
