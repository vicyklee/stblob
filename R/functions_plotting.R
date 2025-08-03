#------------------------------------------------------------------------------#
# Plotting data ####
#------------------------------------------------------------------------------#
#' Save a plot as PDF
#' 
#' @description
#' This function saves a plot as PDF. As [ggplot2::ggsave()] can be buggy with text spacing.
#' @param p plot.
#' @inheritDotParams grDevices::pdf file width height
#' @param ... these include all additional arguments in [grDevices::pdf()].
#' @seealso [grDevices::pdf()]
#' @export

savePDF <- function(p, file, width, height, ...) {
  pdf(file = file, width = width, height = height, ...)
  pdf.options(encoding = 'CP1250')
  print(p)
  invisible(dev.off())
}

#' Plot 2D objective space
#' 
#' @description
#' This function plots 2D objective space.
#' @param pop a list of objects returned from [blob_moo()].
#' @param obj a string vector of length 2 that indicate the objectives.
#' @param colour a string of the attribute to be coloured. It must be "pareto", "r", "batch" or "k_o".
#' @param normalise a logical operator. Should values of objectives be normalised? Default is T.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @param alpha See <[`aes-colour-fill-alpha`][ggplot2::aes_colour_fill_alpha()]>.
#' @return a ggplot.
#' @seealso [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @importFrom magrittr %>%
#' @export

plot_2d_objspace <- function (pop, obj, colour, normalise = T, palette = NULL, alpha = 0.8) {
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
#' @param space a string of the space data plotted in. It must be either "earth" or "euclidean". Default is "earth".
#' @param hull a logical operator. Should convex hulls be drawn? Default is F.
#' @param crs a numeric value of the Coordinate Reference System passed on to [sf::st_as_sf()] and [sf::st_transform()]. Default is 3035.
#' @param buffer a numeric value of the buffer between the data points and edges of the plot. Default is 500000 for `crs = 3035`.
#' @param lab a string vector of length 2. The first one is passed on to [ggplot2::xlab()] and the second [ggplot2::ylab()]. Default is NULL.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @return a ggplot.
#' @seealso [sf::st_as_sf()], [sf::st_transform()], [ggplot2::labs()], [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @importFrom magrittr %>%
#' @export

plot_space <- function (data, clust = NA, space = "earth", hull = F, crs = 3035, buffer = 500000, lab = NULL, palette = NULL) {
  
  clust_levels <- levels(as.factor(clust))
  if (is.null(palette)) {
    clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels)) 
  }
  
  if (space == "euclidean") {
    if (is.null(lab)) lab <- c("x", "y")
    
    data$clust <- as.factor(clust)
    p <- ggplot2::ggplot(data, ggplot2::aes(x = data[,1], y = data[,2], colour = clust, fill = clust)) +
      ggplot2::geom_point(alpha = 0.5, show.legend = F) +
      ggplot2::scale_colour_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_fill_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
      ggplot2::coord_fixed() +
      ggplot2::xlab(lab[1]) +
      ggplot2::ylab(lab[2])
  }
  
  if (space == "earth" & hull == F) {
    if (is.null(lab)) lab <- c("Longitude", "Latitude")
    
    data$clust <- as.factor(clust)
    pts <- sf::st_as_sf(data, coords = c(1,2), crs = 4326)
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
    
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs = crs)
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world, alpha = 0.2) +
      ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::scale_colour_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_fill_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::xlab(lab[1]) +
      ggplot2::ylab(lab[2])
  }
  
  if (space == "earth" & hull == T) {
    if (is.null(lab)) lab <- c("Longitude", "Latitude")
    
    data$clust <- as.factor(clust)
    pts <- sf::st_as_sf(data, coords = c(1,2), crs = 4326)
    pts <- sf::st_transform(pts, crs = crs)
    # Get bounding box of your data
    bb <- sf::st_bbox(pts)
    # Add a buffer
    xlim <- c(bb["xmin"] - buffer, bb["xmax"] + buffer)
    ylim <- c(bb["ymin"] - buffer, bb["ymax"] + buffer)
    
    # Hull
    hulls <- pts %>%
      dplyr::filter(!is.na(clust)) %>%
      dplyr::group_by(clust) %>%
      dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
      sf::st_convex_hull()
    
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs = crs)
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world, alpha = 0.2) +
      ggplot2::geom_sf(data = hulls, ggplot2::aes(colour = clust, fill = clust), lwd = 0.3, alpha = 0.1, show.legend = F) +
      ggplot2::geom_sf(data = pts, alpha = 0.8, ggplot2::aes(colour = clust)) +
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 3035) +
      ggplot2::scale_colour_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_fill_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::xlab(lab[1]) +
      ggplot2::ylab(lab[2])
  }
  return(p)
}

#' Plot data in time
#' 
#' @description
#' This function plots data in time
#' @param data a data frame
#' @param clust a numeric vector of cluster assignemnt.
#' @param orientation a string of orientation of the plot. It must be either "portrait" or "landscape".
#' @param lab a string of label passed on to [ggplot2::xlab()]. Default is NULL.
#' @param palette See `value` of [ggplot2::scale_colour_manual()]. Default is NULL. This uses the palette presets from [MetBrewer::met.brewer()].
#' @return a ggplot.
#' @seealso [ggplot2::labs()], [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @importFrom magrittr %>%
#' @export

plot_time <- function (data, clust, orientation = "portrait", lab = NULL, palette = NULL) {
  clust_levels <- levels(as.factor(clust))
  if (is.null(palette)) {
    clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels)) 
  }
  
  data$clust <- as.factor(clust)
  
  if (is.null(lab)) lab <- "time"
  
  if (orientation == "portrait") {
    p <- ggplot2::ggplot(data, ggplot2::aes(y = clust, x = data[,3], colour = clust, fill = clust)) +
      ggridges::geom_density_ridges(alpha = 0.2, scale = 0.6, show.legend = F) +
      ggplot2::geom_point(shape = "—", size = 3, stroke = 1, alpha = 0.5, show.legend = F) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_fill_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::xlab(lab)
  }
  
  if (orientation == "landscape") {
    p <- ggplot2::ggplot(data, ggplot2::aes(y = clust, x = data[,3], colour = clust, fill = clust)) +
      ggridges::geom_density_ridges(alpha = 0.2, scale = 0.6, show.legend = F) +
      ggplot2::geom_point(shape = "|", size = 3, stroke = 1, alpha = 0.5, show.legend = F) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank()) +
      ggplot2::scale_colour_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::scale_fill_manual(values = setNames(clust_cols, clust_levels)) +
      ggplot2::xlab(lab)
  }
  
  return(p)
}

#' Plot trace summary
#' 
#' @description
#' This function plots 2D objective space.
#' @inheritParams plot_2d_objspace
#' @param colour a string of the attribute to be coloured. It must be "r", "batch" or "k_o".
#' @return a ggplot.
#' @seealso [ggplot2::scale_colour_manual()], [MetBrewer::met.brewer()]
#' @importFrom magrittr %>%
#' @export

plot_trace <- function(pop, alpha = 0.8, colour, palette = NULL) {
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
      p <- p + ggplot2::scale_colour_manual(values = palette,
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
  } else {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = MetBrewer::met.brewer("Hokusai3")) 
    }
    if (colour == "batch") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Egypt", n = length(unique(pop$trace$batch))),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
    if (colour == "k_o") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Archambault", n = length(unique(pop$trace$k_o))),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
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
#' @importFrom magrittr %>%

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
#' @inheritParams plot_2d_objspace
#' @param indicator a string of MOO quality indicator. It must be either "igd", "igd_plus" or "hv".
#' @return a ggplot.
#' @importFrom magrittr %>%
#' @export

plot_moo_quality <- function(pop, indicator) {
  moo_quality <- cbind(batch = as.integer(1:nrow(pop$moo_quality)), pop$moo_quality)
  ggplot2::ggplot(moo_quality) +
    ggplot2::geom_line(ggplot2::aes(x = batch, y = .data[[indicator]]), colour = "#5a97c1") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
}
