require("dplyr")
require("magrittr") # for %>%
require("tidyr") # for pivot_longer()
require("sf") # for geometry
require("ggplot2") # for ggplots
require("Metbrewer") # for palettes supplied when NULL
require("rnaturalearth") # for plotting map

#------------------------------------------------------------------------------#
# Plotting data ####
#------------------------------------------------------------------------------#
savePDF <- function(p, file, width, height, ...) {
  pdf(file = file, width = width, height = height, ...)
  pdf.options(encoding = 'CP1250')
  print(p)
  invisible(dev.off())
}

plot_2d_fitness <- function (pop, obj, colour = "r", normalise = T, palette = NULL, alpha = 0.8) {
  
  objspace <- pop$summary %>%
    dplyr::select(c(dplyr::all_of(obj), "pareto", "r", "batch", "K_o")) %>%
    dplyr::mutate(pareto = as.factor(pareto),
                  batch = as.factor(batch),
                  K_o = as.factor(K_o))
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
    if (colour == "batch" | colour == "K_o" | colour == "pareto") {
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
    if (colour == "K_o") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Archambault", n = length(unique(objspace$K_o))),
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

plot_space <- function (data, clust = NA, space = "euclidean", hull = F, crs = 3035, lab = NULL) {
  
  clust_levels <- levels(as.factor(clust))
  clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels))
  
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
    # # Add a buffer (in degrees, assuming EPSG:4326)
    if (crs == 4036) {
      buffer <- 5
    }
    if (crs == 3035) {
      buffer <- 500000 # 500 km buffer for 3035
    }
    xlim <- c(bb["xmin"] - buffer, bb["xmax"] + buffer)
    ylim <- c(bb["ymin"] - buffer, bb["ymax"] + buffer)
    
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs = st_crs(3035))
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
    pts <- sf::st_transform(pts, crs = 3035)
    # Get bounding box of your data
    bb <- sf::st_bbox(pts)
    # # Add a buffer (in degrees, assuming EPSG:4326)
    # buffer <- 5
    buffer <- 500000 # 500 km buffer for 3035
    xlim <- c(bb["xmin"] - buffer, bb["xmax"] + buffer)
    ylim <- c(bb["ymin"] - buffer, bb["ymax"] + buffer)
    
    # Hull
    hulls <- pts %>%
      dplyr::filter(!is.na(clust)) %>%
      dplyr::group_by(clust) %>%
      dplyr::summarise(geometry = sf::st_combine(geometry)) %>%
      sf::st_convex_hull()
    
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    world <- sf::st_transform(world, crs = 3035)
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

plot_time <- function (data, clust, orientation = "portrait", lab = NULL) {
  
  clust_levels <- levels(as.factor(clust))
  clust_cols <- MetBrewer::met.brewer("Hokusai3", n = length(clust_levels))
  
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

pivot_trace <- function (df) {
  stat <- c("space_wcss", "time_range_mean", "time_range_sd", "time_evenness_mean", "time_evenness_sd", "ari")
  df <- df %>%
    dplyr::select(c(dplyr::all_of(stat), "iter", "r", "run", "batch", "K_o")) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(stat), names_to = "stat", values_to = "value") %>%
    dplyr::mutate(stat = factor(stat, levels = unique(stat)),
                  run = as.factor(run),
                  batch = as.factor(batch),
                  K_o = as.factor(K_o)) 
  return(df)
}

plot_trace <- function(pop, alpha = 0.8, colour, palette = NULL) {
  trace <- pop$trace %>% pivot_trace()
  p <- ggplot2::ggplot(trace) +
    ggplot2::geom_line(ggplot2::aes(x = iter, y = value, group = interaction(stat, run, batch, r, K_o),
                                    colour = .data[[colour]]), alpha = alpha) +
    ggplot2::facet_wrap(~stat, scales = "free", nrow = 2) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())  +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = max(unique(trace$iter))))
  
  if (!is.null(palette)) {
    if (colour == "r") {
      p <- p + ggplot2::scale_colour_gradientn(colours = palette) 
    }
    if (colour == "batch" | colour == "K_o") {
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
    if (colour == "K_o") {
      p <- p + ggplot2::scale_colour_manual(values = MetBrewer::met.brewer("Archambault", n = length(unique(pop$trace$K_o))),
                                            guide = ggplot2::guide_legend(override.aes = list(alpha = 1)))
    }
  }
  return(p)
}

plot_moo_quality <- function(pop, indicator = "hv", colour = "#5a97c1") {
  moo_quality <- cbind(batch = as.integer(1:nrow(pop$moo_quality)), pop$moo_quality)
  ggplot2::ggplot(moo_quality) +
    ggplot2::geom_line(ggplot2::aes(x = batch, y = .data[[indicator]]), colour = colour) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks())
}

