# mesoneo
install.packages(c("here", "future", "patchwork", "stringr"))
remotes::install_github("vicyklee/stblob")

here::i_am("paper/scripts/mesoneo.R")
library(here)
library(stblob)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)

plan(multisession)
ggplot2::theme_set(ggplot2::theme_bw())

res <- readRDS( here("paper/output/mesoneo300_res_k5to6.rds"))
res_k5_summary <- res$summary %>% filter(k==5, pareto ==1, pareto_similar == 0) %>%
  arrange(n_removed, space_wcd, desc(time_wcr), desc(time_wce))
plot_idx <- res_k5_summary$idx
# p_s <- plot_space(data, clust = data$Region, coords = c("Longitude", "Latitude"), alpha = 0.5) +
#   labs(colour = "region") +
#   scale_colour_viridis_d(option = 15) 
p_t <- plot_time(data, clust = data$Region, age = "Age_average") +
  theme(legend.position = "none") +
  xlab("time (ybp)") +
  labs(colour = "region") +
  scale_colour_viridis_d(option = 15) +
  scale_fill_viridis_d(option = 15) +
  scale_x_reverse()
# p_mesoneo <- wrap_plots(p_s, p_t, nrow=1, guides = "collect")

p_list <- list()
for (i in seq_along(plot_idx)) {
  j <- plot_idx[i]
  p_s <- plot_space(data, clust = res$clust[j,], coords = c("Longitude", "Latitude"), alpha = 0.5,
                    hull = T, hull_convex_ratio = 0.8) +
    labs(colour = "cluster") + 
    scale_colour_viridis_d(option = 15, na.value = "grey") +
    scale_fill_viridis_d(option = 15)
  p_t <- plot_time(data, clust = res$clust[j,], age = "Age_average") +
    theme(legend.position = "none") +
    xlab("time (ybp)") +
    scale_colour_viridis_d(option = 15, na.value = "grey") +
    scale_fill_viridis_d(option = 15, na.value = "grey") +
    scale_x_reverse() +
    labs(colour = "region") +
    plot_layout(tag_level = 'new')
  p <- wrap_plots(p_s, p_t, nrow = 1, guides = "collect") &
    labs(caption = paste0("ID = ",j,"    r = ",round(res$summary$r[j], digits = 3))) &
    theme(plot.caption = element_text(hjust = 0.5))
  p_list[[i]] <- p
}

p_clust <- wrap_plots(p_list, nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_clust, here("paper/figures/demo1_res.pdf"), width = 12, height = 4)