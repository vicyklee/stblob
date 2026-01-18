# mesoneo
install.packages(c("here", "future", "patchwork", "stringr"))
remotes::install_github("vicyklee/stblob")

here::i_am("paper/scripts/mesoneo_plots.R")
library(here)
library(stblob)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)

plan(multisession)
ggplot2::theme_set(ggplot2::theme_bw())


#------------------------------------------------------------------------------#
# Archaeological (aDNA) dataset full 
#------------------------------------------------------------------------------#
p_s <- plot_space(data, clust = data$Region, coords = c("Longitude", "Latitude"), alpha = 0.5) +
  labs(colour = "region")
p_t <- plot_time(data, clust = data$Region, age = "Age_average") +
  theme(legend.position = "none") +
  xlab("time (ybp)") +
  labs(colour = "region") +
  scale_x_reverse()

p_mesoneo <- wrap_plots(p_s, p_t, nrow=1, guides = "collect")
save_pdf(p_mesoneo, here("paper/figures/mesoneo.pdf"), width = 8, height = 3)

# Results
res <- readRDS(here("paper/output/mesoneo_res_k5to12_maxna0.rds"))
# Time difference of 1.179481 hours, 96 cores

data <- readr::read_tsv(here("paper/data/mesoneo.tsv")) %>%
  filter(stringr::str_detect(Region, "Europe") | Region == "WesternAsia",
         Country != "Faroes" & Country != "Greenland" & Country != "Iceland",
         Age_average <= 15000)

# PCP plot
p_pcp <- plot_pcp(res, colour = "k", linewidth = 0.5) + xlab("objective") + ylab("objective value")
save_pdf(p_pcp, here("paper/figures/mesoneo_k5to12_pcp.pdf"), width = 6, height = 3)

# plot the first 6 ranks
res_summary <- res$summary %>% filter(pareto==1 & pareto_similar==0) %>%
  arrange(space_wcd, desc(time_wcr), desc(time_wce))

plot_idx <- res_summary$idx[1:6]

p_list <- list()
for (i in seq_along(plot_idx)) {
  j <- plot_idx[i]
  
  p_s <- plot_space(data, clust = res$clust[j,], coords = c("Longitude", "Latitude"),
                    size = 0.5, alpha = 0.5, hull = T, hull_convex_ratio = 0.8) +
    guides(colour=guide_legend(ncol=2)) +
    labs(colour = "cluster",
         caption = paste0("ID = ",j,"    r = ",round(res$summary$r[j], digits = 3))) &
    theme(plot.caption = element_text(hjust = 0.5))
  
  p_t <- plot_time(data, clust = res$clust[j,], age = "Age_average") +
    theme(legend.position = "none") +
    xlab("time (ybp)") +
    scale_x_reverse() +
    labs(colour = "cluster") +
    plot_layout(tag_level = 'new')
  
  p <- wrap_plots(p_s, p_t, nrow = 1, guides = "collect")
  
  p_list[[i]] <- p
}

p_clust <- wrap_plots(p_list, nrow = 3) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_clust, here("paper/figures/mesoneo_k5to12_rank1to6_clust.pdf"), width = 15, height = 9)

plot_idx <- res_summary$idx[7:14]
p_list <- list()
for (i in seq_along(plot_idx)) {
  j <- plot_idx[i]
  
  p_s <- plot_space(data, clust = res$clust[j,], coords = c("Longitude", "Latitude"),
                    size = 0.5, alpha = 0.5, hull = T, hull_convex_ratio = 0.8) +
    guides(colour=guide_legend(ncol=2)) +
    labs(colour = "cluster",
         caption = paste0("ID = ",j,"    r = ",round(res$summary$r[j], digits = 3))) &
    theme(plot.caption = element_text(hjust = 0.5))
  
  p_t <- plot_time(data, clust = res$clust[j,], age = "Age_average") +
    theme(legend.position = "none") +
    xlab("time (ybp)") +
    scale_x_reverse() +
    labs(colour = "cluster") +
    plot_layout(tag_level = 'new')
  
  p <- wrap_plots(p_s, p_t, nrow = 1, guides = "collect")
  
  p_list[[i]] <- p
}

p_clust <- wrap_plots(p_list, nrow = 4) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_clust, here("paper/figures/mesoneo_k5to12_rank7to14_clust.pdf"), width = 15, height = 12)

# convergence
p_trace <- plot_trace(res, alpha = 0.2)
p_mooqi <- wrap_plots(plot_mooquality(res,"hv"),
                      plot_mooquality(res,"igd_plus") + plot_layout(tag_level = 'new'),
                      plot_mooquality(res,"igd") + plot_layout(tag_level = 'new'),
                      nrow = 1)
p_conv <- wrap_plots(p_trace, p_mooqi, nrow = 2) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_conv, here("paper/figures/mesoneo_k5to12_convergence.pdf"), width = 10, height = 4)

# objective space
p_objs1 <- plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.5)
p_objs2 <- plot_objspace(res, obj = res$obj[c(1,3)], colour = "pareto", alpha = 0.5) + plot_layout(tag_level = 'new')
p_objs3 <- plot_objspace(res, obj = res$obj[c(2,3)], colour = "pareto", alpha = 0.5) + plot_layout(tag_level = 'new') 
p_objs <- wrap_plots(p_objs1, p_objs2, p_objs3, nrow = 1, guides = "collect")
save_pdf(p_objs, here("paper/figures/mesoneo_k5to12_objspace.pdf"), width = 9, height = 3)

# all except assignments
p_all <-  wrap_plots(p_objs, p_conv, nrow = 2, heights = c(1,1.5)) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_all, here("paper/figures/mesoneo_k5to12_all.pdf"), width = 10, height = 7.5)

write.csv(res_summary[,c("idx","batch","k", "run", "r", "space_wcd", "time_wcr", "time_wce", "pareto")],
          here("paper/output/mesoneo_k5to12_summary.csv"), row.names = F)
