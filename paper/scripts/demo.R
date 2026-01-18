# Demo
install.packages(c("here", "future", "patchwork"))
remotes::install_github("vicyklee/stblob")

here::i_am("paper/scripts/demo.R")
library(here)
library(stblob)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)

plan(multisession)
ggplot2::theme_set(ggplot2::theme_bw())

#------------------------------------------------------------------------------#
# Simulated dataset 1 ####
#------------------------------------------------------------------------------#
data <- stblob::threeblobs
p_s <- plot_space(data, clust = data$clust, space = "euclidean") + labs(colour = "cluster")
p_t <- plot_time(data, clust = data$clust) + theme(legend.position = "none") + xlab("time")
p_demo1 <- wrap_plots(p_s, p_t, nrow=2, heights = c(1,1), guides = "collect")
save_pdf(p_demo1, here("paper/figures/demo1.pdf"), width = 4, height = 4)

# run the algorithm
set.seed(123)
res <- stblob(data, k = c(2,3), r = c(0.8,1),
              run = 15, batch = 5,
              space_distmethod = "euclidean",
              pareto_similar_ari = 0.8,
              obj = c("space_wcd","-time_wcr","-time_wce"))
saveRDS(res, file = here("paper/output/demo1_res.rds"))

# results
# res <- readRDS(here("paper/output/demo1_res.rds"))
res_summary <- res$summary %>% arrange("space_wcd", desc("time_wcr"), desc("time_wce"))
write.csv(res_summary[,c("idx","batch","k", "run", "r", "space_wcd", "time_wcr", "time_wce", "pareto")],
          here("paper/output/demo1_summary.csv"),row.names = F)

plot_idx <- res$summary %>%
  filter(pareto==1) %>%
  arrange(space_wcd,desc(time_wcr),desc(time_wce)) %>%
  pull(idx)

p_list <- list()
for (i in seq_along(plot_idx)) {
  j <- plot_idx[i]
  p_s <- plot_space(data,clust = res$clust[j,], space = "euclidean", hull = T, hull_convex_ratio = 0.8) +
    theme(legend.position = "none") +
    labs(colour = "cluster")
  if(j %in% res$pareto_idx) {p_s <- p_s + theme(plot.background = element_rect(fill = "#DFF2FF", colour = "#DFF2FF"))} 
  p_t <- plot_time(data,clust = res$clust[j,]) +
    theme(legend.position = "none", plot.caption = element_text(hjust = 0.5)) +
    xlab("time") +
    labs(caption = paste0("ID = ",j,"    r = ",round(res$summary$r[j], digits = 3))) +
    plot_layout(tag_level = 'new')
  if(j %in% res$pareto_idx) {p_t <- p_t + theme(plot.background = element_rect(fill = "#DFF2FF", colour = "#DFF2FF"))} 
  p <- wrap_plots(p_s, p_t, nrow = 2, heights = c(1,1), guides = "collect")
  if(j %in% res$pareto_idx) {p <- p + plot_annotation(theme = theme(plot.background = element_rect(fill = "#DFF2FF")))}
  p_list[[i]] <- p
}

p_clust <- wrap_plots(p_list, nrow = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_clust, here("paper/figures/demo1_clust.pdf"), width = 12, height = 4)

# convergence
p_trace <- plot_trace(res)
p_mooqi <- wrap_plots(plot_mooquality(res,"hv"),
                      plot_mooquality(res,"igd_plus") + plot_layout(tag_level = 'new'),
                      plot_mooquality(res,"igd") + plot_layout(tag_level = 'new'),
                      nrow = 1)
p_conv <- wrap_plots(p_trace, p_mooqi, nrow = 2) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_conv, here("paper/figures/demo1_convergence.pdf"), width = 10, height = 4)

# objective space
p_objs1 <- plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.8) 
p_objs2 <- plot_objspace(res, obj = res$obj[c(1,3)], colour = "pareto", alpha = 0.8) 
p_objs3 <- plot_objspace(res, obj = res$obj[c(2,3)], colour = "pareto", alpha = 0.8) 
p_objs <- wrap_plots(p_objs1, p_objs2, p_objs3, nrow = 1, guides = "collect")
save_pdf(p_objs, here("paper/figures/demo1_objspace.pdf"), width = 9, height = 3)

#------------------------------------------------------------------------------#
# Simulated dataset 1 - WCE^t ####
#------------------------------------------------------------------------------#
data <- stblob::threeblobs
p_s <- plot_space(data, clust = data$clust, space = "euclidean") + labs(colour = "cluster")
p_t <- plot_time(data, clust = data$clust) + theme(legend.position = "none") + xlab("time")
p_demo1 <- wrap_plots(p_s, p_t, nrow=2, heights = c(1,1), guides = "collect")
save_pdf(p_demo1, here("paper/figures/demo1.pdf"), width = 4, height = 4)

# run the algorithm
set.seed(123)
res <- stblob(data, k = 3, r = 0,
              run = 15, batch = 5,
              iter = 20,
              space_distmethod = "euclidean",
              pareto_similar_ari = 0.8,
              obj = c("-time_wcr","-time_wce"),
              filter_intersects = F)
saveRDS(res, file = here("paper/output/demo1_twce_res.rds"))

# results
# res <- readRDS(here("paper/output/demo1_twce_res.rds"))
write.csv(res$summary[res$pareto_idx,c("idx","batch","k", "run", "r", "space_wcd", "time_wcr", "time_wce", "pareto")],
          here("paper/output/demo1_twce_summary.csv"), row.names = F)

p_s <- plot_space(data, clust = res$clust[res$pareto_idx,], space = "euclidean", hull = T, hull_convex_ratio = 0.5) +
  labs(colour = "cluster")
p_t <- plot_time(data, clust = res$clust[res$pareto_idx,]) +
  theme(legend.position = "none", plot.caption = element_text(hjust = 0.5)) +
  xlab("time") +
  labs(caption = paste0("ID = ", res$pareto_idx,"    r = ", round(res$summary$r[res$pareto_idx], digits = 3)))
p_clust <- wrap_plots(p_s, p_t, nrow=2, heights = c(1,1), guides = "collect")

# convergence
p_trace <- plot_trace(res, alpha = 0.2) + theme(legend.position = "none")
p_mooqi <- wrap_plots(plot_mooquality(res,"hv"),
                      plot_mooquality(res,"igd_plus") + plot_layout(tag_level = 'new'),
                      plot_mooquality(res,"igd") + plot_layout(tag_level = 'new'),
                      nrow = 1)
p_conv <- wrap_plots(p_trace, p_mooqi, nrow = 2) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))

# objective space
p_objs <- plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.8) +
  theme(axis.text.x = element_blank())

p_clust[[2]] <- p_clust[[2]] + plot_layout(tag_level = 'new')
p_all <- wrap_plots((p_clust | wrap_plots(p_objs, guides ="collect")) / p_conv) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_all, here("paper/figures/demo1_twce_all.pdf"), width = 8, height = 7.5)

#------------------------------------------------------------------------------#
# Simulated dataset 2 ####
#------------------------------------------------------------------------------#
data <- stblob::threeblobs2
p_s <- plot_space(data, clust = data$clust, space = "euclidean") + labs(colour = "cluster")
p_t <- plot_time(data, clust = data$clust) + theme(legend.position = "none") + xlab("time")
p_demo2 <- wrap_plots(p_s, p_t, nrow=2, heights = c(1,1), guides = "collect")
save_pdf(p_demo2, here("paper/figures/demo2.pdf"), width = 4, height = 4)

# run the algorithm
set.seed(253)
res <- stblob(data, k = c(2,3), r = c(0.8,1),
              run = 15, batch = 5,
              space_distmethod = "euclidean",
              pareto_similar_ari = 0.9,
              obj = c("space_wcd","-time_wcr","-time_wce"))
saveRDS(res,file = here("paper/output/demo2_res.rds"))

# results
# res <- readRDS(here("paper/output/demo2_res.rds"))
res_summary <- res$summary %>% arrange("space_wcd", desc("time_wcr"), desc("time_wce"))
write.csv(res_summary[,c("idx","batch","k", "run", "r", "space_wcd", "time_wcr", "time_wce", "pareto")],
          here("paper/output/demo2_summary.csv"),row.names = F)

p_s <- plot_space(data, clust = res$clust[res$pareto_idx,], space = "euclidean", hull = T, hull_convex_ratio = 0.5) +
  labs(colour = "cluster")
p_t <- plot_time(data, clust = res$clust[res$pareto_idx,]) +
  theme(legend.position = "none", plot.caption = element_text(hjust = 0.5)) +
  xlab("time") +
  labs(caption = paste0("ID = ", res$pareto_idx,"    r = ", round(res$summary$r[res$pareto_idx], digits = 3)))
p_clust <- wrap_plots(p_s, p_t, nrow=2, heights = c(1,1), guides = "collect")

# convergence
p_trace <- plot_trace(res) + theme(legend.position = "none")
p_mooqi <- wrap_plots(plot_mooquality(res,"hv"),
                      plot_mooquality(res,"igd_plus") + plot_layout(tag_level = 'new'),
                      plot_mooquality(res,"igd") + plot_layout(tag_level = 'new'),
                      nrow = 1)
p_conv <- wrap_plots(p_trace, p_mooqi, nrow = 2) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))

# objective space
p_objs1 <- plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.8) +
  theme(axis.text.y = element_blank())
p_objs2 <- plot_objspace(res, obj = res$obj[c(1,3)], colour = "pareto", alpha = 0.8)
p_objs <- wrap_plots(p_objs2, p_objs1, nrow = 2, guides = "collect")

p_clust[[2]] <- p_clust[[2]] + plot_layout(tag_level = 'new')
p_objs[[2]] <- p_objs[[2]] + plot_layout(tag_level = 'new')
p_all <- wrap_plots((p_clust | p_objs)/p_conv) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_all, here("paper/figures/demo2_all.pdf"), width = 8, height = 7.5)


#------------------------------------------------------------------------------#
# Combined simulated datasets figure
#------------------------------------------------------------------------------#
p_demo1[[2]] <- p_demo1[[2]] + plot_layout(tag_level = 'new')
p_demo2[[2]] <- p_demo2[[2]] + plot_layout(tag_level = 'new')
p_demo <- wrap_plots(p_demo1 & theme(legend.position = "none"), p_demo2, guides = "collect") +
  plot_annotation(tag_levels = 'A') & theme(plot.tag=element_text(size = 12, face = "bold"))
save_pdf(p_demo, here("paper/figures/demo.pdf"), width = 7, height = 4)
