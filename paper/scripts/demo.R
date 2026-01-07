# Demo
here::i_am("scripts/demo.R")
library(here)
remotes::install_github("vicyklee/stblob", force = T)
# source(here("R/functions_localsearch.R"))
# source(here("R/functions_moo.R"))
# source(here("R/functions_plotting.R"))
library(stblob)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)

plan(multisession)
ggplot2::theme_set(ggplot2::theme_bw())
#------------------------------------------------------------------------------#
# Simulated dataset 1
#------------------------------------------------------------------------------#
data <- stblob::threeblobs
p_s <- plot_space(data, clust = threeblobs$clust, space = "euclidean") + labs(colour = "Cluster")
p_t <- plot_time(data, clust = threeblobs$clust) + theme(legend.position = "none") + xlab("Time")
p_demo1 <- wrap_plots(p_s,p_t,nrow=2, heights = c(1,1), guides = "collect")
save_pdf(p_demo1, here("figures/demo.pdf"), width = 4, height = 4)

# run the algorithm
set.seed(123)
res <- blob_moo(data, k = c(2,3), r = c(0.8,1),
                run = 15, batch = 5,
                space_distmethod = "euclidean",
                pareto_similar_ari = 0.9,
                obj = c("space_wcd","-time_wcr","-time_wce"))
saveRDS(res,file = here("output/demo_res.rds"))

# results
for (i in res$summary$idx) {
  p_s <- plot_space(data,clust = res$clust[i,], space = "euclidean", hull = T, hull_convex_ratio = 1)
  p_t <- plot_time(data,clust = res$clust[i,]) + theme(legend.position = "none") + xlab("Time")
  print(wrap_plots(p_s, p_t, nrow = 2, guides = "collect"))
}

# convergence
plot_trace(res)
plot_mooquality(res,"hv")
plot_mooquality(res,"igd_plus")
plot_mooquality(res,"igd")

# objective space
plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.8) 
plot_objspace(res, obj = res$obj[c(1,3)], colour = "pareto", alpha = 0.8) 

#------------------------------------------------------------------------------#
# Simulated dataset 2
#------------------------------------------------------------------------------#
data <- stblob::threeblobs2
p_s <- plot_space(data, clust = threeblobs$clust, space = "euclidean") + labs(colour = "Cluster")
p_t <- plot_time(data, clust = threeblobs$clust) + theme(legend.position = "none") + xlab("Time")
p_demo2 <- wrap_plots(p_s,p_t,nrow=2, heights = c(1,1), guides = "collect")
save_pdf(p_demo2, here("figures/demo2.pdf"), width = 4, height = 4)

# run the algorithm
set.seed(253)
res <- blob_moo(data, k = c(2,3), r = c(0.8,1),
                run = 15, batch = 5,
                space_distmethod = "euclidean",
                pareto_similar_ari = 0.9,
                obj = c("space_wcd","-time_wcr","-time_wce"))
saveRDS(res,file = here("output/demo_res.rds"))

# results
for (i in res$summary$idx) {
  p_s <- plot_space(data,clust = res$clust[i,], space = "euclidean", hull = T, hull_convex_ratio = 1)
  p_t <- plot_time(data,clust = res$clust[i,]) + theme(legend.position = "none") + xlab("Time")
  print(wrap_plots(p_s, p_t, nrow = 2, guides = "collect"))
}

# convergence
plot_trace(res)
plot_mooquality(res,"hv")
plot_mooquality(res,"igd_plus")
plot_mooquality(res,"igd")

# objective space
plot_objspace(res, obj = res$obj[c(1,2)], colour = "pareto", alpha = 0.8) 
plot_objspace(res, obj = res$obj[c(1,3)], colour = "pareto", alpha = 0.8) 

#------------------------------------------------------------------------------#
# Archaeological (aDNA) dataset 
#------------------------------------------------------------------------------#
data <- readr::read_tsv(here("data_paper/mesoneo.tsv")) %>%
  filter(stringr::str_detect(Region, "Europe") | Region == "WesternAsia", 
         Age_average <= 15000)

set.seed(123)
data <- data[sample(1:nrow(data), 300),]

# blob_moo
set.seed(123)
res <- blob_moo(data, k = c(5,6), r = c(0.9,1), run = 15, batch = 10,
                coords = c("Longitude","Latitude"), age = "Age_average",
                pareto_similar_ari = 0.8, obj = c("space_wcd", "-time_wcr", "-time_wce", "n_removed"))
saveRDS(res,file = here("output/mesoneo300_res_k5to6.rds"))

set.seed(123)
res <- blob_moo(data, k = c(7,10), r = c(0.9,1), run = 20, batch = 20,
                coords = c("Longitude","Latitude"), age = "Age_average",
                pareto_similar_ari = 0.8)

eval_moo()
saveRDS(res,file = here("output/mesoneo300_res_k7to10.rds"))

for (i in res$pareto_idx) {
  p_s <- plot_space(data,clust = res$clust[i,], space = "earth",
                    coords = c("Longitude","Latitude"), hull = T, hull_convex_ratio = 1) +
    scale_colour_viridis_d(option = 15, na.value = "red") +
    scale_fill_viridis_d(option = 15)
  p_t <- plot_time(data,clust = res$clust[i,], age = "Age_average") +
    scale_colour_viridis_d(option = 15, na.value = "red") +
    scale_fill_viridis_d(option = 15)
  print(wrap_plots(p_s, p_t, nrow = 2, guides = "collect", heights = c(1,1)))
}


plot_mooquality(res, "hv")
plot_objspace(res, res$obj[c(1,2)], colour = "pareto")
plot_objspace(res, res$obj[c(1,3)], colour = "pareto")
plot_trace(res)


