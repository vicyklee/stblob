# Demo
here::i_am("scripts/demo")
library(here)
library(stblob)
library(future)
plan(multisession)
library(ggplot2)
ggplot2::theme_set(ggplot2::theme_bw())
library(patchwork)

data <- stblob::threeblobs

plot_space(threeblobs, clust = threeblobs$clust)
plot_space(threeblobs, clust = threeblobs$)

# blob_moo
set.seed(123)
pop_moo <- blob_moo(data, k = c(2,3), r = c(0.5,1), run = 5, batch = 5, space_distmethod = "euclidean")
for (i in pop_moo$pareto_idx) {
  p_s <- plot_space(data,clust = pop_moo$clust[i,])
  p_t <- plot_time(data,clust = pop_moo$clust[i,])
  wrap_plots(p_s, p_t, nrow = 2, guides = "collect")
}

i = 1
