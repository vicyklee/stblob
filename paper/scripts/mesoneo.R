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

#------------------------------------------------------------------------------#
# Archaeological (aDNA) dataset 
#------------------------------------------------------------------------------#
data <- readr::read_tsv(here("paper/data/mesoneo.tsv")) %>%
  filter(stringr::str_detect(Region, "Europe") | Region == "WesternAsia",
         Country != "Faroes" & Country != "Greenland" & Country != "Iceland",
         Age_average <= 15000)

set.seed(123)
data <- data[sample(1:nrow(data), 300),]

# p_s <- plot_space(data, clust = data$Region, coords = c("Longitude", "Latitude"), alpha = 0.5) +
#   labs(colour = "region") +
#   scale_colour_viridis_d(option = 15) 
# p_t <- plot_time(data, clust = data$Region, age = "Age_average") +
#   theme(legend.position = "none") +
#   xlab("time (ybp)") +
#   labs(colour = "region") +
#   scale_colour_viridis_d(option = 15) +
#   scale_fill_viridis_d(option = 15) +
#   scale_x_reverse()

# p_mesoneo <- wrap_plots(p_s, p_t, nrow=1, guides = "collect")
# save_pdf(p_mesoneo, here("paper/figures/mesoneo300.pdf"), width = 10, height = 4)

# blob_moo
start <- Sys.time()
set.seed(123)
res <- blob_moo(data, k = c(5,6), r = c(0.9,1), run = 20, batch = 30,
                coords = c("Longitude","Latitude"), age = "Age_average",
                pareto_similar_ari = 0.8, obj = c("n_removed", "space_wcd", "-time_wcr", "-time_wce"))
end <- Sys.time()
print(end-start)
saveRDS(res,file = here("paper/output/mesoneo300_res_k5to6.rds"))
# res <- readRDS( here("paper/output/mesoneo300_res_k5to6.rds"))