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
# Archaeological (aDNA) dataset N=1207 k=c(5,12) ####
#------------------------------------------------------------------------------#
data <- readr::read_tsv(here("paper/data/mesoneo.tsv")) %>%
  filter(stringr::str_detect(Region, "Europe") | Region == "WesternAsia",
         Country != "Faroes" & Country != "Greenland" & Country != "Iceland",
         Age_average <= 15000)

start <- Sys.time()
set.seed(123)
res <- stblob(data, k = c(5,12), r = c(0.95,1), run = 50, batch = 100,
              coords = c("Longitude","Latitude"), age = "Age_average",
              pareto_similar_ari = 0.8, obj = c("space_wcd", "-time_wcr", "-time_wce"))
end <- Sys.time()
print(end-start)
saveRDS(res,file = here("paper/output/mesoneo_res_k5to12.rds"))
