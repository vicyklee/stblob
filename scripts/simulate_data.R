# Simulate datasets
here::i_am("scripts/simulate_data.R")
library(here)
library(stblob)
#------------------------------------------------------------------------------#
# 3 blobs ####
#------------------------------------------------------------------------------#
set.seed(12)
space <- gen_gaussian_data(n = 3,
                           size = c(100, 100, 100),
                           center = as.matrix(data.frame(x = c(-10,-5,10), y = c(-3,3,0))),
                           sigma = matrix(c(2,0,0,2), nrow = 2))

age <- age <- c(seq(1,300,1)[-seq(1,300,3)][1:100],
                seq(1,300,1)[-seq(1,300,3)][101:200],
                seq(1,300,3))

threeblobs <- cbind(space, age)
threeblobs <- threeblobs[ , c(1,2,4,3)]
usethis::use_data(threeblobs)

#------------------------------------------------------------------------------#
# Others ####
#------------------------------------------------------------------------------#
pathbased <- readr::read_tsv(here("data/pathbased.txt"),col_names = F)
names(pathbased) <- c("x","y","clust")
age <- rep(0,nrow(pathbased))
pathbased <- cbind(pathbased, age)
pathbased <- pathbased[ ,c(1,2,4,3)]
plot_space(pathbased, space = "euclidean",clust = pathbased$clust)
usethis::use_data(pathbased)

jain <- readr::read_tsv(here("data/jain.txt"),col_names = F)
set.seed(123)
age <- as.numeric(sample(1:10000, nrow(jain), replace = T))
names(jain) <- c("x","y","clust")
jain <- cbind(jain, age)
jain <- jain[ ,c(1,2,4,3)]
plot_space(jain, space = "euclidean", clust = jain$clust)
plot_time(jain, jain$clust)
usethis::use_data(jain)

aggregation <- readr::read_tsv(here("data/aggregation.txt"),col_names = F)
set.seed(123)
age <- as.numeric(sample(1:10000, nrow(aggregation), replace = T))
names(aggregation) <- c("x","y","clust")
aggregation <- cbind(aggregation, age)
aggregation <- aggregation[ ,c(1,2,4,3)]
plot_space(aggregation, space = "euclidean", clust = aggregation$clust)
plot_time(aggregation, clust = aggregation$clust)
usethis::use_data(aggregation)

asymmetric <- readr::read_delim(here("data/asymmetric.txt"),col_names = F, delim = " ")
asymmetric_pa <- readr::read_delim(here("data/asymmetric.pa.txt"), col_names = F, delim = " ", skip = 4)
set.seed(123)
age <- as.numeric(sample(1:10000, nrow(asymmetric), replace = T))
asymmetric <- cbind(asymmetric, age, asymmetric_pa)
names(asymmetric) <- c("x","y","age","clust")
plot_space(asymmetric, space = "euclidean",clust = asymmetric$clust)
plot_time(asymmetric, clust = asymmetric$clust)
usethis::use_data(asymmetric)

# factoextra
multishapes <- factoextra::multishapes
noisycircles <- subset(multishapes, shape == 1 | shape == 2)
set.seed(123)
age <- as.numeric(sample(1:10000, nrow(noisycircles), replace = T))
names(noisycircles) <- c("x","y","clust")
noisycircles <- cbind(noisycircles, age)
noisycircles <- noisycircles[ ,c(1,2,4,3)]
plot_space(noisycircles, space = "euclidean", clust = noisycircles$clust)
plot_time(noisycircles, clust = noisycircles$clust)
usethis::use_data(noisycircles)