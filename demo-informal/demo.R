devtools::install()
library(stblob)
# show progress bar
progressr::handlers(global = TRUE)

# Demo dataset ----------------------------------------------------------------
threeblobs_subset
# add an outlier
# threeblobs
# threeblobs2 <- rbind(threeblobs_subset, c(301, 1, -5, 100, 3))
threeblobs2 <- threeblobs_subset
# threeblobs2$type <- rep(c(1,2,3), length.out = 301)

# Visualise dataset -----------------------------------------------------------
plot(threeblobs2[c("x","y")], col = heat.colors(300)[threeblobs2$z],
     main = "3 spatial blobs\ncoloured by ages")
plot(threeblobs2[c("x","y")], col = topo.colors(3)[threeblobs2$type],
     main = "3 spatial blobs\ncoloured by data types")


# stblob ----------------------------------------------------------------------
moo <- stblob(data = threeblobs2,
              k = 2:3,
              w_space = c(0.8,1),
              optim_type_diversity = TRUE,
              run = 10,
              batch = 10,
              space_distmethod = "euclidean")

# filter_paretof --------------------------------------------------------------
moo$paretof_idx
paretof_idx_subset <- filter_paretof(moo, paretof_ari = 0.8, z_dist = 0.1)

plot(threeblobs2[c("x","y")], col = moo$clust[moo$paretof_idx[2],],
     main = "3 spatial blobs\ncoloured by clusters")







# stblob_lsearch --------------------------------------------------------------
set.seed(145)
sol <- stblob_lsearch(data = threeblobs2,
                      k = 3,
                      w_space = 0,
                      optim_type_diversity = TRUE,
                      w_type = 1,
                      coords = c("x","y"),
                      age = "z",
                      type = "type",
                      space_distmethod = "euclidean",
                      filter_intersects = FALSE,
                      filter_clustsize = FALSE,
                      outlier_removal = TRUE)

# colour the NA with outlier_removal = TRUE
sol$clust[is.na(sol$clust)] <- 4
plot(sol$data[c("x","y")], col = sol$clust,
     main = "3 spatial blobs\ncoloured by clusters")

# stblob_pop ------------------------------------------------------------------
set.seed(1314)
pop <- stblob_populate(data = threeblobs2,
                       k = 2:4,
                       w_space = c(0.5,1),
                       optim_type_diversity = TRUE,
                       run = 5,
                       batch = 5,
                       space_distmethod = "euclidean")

# combine_pop -----------------------------------------------------------------
set.seed(234)
pop2 <- stblob_populate(data = threeblobs2,
                        k = 2:5,
                        w_space = c(0.8,1),
                        optim_type_diversity = TRUE,
                        run = 5,
                        batch = 3,
                        space_distmethod = "euclidean")

set.seed(123)
pop3 <- stblob_populate(data = threeblobs2,
                        k = 3:5,
                        w_space = c(0.7,1),
                        optim_type_diversity = TRUE,
                        run = 5,
                        batch = 3,
                        space_distmethod = "euclidean")

pop_combined <- combine_pop(pop,pop2,pop3)

# eval_moo --------------------------------------------------------------------
moo2 <- eval_moo(pop_combined)

