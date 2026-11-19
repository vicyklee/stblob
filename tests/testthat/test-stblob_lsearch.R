test_that("stblob_lsearch works", {
  # simplest working case
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0.8)
  expect_s3_class(sol,  class = "stblob_sol")
  
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0.8,
                        optim_type_diversity = TRUE)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # specifying w_time
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 2,
                        w_space = 0.6,
                        optim_type_diversity = TRUE,
                        w_time = 0.2)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # specifying w_type
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 2,
                        w_space = 0.6,
                        optim_type_diversity = TRUE,
                        w_type = 0.2)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # specifying all weights
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 2,
                        w_space = 0.6,
                        optim_type_diversity = TRUE,
                        w_time = 0.3,
                        w_type = 0.1)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # specifying w_space > 1
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 2,
                              w_space = 2,
                              optim_type_diversity = TRUE),
               "'w_space' must be between 0 and 1")
  
  # specifying sum of weights > 1
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 2,
                              w_space = 0.5,
                              optim_type_diversity = TRUE,
                              w_time = 0.6),
               "'w_space', 'w_time' and 'w_type' must sum up to 1")
  
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 2,
                              w_space = 0.5,
                              optim_type_diversity = TRUE,
                              w_time = 0.6,
                              w_type = 0.2),
               "'w_space', 'w_time' and 'w_type' must sum up to 1")
  
  # specifying weight not used
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 2,
                        w_space = 0.6,
                        optim_type_diversity = FALSE,
                        w_type = 0.1)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # specifying k < 2
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 1,
                              w_space = 0.5),
               "'k' must be at least 2 and smaller than the number of samples")
  
  # specifying k >= N
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 30,
                              w_space = 0.5),
               "'k' must be at least 2 and smaller than the number of samples")
  
  # turn on outlier removal
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 1,
                        outlier_removal = TRUE,
                        outlier_iqrm = 1.5)
  expect_s3_class(sol,  class = "stblob_sol")
  expect_equal(sol$summary$n_outliers, 1)
  
  # turn off sf_use_s2
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0.8,
                        sf_use_s2 = FALSE)
  expect_s3_class(sol,  class = "stblob_sol")
  
  # filter_intersects
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0)
  expect_equal(sol$status,  2)
  
  set.seed(123)
  expect_message(stblob_lsearch(data = threeblobs_s,
                                k = 3,
                                w_space = 0),
                 "Status code 2: Intersecting clusters")
  
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0,
                        filter_intersects = FALSE)
  expect_equal(sol$status,  0)
  
  # filter_clustsize
  set.seed(55)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 7,
                        w_space = 1)
  expect_equal(sol$status, 3)
  
  set.seed(55)
  expect_message(stblob_lsearch(data = threeblobs_s,
                                k = 7,
                                w_space = 1),
                 "Status code 3: Cluster size below threshold")
  
  set.seed(55)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 7,
                        w_space = 1,
                        filter_clustsize = FALSE)
  expect_equal(sol$status, 0)
  
  # space_distmat mode
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0.8,
                        space_distmat = threeblobs_s_sdistm)
  expect_s3_class(sol,  class = "stblob_sol")
  
  set.seed(123)
  sol <- stblob_lsearch(data = threeblobs_s,
                        k = 3,
                        w_space = 0.8,
                        space_distmat = as.dist(threeblobs_s_sdistm))
  expect_s3_class(sol,  class = "stblob_sol")
  
  # space_dist wrong dim
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 3,
                              w_space = 0.8,
                              space_distmat = threeblobs_s_sdistm[1:20,1:20]),
               "The dimension of 'space_distmat' does not correspond to the number of rows of 'data'")
  
  # specifying the wrong columns
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 3,
                              w_space = 0.8,
                              coords = c("g","y")),
               "'coords' does not match any column names")
  
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 3,
                              w_space = 0.8,
                              age = "g"),
               "'age' does not match any column names")
  
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 3,
                              w_space = 0.8,
                              optim_type_diversity = TRUE,
                              type = "g"),
               "'type' does not match any column names")
  
  expect_error(stblob_lsearch(data = threeblobs_s,
                              k = 3,
                              w_space = 0.8,
                              coords = c("x","x")),
               "Same 'data' column used twice")
})
