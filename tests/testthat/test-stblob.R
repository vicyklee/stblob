test_that("stblob works", {
  # simplest working case
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = c(0.5,1))
  expect_s3_class(moo, "stblob_moo")
  
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = c(0.5,1),
                optim_type_diversity = TRUE)
  expect_s3_class(moo, "stblob_moo")
  
  # fixed weights
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = 0.7,
                optim_type_diversity = TRUE,
                w_time = 0.2)
  expect_s3_class(moo, "stblob_moo")
  
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = 0.7,
                optim_type_diversity = TRUE,
                w_type = 0.2)
  expect_s3_class(moo, "stblob_moo")
  
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = 0.7,
                optim_type_diversity = TRUE,
                w_time = 0.1,
                w_type = 0.2)
  expect_s3_class(moo, "stblob_moo")
  
  # specifying the fixed weights such that sum > 1
  expect_error(stblob(threeblobs_s,
                      k = 2:3,
                      w_space  = 0.6,
                      w_time = 0.5),
               "'w_space' and 'w_time' must sum up to 1")
  
  expect_error(stblob(threeblobs_s,
                      k = 2:3,
                      w_space  = 0.6,
                      optim_type_diversity = TRUE,
                      w_type = 0.5),
               "'w_space', 'w_time' and 'w_type' must sum up to 1")
  
  # specifying a range for w_time 
  expect_error(stblob(threeblobs_s,
                      k = 2:3,
                      w_space  = c(0.5,1),
                      w_time = c(0.2,1)),
               "'w_time' must be a numeric value of length 1")
  
  # specifying a range for w_type 
  expect_error(stblob(threeblobs_s,
                      k = 2:3,
                      w_space  = c(0.5,1),
                      optim_type_diversity = TRUE,
                      w_type = c(0.2,1)),
               "'w_type' must be a numeric value of length 1")
  
  # specifying w_space > 1
  expect_error(stblob(threeblobs_s,
                      k = 2:3,
                      w_space  = c(0.5,2)),
               "'w_space' must be between 0 and 1")
  
  # specifying k < 2
  expect_error(stblob(data = threeblobs_s,
                      k = 1:2,
                      w_space = c(0.5,1)),
               "'k' must be at least 2 and smaller than the number of samples")
  
  # specifying k >= N
  expect_error(stblob(data = threeblobs_s,
                      k = 2:30,
                      w_space = c(0.5,1)),
               "'k' must be at least 2 and smaller than the number of samples")
  
  # turn on outlier removal
  set.seed(123)
  moo <- stblob(data = threeblobs_s,
                k = 2:3,
                w_space = c(0.5,1),
                outlier_removal = TRUE,
                outlier_iqrm = 1.5)
  expect_s3_class(moo,  class = "stblob_moo")
  expect_true(any(moo$summary$n_outliers > 1))
  
  # turn off sf_use_s2
  set.seed(123)
  moo <- stblob(data = threeblobs_s,
                k = 2:3,
                w_space = c(0.5,1),
                sf_use_s2 = FALSE)
  expect_s3_class(moo,  class = "stblob_moo")
  
  # filter_intersects
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
  moo <- stblob(data = threeblobs_s,
                k = 2:3,
                w_space = c(0.5,1),
                space_distmat = threeblobs_s_sdistm)
  expect_s3_class(moo,  class = "stblob_moo")
  
  set.seed(123)
  moo <- stblob(data = threeblobs_s,
                k = 2:3,
                w_space = c(0.5,1),
                space_distmat = as.dist(threeblobs_s_sdistm))
  expect_s3_class(moo,  class = "stblob_moo")
  
  # space_dist wrong dim
  expect_error(stblob(data = threeblobs_s,
                      k = 2:3,
                      w_space = c(0.5,1),
                      space_distmat = threeblobs_s_sdistm[1:20,1:20]),
               "The dimension of 'space_distmat' does not correspond to the number of rows of 'data'")
  
  # specifying the wrong columns
  expect_error(stblob(data = threeblobs_s,
                      k = 2:3,
                      w_space = c(0.5,1),
                      coords = c("g","y")),
               "'coords' does not match any column names")
  
  expect_error(stblob(data = threeblobs_s,
                      k = 2:3,
                      w_space = c(0.5,1),
                      age = "g"),
               "'age' does not match any column names")
  
  expect_error(stblob(data = threeblobs_s,
                      k = 2:3,
                      w_space = c(0.5,1),
                      optim_type_diversity = TRUE,
                      type = "g"),
               "'type' does not match any column names")
  
  expect_error(stblob(data = threeblobs_s,
                      k = 2:3,
                      w_space = c(0.5,1),
                      coords = c("x","x")),
               "Same 'data' column used twice")
})
