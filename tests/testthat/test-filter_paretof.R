test_that("filter_paretof works", {
  set.seed(123)
  moo <- stblob(threeblobs_s,
                k = 2:3,
                w_space  = c(0.5,1),
                optim_type_diversity = TRUE)
  expect_s3_class(moo, "stblob_moo")
  
  idx_pf <- filter_paretof(moo)
  expect_vector(idx_pf)
  expect_type(idx_pf, "integer")
  
  # expect no effect
  idx_pf <- filter_paretof(moo, 1, 0)
  expect_equal(idx_pf, moo$summary$idx_pf[!is.na(moo$summary$idx_pf)])
  
  # expect NULL
  moo$clust_pf <- NULL
  expect_message(idx_pf <- filter_paretof(moo), "No solution to filter")
  expect_null(idx_pf)
  
  # expect 1 solution
  moo$clust_pf <- c(1,1,2,3,2,1,1)
  expect_message(filter_paretof(moo), "Only one soluton on the Pareto front: returning 'idx_pf' of 1")
})
