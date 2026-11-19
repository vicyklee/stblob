test_that("compute_distmat works", {
  distmat <- compute_distmat(threeblobs_s[c("x","y")])
  expect_type(distmat, "double")
  expect_all_equal(dim(distmat), nrow(threeblobs_s))
})
