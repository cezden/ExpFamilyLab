context("Distributions - Exponential")

test_that("basic tests",{
  exp.spec.obj <- dist_Exponential()
  expect_is(exp.spec.obj, "N_ExpFam_dist")

  test_basic_invariance.ExpFam_dist(
    exp.spec.obj,
    valid.eta = seq(from = -10, to = -0.1, by = 0.1),
    valid.theta = seq(from = 0.1, to = 10, by = 0.1),
    valid.values = seq(from = 0, to = 20, by = 0.1)
  )

  expect_equal(exp.spec.obj$N.density.theta(x = 0:10, theta = 0.6), dexp(0:10, rate = 0.6))   #density for rate = 0.6

})

