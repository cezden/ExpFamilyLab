context("Distributions - Bernoulli")

test_that("basic tests",{
  bern.obj <- ExpFam_dist_Bernoulli()
  expect_is(bern.obj, "ExpFam_dist")

  bern.obj <- dist_Bernoulli()

  test_basic_invariance(
    bern.obj,
    valid.eta = seq(from = -10, to = 10, by = 0.1),
    valid.theta = seq(from = 0.1, to = 0.9, by = 0.1),
    valid.values = c(0, 1)
  )



  expect_equal(ExpFam_density_theta(bern.obj, 0.6)(1), 0.6)   #density for value 1

  expect_equal(ExpFam_density_eta(bern.obj, 0)(1),  0.5)   #density for value 1

  theta.test <- 0.6
  d.t.06 <- ExpFam_density_theta(bern.obj, theta.test)
  expect_true(all(c(d.t.06(0) == 1 - theta.test, d.t.06(1) == theta.test)))

  d.t.07 <- ExpFam_density_theta(bern.obj, 0.7)
  expect_true(all(c(d.t.07(0) == 1 - 0.7, d.t.07(1) == 0.7)))

  eta.test <- 0.1
  eta.test.theta <- bern.obj$theta.from.eta(eta.test)
  d.e.0 <- ExpFam_density_eta(bern.obj, eta.test)
  expect_true(sum((c(d.e.0(0) - (1 - eta.test.theta), d.e.0(1) - eta.test.theta))^2) < 1e-16)

})

