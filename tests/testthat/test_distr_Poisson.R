context("Distributions - Poisson")

test_that("basic tests",{
  poiss.obj <- ExpFam_dist_Poisson()
  expect_is(poiss.obj, "ExpFam_dist")

  poiss.spec.obj <- dist_Poisson()
  expect_is(poiss.spec.obj, "ExpFam_dist")
  expect_is(poiss.spec.obj, "N_ExpFam_dist")

  test_basic_invariance.ExpFam_dist(
    poiss.spec.obj,
    valid.eta = seq(from = -10, to = 5, by = 0.1),
    valid.theta = seq(from = 0.1, to = 10, by = 0.1),
    valid.values = seq(from = 0, to = 20, by = 1)
    )

  expect_equal(poiss.spec.obj$N.density.theta(x = 0:10, theta = 0.6), dpois(0:10, 0.6))   #density for lambda = 0.6
  expect_equal(poiss.spec.obj$N.density.eta(x = 0:10, eta = log(0.6)), dpois(0:10, 0.6))   #density for lambda = 0.6


  expect_equal(ExpFam_density_theta(poiss.obj, 0.6)(1), dpois(1, 0.6))   #density for lambda = 0.6
  expect_equal(ExpFam_density_theta(poiss.obj, 0.6)(1:3), dpois(1:3, 0.6))   #density for lambda = 0.6
  expect_equal(ExpFam_density_theta(poiss.obj, 1.6)(1), dpois(1, 1.6))   #density for lambda = 1.6
  expect_equal(ExpFam_density_theta(poiss.obj, 1.6)(1:3), dpois(1:3, 1.6))   #density for lambda = 1.6

  test.eta <- log(0.6)
  expect_equal(ExpFam_density_eta(poiss.obj, test.eta)(1), dpois(1, exp(test.eta)))   #density for lambda = 0.6
  expect_equal(ExpFam_density_eta(poiss.obj, test.eta)(1:3), dpois(1:3, exp(test.eta)))   #density for lambda = 0.6
  test.eta <- log(1.6)
  expect_equal(ExpFam_density_eta(poiss.obj, test.eta)(1), dpois(1, 1.6))   #density for lambda = 1.6
  expect_equal(ExpFam_density_eta(poiss.obj, test.eta)(1:3), dpois(1:3, 1.6))   #density for lambda = 1.6

  expect_equal(ExpFam_density_eta(poiss.obj, 0)(1),  dpois(1, exp(0)))   #density for value 1

  # theta.test <- 0.6
  # d.t.06 <- ExpFam_density_theta(bern.obj, theta.test)
  # expect_true(all(c(d.t.06(0) == 1 - theta.test, d.t.06(1) == theta.test)))
  #
  # d.t.07 <- ExpFam_density_theta(bern.obj, 0.7)
  # expect_true(all(c(d.t.07(0) == 1 - 0.7, d.t.07(1) == 0.7)))
  #
  # eta.test <- 0.1
  # eta.test.theta <- bern.obj$theta.from.eta(eta.test)
  # d.e.0 <- ExpFam_density_eta(bern.obj, eta.test)
  # expect_true(sum((c(d.e.0(0) - (1 - eta.test.theta), d.e.0(1) - eta.test.theta))^2) < 1e-16)

})

