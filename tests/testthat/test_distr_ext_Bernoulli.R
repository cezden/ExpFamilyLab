context("Distributions (ext) - Bernoulli")

test_that("basic tests",{
  bern.obj <- dist_ext_Bernoulli()
  expect_is(bern.obj, "ExpFam_dist_ext")

  bern.canonical <- bern.obj[["canonical"]]
  bern.mean <- bern.obj[["mean"]]

  dens1 <- ExpFam_density(bern.canonical, eta = 0.7, num.opt = TRUE)
  dens2 <- ExpFam_density(bern.canonical, eta = 0.7, num.opt = FALSE)

  expect_equal(dens1(c(0,1)), dens2(c(0,1)))

  bern.par.mean <- ExpFam_bind_parametrization(bern.obj, "mean")
  test_basic_invariance.ExpFam_dist(
    bern.par.mean,
    valid.eta = seq(from = -10, to = 10, by = 0.1),
    valid.values = c(0, 1),
    valid.theta = seq(from = 0.1, to = 0.9, by = 0.1)
  )


  test_basic_invariance.ExpFam_dist_ext(
    bern.obj,
    valid.eta = seq(from = -10, to = 10, by = 0.1),
    valid.values = c(0, 1),
    list(
      "mean" = list(valid.theta = seq(from = 0.1, to = 0.9, by = 0.1))
    )
  )



  expect_equal(ExpFam_density_theta(bern.par.mean, 0.6)(1), 0.6)   #density for value 1

  expect_equal(ExpFam_density_eta(bern.par.mean, 0)(1),  0.5)   #density for value 1

  theta.test <- 0.6
  d.t.06 <- ExpFam_density_theta(bern.par.mean, theta.test)
  expect_true(all(c(d.t.06(0) == 1 - theta.test, d.t.06(1) == theta.test)))

  d.t.07 <- ExpFam_density_theta(bern.par.mean, 0.7)
  expect_true(all(c(d.t.07(0) == 1 - 0.7, d.t.07(1) == 0.7)))

  eta.test <- 0.1
  eta.test.theta <- bern.par.mean$theta.from.eta(eta.test)
  d.e.0 <- ExpFam_density_eta(bern.par.mean, eta.test)
  expect_true(sum((c(d.e.0(0) - (1 - eta.test.theta), d.e.0(1) - eta.test.theta))^2) < 1e-16)

})

