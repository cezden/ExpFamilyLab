context("Distributions (ext) - Exponential")

test_that("basic tests",{
  distr.obj <- dist_ext_Exponential()

  test_basic_invariance.ExpFam_dist_ext(
    distr.obj,
    valid.eta = seq(from = -10, to = -0.1, by = 0.1),
    valid.values = seq(from = 0.01, to = 10, by = 0.5),
    list(
      "mean" = list(valid.theta = seq(from = 0.1, to = 2, by = 0.1)),
      "rate" = list(valid.theta = seq(from = 0.1, to = 2, by = 0.1))
    )
  )

  distr.par.mean <- bind_parametrization(distr.obj, "mean")

  test.mean <- 0.7
  dens1 <- ExpFam_density_theta(distr.par.mean, theta = test.mean, num.opt = TRUE)
  dens2 <- ExpFam_density_theta(distr.par.mean, theta = test.mean, num.opt = FALSE)

  test.seq <- seq(from = 0.1, to = 25, by = 0.1)

  dens1.v <- dens1(test.seq)
  dens2.v <- dens2(test.seq)

  expect_equal(dens1.v, dens2.v)

  expect_equal(dens1.v, dexp(test.seq, rate = 1/test.mean))

})

