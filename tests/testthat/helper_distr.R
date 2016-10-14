test_basic_invariance <- function(dist.obj, valid.eta, valid.theta, valid.values) {
  expect_is(dist.obj, "ExpFam_dist")
  if (has_stat_GLM(dist.obj)) {
    expect_true(dist.obj$stats.glm$validmu(valid.theta))
    expect_true(dist.obj$stats.glm$valideta(valid.eta))
  }
  expect_true(all(dist.obj$h(valid.values) > 0))

  # theta -> eta -> theta invariance
  test.eta1 <- dist.obj$eta.from.theta(valid.theta)
  test.theta <- dist.obj$theta.from.eta(test.eta1)
  expect_equal(test.theta, valid.theta)

  test.theta1 <- dist.obj$theta.from.eta(valid.eta)
  test.eta <- dist.obj$eta.from.theta(test.theta1)
  expect_equal(test.eta, valid.eta)

  if (has_stat_GLM(dist.obj)) {
    #external checking for validity of generated vectors
    expect_true(dist.obj$stats.glm$validmu(test.theta1))
    expect_true(dist.obj$stats.glm$valideta(test.eta1))
  }

  if (is(dist.obj, "N_ExpFam_dist")) {
    #checking numerical approx
    test.eta <- sample(valid.eta, 1)
    expect_equal(
      ExpFam_density_eta(dist.obj, test.eta)(valid.values),
      dist.obj$N.density.eta(x = valid.values, eta = test.eta)
    )
    test.theta <- sample(valid.theta, 1)
    expect_equal(
      ExpFam_density_theta(dist.obj, test.theta)(valid.values),
      dist.obj$N.density.theta(x = valid.values, theta = test.theta)
    )

  }


}
