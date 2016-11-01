test_vec_num_ok <- function(x) {
  cond.ok.full <- !is.na(x) & !is.nan(x) & is.finite(x)

  expect_true(all(!is.na(x)), info = "NA in vector")
  expect_true(all(!is.nan(x)), info = "NaN in vector")
  expect_true(all(is.finite(x)), info = "infinites in vector")

  all(cond.ok.full)
}

test_basic_invariance.ExpFam_dist <- function(dist.obj, valid.eta, valid.theta, valid.values) {
  expect_is(dist.obj, "ExpFam_dist")
  if (has_stat_GLM(dist.obj)) {
    expect_true(dist.obj$stats.glm$validmu(valid.theta), info = "provided valid.theta is invalid")
    expect_true(dist.obj$stats.glm$valideta(valid.eta), info = "provided valid.eta is invalid")
  }
  expect_true(all(dist.obj$h(valid.values) > 0), info = "provided valid.values contain elems with P=0 ")

  # theta -> eta -> theta invariance
  test.eta1 <- dist.obj$eta.from.theta(valid.theta)
  expect_true(all(!is.na(test.eta1) & !is.nan(test.eta1) & is.finite(test.eta1)))
  test.theta <- dist.obj$theta.from.eta(test.eta1)
  expect_equal(test.theta, valid.theta)

  test.theta1 <- dist.obj$theta.from.eta(valid.eta)
  expect_true(all(!is.na(test.theta1) & !is.nan(test.theta1) & is.finite(test.theta1)))
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
    dist.obj.inst <- ExpFam_density_eta(dist.obj, test.eta)
    test.eta.density <- dist.obj.inst(valid.values)
    expect_true(test_vec_num_ok(test.eta.density), info = "ExpFam_density_eta returned invalid results")
    test.eta.density.N <- dist.obj$N.density.eta(x = valid.values, eta = test.eta)
    expect_true(test_vec_num_ok(test.eta.density.N), info = "N.density.eta returned invalid results")
    expect_equal(
      test.eta.density,
      test.eta.density.N
    )

    test.theta <- sample(valid.theta, 1)
    dist.obj.inst <- ExpFam_density_theta(dist.obj, test.theta)
    test.theta.density <- dist.obj.inst(valid.values)
    expect_true(test_vec_num_ok(test.theta.density), info = "ExpFam_density_theta returned invalid results")
    test.theta.density.N <- dist.obj$N.density.theta(x = valid.values, theta = test.theta)
    expect_true(test_vec_num_ok(test.theta.density.N), info = "N.density.theta returned invalid results")

    expect_equal(
      test.theta.density,
      test.theta.density.N
    )

    # theta --> eta, dens.theta(values) == dens.eta(values)
    test.theta <- sample(valid.theta, 1)
    dist.obj.inst1 <- ExpFam_density_theta(dist.obj, test.theta)
    test.theta.density <- dist.obj.inst1(valid.values)
    expect_true(test_vec_num_ok(test.theta.density), info = "ExpFam_density_theta returned invalid results")
    test.eta <- dist.obj$eta.from.theta(test.theta)
    dist.obj.inst2 <- ExpFam_density_eta(dist.obj, test.eta)
    test.eta.density <- dist.obj.inst2(valid.values)
    expect_equal(test.eta.density, test.theta.density)

  }


}

test_basic_invariance.ExpFam_dist_ext <- function(dist.obj, valid.eta, valid.values, param.params) {
  expect_is(dist.obj, "ExpFam_dist_ext")

  par.can <- dist.obj[["canonical"]]
  expect_true(par.can$eta.in.domain(valid.eta), info = "provided valid.eta is invalid")


  for (pname in names(param.params)) {
    dist.obj.par.raw <- dist.obj[[pname]]
    valid.theta.claimed <- param.params[[pname]]$valid.theta

    expect_true(
      dist.obj.par.raw$theta.in.domain(valid.theta.claimed),
      info = "provided valid.theta is invalid"
      )

    for (num.opt in c(TRUE, FALSE)) {
      dist.obj.par <- ExpFam_bind_parametrization(dist.obj, pname, num.opt = num.opt)
      test_basic_invariance.ExpFam_dist(
        dist.obj.par,
        valid.eta = valid.eta,
        valid.values = valid.values,
        valid.theta = valid.theta.claimed
      )
    }
  }
}




