context("Objects - density")

test_that("basic tests, univariate",{
  bern.obj <- dist_Bernoulli()
})

test_that("basic tests, multivariate",{
  beta.obj <- dist_Beta()

  tstvec <- 1:9/10

  S.tstvec <- beta.obj$S(tstvec)
  expect_equal(nrow(S.tstvec), length(tstvec))
  expect_equal(ncol(S.tstvec), beta.obj$eta.dim)

  tsteta <- c(2, 4)
  tstdens <- ExpFam_density_eta(beta.obj, tsteta)(tstvec)
  expect_equal(tstdens, dbeta(tstvec, tsteta[1], tsteta[2]))

  tst.1 <- sample(tstvec, 1)
  expect_equal(ExpFam_density_eta(beta.obj, tsteta)(tst.1), dbeta(tst.1, tsteta[1], tsteta[2]))

  tst.2 <- sample(tstvec, beta.obj$eta.dim)
  expect_equal(ExpFam_density_eta(beta.obj, tsteta)(tst.2), dbeta(tst.2, tsteta[1], tsteta[2]))


})

