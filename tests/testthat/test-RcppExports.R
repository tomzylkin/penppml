test_that("Faster matrix multiplication works", {
  A <- matrix(rnorm(n = 20), nrow = 5)
  B <- matrix(rnorm(n = 20), nrow = 4)

  expect_equal(A %*% B, eigenMatMult(A, B))
  # TODO: learn why the following line doesn't work (something to do with the Eigen::Map class)
  #expect_equal(A %*% B, eigenMapMatMult(A, B))
})

# TODO: learn why the following test doesn't work (something to do with the Eigen::Map class)

# test_that("Faster A'A works", {
#  A <- matrix(1:20, nrow = 4)
#  expect_equal(t(A) %*% A, AtA(A))
#})

test_that("Fast OLS works", {

  V <- matrix(runif(n = 25, min = -1, max = 1), nrow = 5, ncol = 5)
  V <- t(V) %*% V

  X <- MASS::mvrnorm(n = 1000, mu = rep(0, times = 5), Sigma = V)
  beta <- c(0, 1, 0, 2, 3)
  e <- rnorm(n = 1000)

  y <- X %*% beta + e

  slow_ols <- lm.fit(x = X, y = y)

  expect_equal(length(fastolsCpp(X, y)), ncol(X))
  expect_true(is.numeric(fastolsCpp(X, y)))
  expect_equal(fastolsCpp(X, y), as.vector(solve(t(X) %*% X) %*% t(X) %*% y))
  expect_equal(fastolsCpp(X, y), as.vector(slow_ols$coefficients))
})

test_that("Fast ridge regression works", {

  V <- matrix(runif(n = 25, min = -1, max = 1), nrow = 5, ncol = 5)
  V <- t(V) %*% V

  X <- MASS::mvrnorm(n = 1000, mu = rep(0, times = 5), Sigma = V)
  beta <- c(0, 1, 0, 2, 3)
  e <- rnorm(n = 1000)

  y <- X %x% beta + e

  expect_equal(length(fastridgeCpp(X, y, 0)), ncol(X))
  expect_true(is.numeric(fastridgeCpp(X, y, 0)))
  expect_equal(fastridgeCpp(X, y, 0), fastolsCpp(X, y))
})

test_that("Fast standard deviation works", {

  V <- matrix(runif(n = 25, min = -1, max = 1), nrow = 5, ncol = 5)
  V <- t(V) %*% V

  X <- MASS::mvrnorm(n = 1000, mu = rep(0, times = 5), Sigma = V)
  beta <- c(0, 1, 0, 2, 3)
  e <- rnorm(n = 1000)

  y <- X %x% beta + e

  w = rep(1, nrow(X))

  expect_equal(length(faststddev(X, w)), ncol(X))
  expect_true(is.numeric(faststddev(X, w)))
})

