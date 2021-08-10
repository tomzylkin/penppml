test_that("hdfeppml.int works", {
  V <- matrix(runif(n = 25, min = -1, max = 1), nrow = 5, ncol = 5)
  V <- t(V) %*% V

  X <- MASS::mvrnorm(n = 1000, mu = rep(0, times = 5), Sigma = V)
  X <- X^2                   # Our functions can't handle negative dependent variables.
  beta <- c(0, 1, 0, 2, 3)
  e <- rnorm(n = 1000) ^ 2

  y <- as.vector(X %*% beta + e)

  fes <- list(factor(sample(1:50, size = 1000, replace = TRUE)))

  reg <- hdfeppml_int(y = y, x = X, fes = fes)

  expect_true(is.numeric(reg$coefficients))
})
