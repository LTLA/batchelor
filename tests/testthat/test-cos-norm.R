# Checks the application of cosine normalization.
# require(Scratch); require(testthat); source("test-cos-norm.R")

set.seed(10001)
test_that("Cosine normalization is correct", {
    X <- matrix(rnorm(10000), ncol=100)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, cosineNorm(X))

    X <- matrix(rpois(20000, lambda=5), ncol=1000)
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, cosineNorm(X))
})

set.seed(10002)
test_that("Cosine normalization behaves with silly inputs", {
    X <- matrix(rnorm(10000), ncol=100)
    expect_equal(cosineNorm(X[,0]), matrix(0, nrow(X), 0))
    expect_equal(cosineNorm(X[0,]), matrix(0, 0, ncol(X)))
    expect_equal(cosineNorm(X[,0], mode="l2norm"), numeric(0))
    expect_equal(cosineNorm(X[0,], mode="l2norm"), numeric(ncol(X)))
})
