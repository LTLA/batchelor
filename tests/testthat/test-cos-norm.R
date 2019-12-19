# Checks the application of cosine normalization.
# require(batchelor); require(testthat); source("setup.R"); source("test-cos-norm.R")

set.seed(10001)
test_that("Cosine normalization is correct", {
    X <- matrix(rnorm(10000), ncol=100)
    cellnorm <- sqrt(colSums(X^2))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, cosineNorm(X))

    X <- matrix(rpois(20000, lambda=5), ncol=1000)
    cellnorm <- sqrt(colSums(X^2))
    ref <- X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
    expect_equal(ref, cosineNorm(X))

    # Recovering both or either works.
    out <- cosineNorm(X, mode="all")
    expect_equal(out$matrix, ref)
    expect_equal(out$l2norm, cellnorm)

    out <- cosineNorm(X, mode="l2norm")
    expect_equal(out, cellnorm)
})

set.seed(100011)
test_that("Cosine normalization preserves sparsity", {
    library(Matrix)
    X <- rsparsematrix(100, 20, 0.2)
    out <- cosineNorm(X)
    expect_s4_class(out, "dgCMatrix")

    ref <- cosineNorm(as.matrix(X))
    expect_equivalent(ref, as.matrix(out))
    
    # We get the l2 norm correctly.
    out <- cosineNorm(X, mode="l2norm")
    ref <- cosineNorm(as.matrix(X), mode="l2norm")
    expect_equivalent(ref, out)
})

test_that("Cosine normalization avoids division by zero", {
    X <- matrix(0, 100, 20)
    out <- cosineNorm(X, mode="all")
    expect_equal(out$matrix, X)
    expect_equal(out$l2norm, numeric(ncol(X)))
})

set.seed(100012)
test_that("Cosine normalization preserves Delayed'ness", {
    X <- DelayedArray(matrix(rnorm(2000), nrow=100))
    out <- cosineNorm(X)
    ref <- cosineNorm(as.matrix(X))
    expect_equal(as.matrix(out), ref)

    # We get the l2 norm correctly.
    out <- cosineNorm(X, mode="l2norm")
    ref <- cosineNorm(as.matrix(X), mode="l2norm")
    expect_equivalent(ref, out)
})

set.seed(100013)
test_that("Cosine normalization preserves dimension names", {
    X <- matrix(rnorm(2000), nrow=100)
    dimnames(X) <- list(sprintf("GENE_%i", seq_len(nrow(X))), 
        sprintf("CELL_%i", seq_len(ncol(X))))
    out <- cosineNorm(X)
    expect_identical(dimnames(X), dimnames(out))

    Y <- DelayedArray(X)
    expect_identical(dimnames(X), dimnames(Y))
    out <- cosineNorm(Y)
    expect_identical(dimnames(Y), dimnames(out))
})

set.seed(100014)
test_that("Cosine normalization works with subsetting", {
    X <- matrix(rnorm(2000), nrow=100)
    expect_identical(cosineNorm(X, subset.row=1:10), cosineNorm(X[1:10,]))
    expect_identical(cosineNorm(X, subset.row=1:10, mode="l2norm"), cosineNorm(X[1:10,], mode="l2norm"))
})

set.seed(10002)
test_that("Cosine normalization behaves with silly inputs", {
    X <- matrix(rnorm(10000), ncol=100)
    expect_equal(cosineNorm(X[,0]), matrix(0, nrow(X), 0))
    expect_equal(cosineNorm(X[0,]), matrix(0, 0, ncol(X)))
    expect_equal(cosineNorm(X[,0], mode="l2norm"), numeric(0))
    expect_equal(cosineNorm(X[0,], mode="l2norm"), numeric(ncol(X)))
})
