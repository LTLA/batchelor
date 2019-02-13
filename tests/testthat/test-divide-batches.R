# Tests the divideIntoBatches function.
# library(batchelor); library(testthat); source("test-divide-batches.R")

set.seed(1000)
test_that("divideIntoBatches works correctly", {
    A <- matrix(runif(2000), 50, 40)

    # Testing by row.
    b <- sample(LETTERS[1:5], nrow(A), replace=TRUE)
    X <- divideIntoBatches(A, b, byrow=TRUE) 
    expect_identical(names(X$batches), sort(unique(b)))

    Y <- do.call(rbind, X$batches)
    expect_identical(Y[X$reorder,,drop=FALSE], A)

    # Testing by column.
    b <- sample(LETTERS[15:6], ncol(A), replace=TRUE)
    X <- divideIntoBatches(A, b)
    expect_identical(names(X$batches), sort(unique(b)))

    Y <- do.call(cbind, X$batches)
    expect_identical(Y[,X$reorder,drop=FALSE], A)

    # Preserves dimnames.
    B <- A
    colnames(B) <- seq_len(ncol(B))
    rownames(B) <- seq_len(nrow(B))
    X <- divideIntoBatches(B, b)
    Y <- do.call(cbind, X$batches)
    expect_identical(Y[,X$reorder,drop=FALSE], B)
})

set.seed(1000)
test_that("divideIntoBatches handles restriction correctly", { 
    A <- matrix(runif(2000), 50, 40)

    # Testing by row.
    b <- sample(LETTERS[1:5], nrow(A), replace=TRUE)
    restrict <- sample(nrow(A), nrow(A)/2)
    X <- divideIntoBatches(A, b, byrow=TRUE, restrict=restrict) 

    expect_identical(X$restricted, lapply(split(seq_len(nrow(A)) %in% restrict, b), which))

    # Testing by row.
    b <- sample(LETTERS[1:5], ncol(A), replace=TRUE)
    restrict <- sample(nrow(A), ncol(A)/2)
    X <- divideIntoBatches(A, b, restrict=restrict) 

    expect_identical(X$restricted, lapply(split(seq_len(ncol(A)) %in% restrict, b), which))

    # Dummy tests.
    X <- divideIntoBatches(A, b, restrict=NULL)
    expect_identical(X$restricted, NULL)

    expect_error(X <- divideIntoBatches(A, b, restrict=integer(0)), "no cells") 
    expect_error(X <- divideIntoBatches(A, b, restrict=b==b[1]), "no cells") 
})

set.seed(1002)
test_that("divideIntoBatches fails correctly", {
    A <- matrix(runif(2000), 50, 40)

    # Testing for errors.
    expect_error(divideIntoBatches(A, NULL), "must be specified")
    expect_error(divideIntoBatches(A, 1), "not the same")
    expect_error(divideIntoBatches(A, 1, byrow=TRUE), "not the same")

    # Testing silly inputs.
    X <- divideIntoBatches(A, rep(1, ncol(A)))
    expect_identical(X$batches[[1]], A)
    X <- divideIntoBatches(A[,0], integer(0))
    expect_identical(length(X$batches), 0L)
    expect_identical(X$reorder, integer(0))
})
