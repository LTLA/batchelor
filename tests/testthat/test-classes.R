# Tests the various classes.
# library(testthat); library(batchelor); source("test-classes.R")

test_that("class construction works as expected", {
    m <- list(1,list(3,4))
    out <- NoCorrectParam(merge.order=m)
    expect_identical(out$merge.order, m) # protects list from SimpleList constructor.

    out <- NoCorrectParam(merge.order=m, d=20)
    expect_identical(out$merge.order, m)
    expect_identical(20, out$d) 
})
