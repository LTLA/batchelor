# This tests the functions related to multiBatchPCA.
# library(Scratch); library(testthat); source("test-multi-pca.R")

expect_equal_besides_sign <- function(left, right, ...) {
    ratio <- left/right
    right2 <- sweep(right, 2, sign(ratio[1,]), FUN="*")
    expect_equal(left, right2, ...)
}

set.seed(1200001)
test_that("multi-sample PCA works as expected", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)

    # Checking that we get the same result independent of the multiples of the data. 
    ref <- scran:::.multi_pca(list(test1, test2), d=5)
    out <- scran:::.multi_pca(list(test1, cbind(test2, test2)), d=5)
    expect_equal_besides_sign(ref[[1]], out[[1]])
    expect_equal_besides_sign(rbind(ref[[2]], ref[[2]]), out[[2]])
    expect_identical(ncol(ref[[1]]), 5L)
    expect_identical(ncol(ref[[2]]), 5L)

    # Another check. 
    ref <- scran:::.multi_pca(list(test1, test2), d=3)
    out <- scran:::.multi_pca(list(cbind(test1, test1, test1), test2), d=3)
    expect_equal_besides_sign(ref[[2]], out[[2]])
    expect_equal_besides_sign(rbind(ref[[1]], ref[[1]], ref[[1]]), out[[1]])
    expect_identical(ncol(ref[[1]]), 3L)
    expect_identical(ncol(ref[[2]]), 3L)

    # Checking with equal numbers of cells - should be equivalent to cbind'd PCA. 
    test3 <- matrix(rnorm(1000), nrow=10)
    out <- scran:::.multi_pca(list(test1, test3), d=4)
    ref <- prcomp(t(cbind(test1, test3)), rank.=4)
    expect_equal_besides_sign(rbind(out[[1]], out[[2]]), unname(ref$x))
    
    # Checking that the distances match up to the original values when we use full rank.
    out <- scran:::.multi_pca(list(test1, test2), d=nrow(test1))

    everything <- rbind(t(test1), t(test2))
    ref.dist <- as.matrix(dist(everything))
    out.dist <- as.matrix(dist(rbind(out[[1]], out[[2]])))
    expect_equal(ref.dist, out.dist)
    
    out2 <- scran:::.multi_pca(list(test1, cbind(test2, test2)), d=nrow(test1))
    out.dist2 <- as.matrix(dist(do.call(rbind, as.list(out2))))
    everything2 <- rbind(everything, t(test2))
    ref.dist2 <- as.matrix(dist(everything2)) 
    expect_equal(out.dist2, ref.dist2)
})

set.seed(1200001)
test_that("multi-sample PCA works with high-level options", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)

    # Subsetting works correctly.
    i <- sample(nrow(test1), 5)
    ref <- multiBatchPCA(test1, test2, d=3, subset.row=i)
    out <- multiBatchPCA(test1[i,], test2[i,], d=3)
    expect_equal(ref, out)

    # Respects names.
    out <- multiBatchPCA(A=test1, V=test2, d=3)
    expect_identical(names(out), c("A", "V"))

    # Behaves with SCE objects as input.
    sce1 <- SingleCellExperiment(list(logcounts=test1))
    sce2 <- SingleCellExperiment(list(logcounts=test2))
    ref <- multiBatchPCA(sce1, sce2, d=4)
    out <- multiBatchPCA(test1, test2, d=4)
    expect_equal(ref, out)

    # Behaves with spikes as input.
    isp <- sample(nrow(test1), 2)
    isSpike(sce1, "ERCC") <- isp
    isSpike(sce2, "ERCC") <- isp
    ref <- multiBatchPCA(sce1, sce2, d=5)
    out <- multiBatchPCA(test1[-isp,], test2[-isp,], d=5)
    expect_equal(ref, out)

    # Spikes and subsetting interact correctly.
    ref <- multiBatchPCA(sce1, sce2, d=2, subset.row=i)
    out <- multiBatchPCA(sce1[i,], sce2[i,], d=2)
    expect_equal(ref, out)

    # Throws a variety of useful errors.
    expect_error(multiBatchPCA(sce1, test2), "cannot mix")
    expect_error(multiBatchPCA(), "at least one batch")
    expect_error(multiBatchPCA(test1, test2[0,,drop=FALSE]), "not the same")
    sce1x <- clearSpikes(sce1)
    expect_error(multiBatchPCA(sce1x, sce2), "spike-in sets")
})
