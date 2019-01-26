# This tests the functions related to multiBatchPCA.
# library(batchelor); library(testthat); source("test-multi-pca.R")

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
    ref <- multiBatchPCA(test1, test2, d=5)
    out <- multiBatchPCA(test1, cbind(test2, test2), d=5)
    expect_equal_besides_sign(ref[[1]], out[[1]])
    expect_equal_besides_sign(rbind(ref[[2]], ref[[2]]), out[[2]])
    expect_identical(ncol(ref[[1]]), 5L)
    expect_identical(ncol(ref[[2]]), 5L)

    # Another check. 
    ref <- multiBatchPCA(test1, test2, d=3)
    out <- multiBatchPCA(cbind(test1, test1, test1), test2, d=3)
    expect_equal_besides_sign(ref[[2]], out[[2]])
    expect_equal_besides_sign(rbind(ref[[1]], ref[[1]], ref[[1]]), out[[1]])
    expect_identical(ncol(ref[[1]]), 3L)
    expect_identical(ncol(ref[[2]]), 3L)

    # Checking with equal numbers of cells - should be equivalent to cbind'd PCA. 
    test3 <- matrix(rnorm(1000), nrow=10)
    out <- multiBatchPCA(test1, test3, d=4)
    ref <- prcomp(t(cbind(test1, test3)), rank.=4)
    expect_equal_besides_sign(rbind(out[[1]], out[[2]]), unname(ref$x))
    
    # Checking that the distances match up to the original values when we use full rank.
    out <- multiBatchPCA(test1, test2, d=nrow(test1))

    everything <- rbind(t(test1), t(test2))
    ref.dist <- as.matrix(dist(everything))
    out.dist <- as.matrix(dist(rbind(out[[1]], out[[2]])))
    expect_equal(ref.dist, out.dist)
    
    out2 <- multiBatchPCA(test1, cbind(test2, test2), d=nrow(test1))
    out.dist2 <- as.matrix(dist(do.call(rbind, as.list(out2))))
    everything2 <- rbind(everything, t(test2))
    ref.dist2 <- as.matrix(dist(everything2)) 
    expect_equal(out.dist2, ref.dist2)

    # Throws an error on mismatched dimensions.
    expect_error(multiBatchPCA(test1, test2[0,,drop=FALSE]), "not the same")
})

set.seed(12000011)
test_that("multi-sample PCA works with subsetting", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)

    i <- sample(nrow(test1), 5)
    ref <- multiBatchPCA(test1, test2, d=3, subset.row=i)
    out <- multiBatchPCA(test1[i,], test2[i,], d=3)
    expect_equal(ref, out)

    l <- rbinom(nrow(test1), 1, 0.5)==1
    ref <- multiBatchPCA(test1, test2, d=3, subset.row=l)
    out <- multiBatchPCA(test1[l,], test2[l,], d=3)
    expect_equal(ref, out) 
})

set.seed(12000012)
test_that("multi-sample PCA works with SCEs", {
    test1 <- matrix(rnorm(1000), nrow=8)
    test2 <- matrix(rnorm(2000), nrow=8)
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
    i <- sample(nrow(test1), 5)
    ref <- multiBatchPCA(sce1, sce2, d=2, subset.row=i)
    out <- multiBatchPCA(sce1[i,], sce2[i,], d=2)
    expect_equal(ref, out)

    # Throws a variety of useful errors.
    expect_error(multiBatchPCA(sce1, test2), "cannot mix")
    sce1x <- clearSpikes(sce1)
    expect_error(multiBatchPCA(sce1x, sce2), "spike-in sets")
})

set.seed(1200002)
test_that("multi-sample PCA rotation vector inference works", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)
    N <- nrow(test1)

    ref <- multiBatchPCA(test1, test2, d=5)
    others <- c(7, 3, 1, 2)
    expanded <- c(seq_len(N), others)
    out <- multiBatchPCA(test1[expanded,], test2[expanded,], d=5, subset.row=1:N, rotate.all=TRUE)

    expect_equal(ref[[1]], out[[1]])
    expect_equal(ref[[2]], out[[2]])
    expect_equal(metadata(ref)$rotation[expanded,], metadata(out)$rotation)

    # Also works with a logical vector for subsetting. 
    alt <- multiBatchPCA(test1[expanded,], test2[expanded,], d=5, subset.row=seq_along(expanded) <= N, rotate.all=TRUE)
    expect_equal(alt, out)
})

set.seed(12000021)
test_that("multi-sample PCA works with a single input", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)
    ref <- multiBatchPCA(test1, test2, d=5)

    combined <- cbind(test1, test2)
    batch <- rep(1:2, c(ncol(test1), ncol(test2)))
    out <- multiBatchPCA(combined, batch=batch, d=5)
    expect_equal(rbind(ref[[1]], ref[[2]]), out[[1]])
    expect_equal(metadata(out), metadata(ref))

    # Handles reordering.
    o <- sample(length(batch))
    alt <- multiBatchPCA(combined[,o,drop=FALSE], batch=batch[o], d=5)
    expect_equal(out[[1]][o,], alt[[1]])
    expect_equal(metadata(out), metadata(alt))

    # Expected errors.
    expect_error(multiBatchPCA(), "at least one batch")
    expect_error(multiBatchPCA(test1), "'batch' must be specified")
})

set.seed(1200003)
test_that("multi-sample PCA preserves dimension names in the output", {
    test1 <- matrix(rnorm(1000), nrow=10)
    test2 <- matrix(rnorm(2000), nrow=10)

    # Respects batch names themselves.
    out <- multiBatchPCA(A=test1, V=test2, d=3)
    expect_identical(names(out), c("A", "V"))

    # Preserves row names in rotation vectors.
    rownames(test1) <- rownames(test2) <- LETTERS[sample(nrow(test1))]
    out <- multiBatchPCA(test1, test2, d=3)
    expect_identical(rownames(test1), rownames(metadata(out)$rotation))

    out <- multiBatchPCA(test1, test2, d=3, subset.row=2:6)
    expect_identical(rownames(test1)[2:6], rownames(metadata(out)$rotation))

    out <- multiBatchPCA(test1, test2, d=3, subset.row=2:6, rotate.all=TRUE)
    expect_identical(rownames(test1), rownames(metadata(out)$rotation))

    # Preserves column names.
    test1 <- matrix(rnorm(1000), nrow=10)
    colnames(test1) <- sprintf("CELL_%i", seq_len(ncol(test1)))
    test2 <- matrix(rnorm(2000), nrow=10)
    colnames(test2) <- sprintf("CELL2_%i", seq_len(ncol(test2)))
    
    out <- multiBatchPCA(test1, test2, d=3)
    expect_identical(rownames(out[[1]]), colnames(test1))
    expect_identical(rownames(out[[2]]), colnames(test2))

    # ... even when they're non-unique.
    test2 <- matrix(rnorm(2000), nrow=10)
    colnames(test2) <- sprintf("CELL_%i", seq_len(ncol(test2)))
    out <- multiBatchPCA(test1, test2, d=3)
    expect_identical(rownames(out[[1]]), colnames(test1))
    expect_identical(rownames(out[[2]]), colnames(test2))
})


