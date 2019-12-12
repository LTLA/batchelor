# Tests the behaviour of the regressBatches function.
# library(batchelor); library(testthat); source("test-regress-batch.R")

expect_equal_matrix <- function(x, y) {
    # A ResidualMatrix isn't technically the same
    # after subsetting, so we need to compare the
    # matrices rather than the S4 objects.
    expect_equal(as.matrix(x), as.matrix(y))
}

set.seed(130000)
test_that("regressBatches works correctly", {
    means <- 2^rgamma(1000, 2, 1)
    A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2

    ave1 <- rowMeans(A1)
    ave2 <- rowMeans(A2)
    ref <- pmin(ave1, ave2)

    B1 <- log2(A1 + 1)
    B2 <- log2(A2 + 1)

    # Checking against the reference calculations.
    corrected <- regressBatches(B1, B2)
    expect_identical(as.integer(corrected$batch), rep(1:2, c(ncol(B1), ncol(B2))))
    expect_s4_class(assay(corrected), "ResidualMatrix")

    left <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==1)]
    right <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==2)]
    expect_equal_matrix(left, B1 - rowMeans(B1))
    expect_equal_matrix(right, B2 - rowMeans(B2))

    # Behaves with subsetting.
    keep <- sample(nrow(B1), nrow(B1)/2)
    subsetted <- regressBatches(B1, B2, subset.row=keep)
    corrected <- regressBatches(B1[keep,], B2[keep,])
    expect_equal(subsetted, corrected)
}) 

set.seed(130000)
test_that("regressBatches behaves correctly with a design matrix", {
    means <- 2^rgamma(1000, 2, 1)
    A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2

    B1 <- log2(A1 + 1)
    B2 <- log2(A2 + 1)
    ref <- regressBatches(B1, B2)

    b <- rep(1:2, c(ncol(A1), ncol(A2)))
    corrected <- regressBatches(B1, B2, design=model.matrix(~factor(b)))
    expect_equal(as.matrix(assay(ref)), as.matrix(assay(corrected)))

    # Testing against continuous covariates.
    combined <- cbind(B1, B2)

    stuff <- rnorm(ncol(combined))
    corrected2 <- regressBatches(B1, B2, design=model.matrix(~b))
    expect_equivalent(as.matrix(assay(corrected2)),
        t(lm.fit(t(combined), x=model.matrix(~b))$residual))

    # Throws an error.
    expect_error(regressBatches(B1, B2, design=cbind(1)), "total number")
})

set.seed(130001)
test_that("regressBatches works correctly with SCE inputs", {
    A1 <- matrix(rpois(10000, lambda=10), nrow=100) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=10), nrow=100) # Batch 2
    B1 <- log2(A1+1)
    B2 <- log2(A2+1)

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    out <- regressBatches(sce1, sce2)
    ref <- regressBatches(B1, B2)
    expect_equal(ref, out)

    # Subsetting works correctly.
    i <- rbinom(nrow(B1), 1, 0.5)==1L
    out <- regressBatches(sce1, sce2, subset.row=i)
    ref <- regressBatches(sce1[i,], sce2[i,])
    expect_equal(ref, out)
})

set.seed(1300011)
test_that("regressBatches reports names correctly", {
    A1 <- matrix(rpois(10000, lambda=10), nrow=100) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=10), nrow=100) # Batch 2
    B1 <- log2(A1+1)
    B2 <- log2(A2+1)
    
    rownames(B1) <- rownames(B2) <- sprintf("GENE_%i", seq_len(nrow(A1)))
    out <- regressBatches(B1, B2)
    expect_identical(rownames(out), rownames(B1))
    expect_identical(colnames(out), NULL)

    colnames(B1) <- sprintf("Cell_%i", seq_len(ncol(B1)))
    out <- regressBatches(B1, B2)
    expect_identical(colnames(out), c(colnames(B1), character(ncol(B2))))

    colnames(B2) <- sprintf("Yay_%i", seq_len(ncol(B2)))
    out <- regressBatches(B1, B2)
    expect_identical(colnames(out), c(colnames(B1), colnames(B2)))

    # Handles names of batches.
    out2 <- regressBatches(A=B1, B=B2)
    expect_identical(LETTERS[out$batch], out2$batch)

    # Works correctly upon subsetting.
    out <- regressBatches(B1, B2, subset.row=1:10)
    ref <- regressBatches(B1[1:10,], B2[1:10,])
    expect_identical(assay(out), assay(ref))
    expect_identical(out$batch, ref$batch)
})

set.seed(130002)
test_that("regressBatches works with within-object batches", {
    A1 <- matrix(rpois(10000, lambda=5), nrow=100) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=5), nrow=100) # Batch 2
    A3 <- matrix(rpois(15000, lambda=5), nrow=100) # Batch 2
    B1 <- log2(A1+1)
    B2 <- log2(A2+1)
    B3 <- log2(A3+1)

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    sce3 <- SingleCellExperiment(list(logcounts=B3))
    combined <- cbind(sce1, sce2, sce3)
    batches <- rep(1:3, c(ncol(sce1), ncol(sce2), ncol(sce3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    ref <- regressBatches(B1, B2, B3)
    out <- regressBatches(combined, batch=batches)
    expect_equal_matrix(assay(ref)[,shuffle], assay(out))
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # Same for matrix input.        
    out2 <- regressBatches(logcounts(combined), batch=batches)
    expect_equal(out, out2)

    # Preserves column names.
    combined2 <- combined
    colnames(combined2) <- sprintf("Cell_%i", shuffle)
    out3 <- regressBatches(combined2, batch=batches)
    expect_equal(colnames(combined2), colnames(out3))
})

set.seed(1300021)
test_that("regressBatches works with restrictions", {
    A1 <- matrix(rpois(10000, lambda=5), nrow=100) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=5), nrow=100) # Batch 2
    B1 <- log2(A1+1)
    B2 <- log2(A2+1)
    ref <- regressBatches(B1, B2)

    C1 <- cbind(B1, B1[,1:10])
    C2 <- cbind(B2, B2[,1:20])
    extras1 <- ncol(B1) + 1:10
    extras2 <- ncol(B2) + 1:20
    test <- regressBatches(C1, C2, restrict=list(-extras1, -extras2))

    expect_equal_matrix(assay(ref), assay(test)[,c(-extras1, -(ncol(C1) + extras2))])
    expect_equal_matrix(assay(ref)[,1:10], assay(test)[,extras1])
    expect_equal_matrix(assay(ref)[,ncol(B1) + 1:20], assay(test)[,ncol(C1) + extras2])

    expect_error(regressBatches(C1, C2, restrict=list(integer(0), integer(0))), "no cells remaining")

    # Repeating with a single batch.
    D <- cbind(C1, C2)
    batches <- rep(LETTERS[1:2], c(ncol(C1), ncol(C2)))

    to_remove <- logical(ncol(D))
    to_remove[c(extras1, ncol(C1) + extras2)] <- TRUE
    test <- regressBatches(D, batch=batches, restrict=list(!to_remove))

    ref <- regressBatches(D[,!to_remove], batch=batches[!to_remove])
    expect_equal_matrix(assay(ref), assay(test)[,!to_remove])
    expect_identical(ref$batch, test$batch[!to_remove])
    expect_equal_matrix(assay(ref)[,1:10], assay(test)[,extras1])
    expect_equal_matrix(assay(ref)[,ncol(B1) + 1:20], assay(test)[,ncol(C1) + extras2])

    # With shuffling of the cells.
    shuffle <- sample(ncol(D))
    test2 <- regressBatches(D[,shuffle], batch=batches[shuffle], restrict=list(!to_remove[shuffle]))
    expect_equal_matrix(assay(test)[,shuffle], assay(test2))
})

set.seed(130003)
test_that("regressBatches fails on silly inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    # Throws errors properly with no genes or no cells.
    expect_error(regressBatches(), "at least two batches")
    expect_error(regressBatches(B1), "'batch' must be specified")

    # Throws errors upon row checks.
    expect_error(regressBatches(B1[1:10,], B2), "number of rows is not the same")
    xB1 <- B1
    xB2 <- B2
    rownames(xB1) <- sample(nrow(B1))
    rownames(xB2) <- sample(nrow(B2))
    expect_error(regressBatches(xB1, xB2), "row names are not the same")
})
