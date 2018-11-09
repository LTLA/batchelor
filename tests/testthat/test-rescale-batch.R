# Tests the behaviour of the rescaleBatches function.
# library(batchelor); library(testthat); source("test-rescale-batch.R")

set.seed(130000)
test_that("rescaleBatches works correctly", {
    means <- 2^rgamma(1000, 2, 1)
    B1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    B2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2

    library(scater)
    ave1 <- calcAverage(B1)
    ave2 <- calcAverage(B2)
    ref <- pmin(ave1, ave2)

    corrected <- rescaleBatches(B1, B2)
    expect_identical(as.integer(corrected$batch), rep(1:2, c(ncol(B1), ncol(B2))))

    left <- assay(corrected, withDimnames=FALSE)[,which(corrected$batch==1)]
    right <- assay(corrected, withDimnames=FALSE)[,which(corrected$batch==2)]
    expect_equal(left, B1 * ref/ave1)
    expect_equal(right, B2 * ref/ave2)

    # Behaves with subsetting.
    keep <- sample(nrow(B1), nrow(B1)/2)
    ref <- rescaleBatches(B1, B2, subset.row=keep)
    corrected <- rescaleBatches(B1[keep,], B2[keep,])
    expect_equal(ref, corrected)
    
    # Behaves with lots of zeroes.
    B1[1:10,] <- 0
    B2[5:15,] <- 0
    corrected <- rescaleBatches(B1, B2)
    expect_equal(assay(corrected, withDimnames=FALSE)[1:15,], matrix(0, 15, ncol(B1)+ncol(B2)))
})

set.seed(130001)
test_that("rescaleBatches works correctly with SCE inputs", {
    B1 <- matrix(rpois(10000, lambda=10), nrow=100) # Batch 1 
    B2 <- matrix(rpois(20000, lambda=10), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(counts=B1))
    sce2 <- SingleCellExperiment(list(counts=B2))
    out <- rescaleBatches(sce1, sce2)
    ref <- rescaleBatches(B1, B2)
    expect_equal(ref, out)

    # Behaves with spikes as input.
    isp <- rbinom(nrow(B1), 1, 0.1)==1L
    isSpike(sce1, "ERCC") <- isp
    isSpike(sce2, "ERCC") <- isp
    out <- rescaleBatches(sce1, sce2)
    ref <- rescaleBatches(B1[!isp,], B2[!isp,])
    expect_equal(ref, out)

    # Spikes and subsetting interact correctly
    i <- rbinom(nrow(B1), 1, 0.5)==1L
    out <- rescaleBatches(sce1, sce2, subset.row=i)
    ref <- rescaleBatches(sce1[i,], sce2[i,])
    expect_equal(ref, out)

    ref2 <- rescaleBatches(B1[!isp & i,], B2[!isp & i,])
    expect_equal(ref2, out)
})

set.seed(130002)
test_that("rescaleBatches works with within-object batches", {
    B1 <- matrix(rpois(10000, lambda=5), nrow=100) # Batch 1 
    B2 <- matrix(rpois(20000, lambda=5), nrow=100) # Batch 2
    B3 <- matrix(rpois(15000, lambda=5), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(counts=B1))
    sce2 <- SingleCellExperiment(list(counts=B2))
    sce3 <- SingleCellExperiment(list(counts=B3))
    combined <- cbind(sce1, sce2, sce3)
    batches <- rep(1:3, c(ncol(sce1), ncol(sce2), ncol(sce3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    ref <- rescaleBatches(B1, B2, B3)
    out <- rescaleBatches(combined, batch=batches)
    expect_equal(assay(ref)[,shuffle], assay(out))
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # Same for matrix input.        
    out2 <- rescaleBatches(counts(combined), batch=batches)
    expect_equal(out, out2)
})

set.seed(130003)
test_that("rescaleBatches fails on silly inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    # Throws errors properly with no genes or no cells.
    expect_error(rescaleBatches(), "at least two batches")
    expect_error(rescaleBatches(B1), "'batch' must be specified")

    # SCE vs matrix errors.    
    expect_error(rescaleBatches(SingleCellExperiment(list(logcounts=B1)), B2), "cannot mix")
 
    # Throws errors upon row checks.
    expect_error(rescaleBatches(B1[1:10,], B2), "number of rows is not the same")
    xB1 <- B1
    xB2 <- B2
    rownames(xB1) <- sample(nrow(B1))
    rownames(xB2) <- sample(nrow(B2))
    expect_error(rescaleBatches(xB1, xB2), "row names are not the same")
})


