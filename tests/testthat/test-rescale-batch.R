# Tests the behaviour of the rescaleBatches function.
# library(batchelor); library(testthat); source("test-rescale-batch.R")

set.seed(130000)
test_that("rescaleBatches works correctly", {
    means <- 2^rgamma(1000, 2, 1)
    A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2

    ave1 <- rowMeans(A1)
    ave2 <- rowMeans(A2)
    ref <- pmin(ave1, ave2)

    B1 <- log2(A1 + 1)
    B2 <- log2(A2 + 1)

    # Checking against the reference calculations.
    corrected <- rescaleBatches(B1, B2)
    expect_identical(as.integer(corrected$batch), rep(1:2, c(ncol(B1), ncol(B2))))

    left <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==1)]
    right <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==2)]
    expect_equal(left, log2(A1 * ref/ave1 + 1))
    expect_equal(right, log2(A2 * ref/ave2 + 1))

    # Behaves with lots of zeroes.
    B1[1:10,] <- 0
    B2[5:15,] <- 0
    corrected <- rescaleBatches(B1, B2)
    expect_equal(assay(corrected, withDimnames=FALSE)[1:15,], matrix(0, 15, ncol(B1)+ncol(B2)))
})

set.seed(1300001)
test_that("rescaleBatches responds to the various options", {
    means <- 2^rgamma(1000, 2, 1)
    A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2

    ave1 <- rowMeans(A1)
    ave2 <- rowMeans(A2)
    ref <- pmin(ave1, ave2)

    B1 <- log2(A1 + 1)
    B2 <- log2(A2 + 1)

    # Behaves with subsetting.
    keep <- sample(nrow(B1), nrow(B1)/2)
    subsetted <- rescaleBatches(B1, B2, subset.row=keep)
    corrected <- rescaleBatches(B1[keep,], B2[keep,])
    expect_equal(subsetted, corrected)

    # Behaves with different log-bases.
    rescaled <- rescaleBatches(B1, B2)
    corrected <- rescaleBatches(B1/log2(10), B2/log2(10), log.base=10)
    expect_equal(assay(rescaled)/log2(10), assay(corrected))

    # Behaves with different pseudo.counts.
    C1 <- log10(A1 + 3.2)
    C2 <- log10(A2 + 3.2)
    corrected <- rescaleBatches(C1, C2, pseudo.count=3.2, log.base=10)
    left <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==1)]
    right <- assay(corrected, withDimnames=FALSE)[,BiocGenerics::which(corrected$batch==2)]
    expect_equal(left, log10(A1 * ref/ave1 + 3.2))
    expect_equal(right, log10(A2 * ref/ave2 + 3.2))
}) 

set.seed(1300002)
test_that("rescaleBatches preserves sparsity if possible", {
    means <- 2^rgamma(1000, 2, 1)/10
    A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=means*runif(1000, 0, 2)), ncol=100) # Batch 2
    B1 <- log2(A1 + 1)
    B2 <- log2(A2 + 1)

    library(Matrix)
    B1s <- as(B1, "dgCMatrix")
    B2s <- as(B2, "dgCMatrix")

    rescaled <- rescaleBatches(B1s, B2s)
    expect_s4_class(assay(rescaled), "dgCMatrix")

    ref <- rescaleBatches(B1, B2) 
    expect_equivalent(assay(ref), as.matrix(assay(rescaled)))
})

set.seed(130001)
test_that("rescaleBatches works correctly with SCE inputs", {
    A1 <- matrix(rpois(10000, lambda=10), nrow=100) # Batch 1 
    A2 <- matrix(rpois(20000, lambda=10), nrow=100) # Batch 2
    B1 <- log2(A1+1)
    B2 <- log2(A2+1)

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
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

    ref <- rescaleBatches(B1, B2, B3)
    out <- rescaleBatches(combined, batch=batches)
    expect_equal(assay(ref)[,shuffle], assay(out))
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # Same for matrix input.        
    out2 <- rescaleBatches(logcounts(combined), batch=batches)
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


