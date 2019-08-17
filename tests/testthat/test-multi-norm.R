# This tests the multi-block normalization calculations.
# require(batchelor); require(testthat); source("test-multi-norm.R")

set.seed(20010)
ncells <- 200
ngenes <- 1000
means <- 2^runif(ngenes, -1, 5)
dummy <- matrix(rnbinom(ngenes*ncells, mu=means, size=5), ncol=ncells, nrow=ngenes)
rownames(dummy) <- paste0("X", seq_len(ngenes))

X <- SingleCellExperiment(list(counts=dummy))
sizeFactors(X) <- runif(ncol(X))

set.seed(20011)
test_that("multiBatchNorm works properly", {
    X2 <- X
    counts(X2) <- counts(X2) * 2L
    X3 <- X
    counts(X3) <- counts(X3) * 3L

    # Checking that it works in the vanilla case.
    out <- multiBatchNorm(X, X2, X3)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[2]])/2)
    expect_equal(sizeFactors(out[[1]]), sizeFactors(out[[3]])/3)
    expect_equal(logcounts(out[[1]]), logcounts(out[[2]]))
    expect_equal(logcounts(out[[1]]), logcounts(out[[3]]))

    # Checking that the order does not matter.
    re.out <- multiBatchNorm(X3, X, X2)
    expect_equal(out[c(3,1,2)], re.out)

    # Single batch input just returns the same object as logNormCounts().
    solo.out <- multiBatchNorm(X3)
    expect_equal(solo.out[[1]], scater::logNormCounts(X3))

    # Reverts to the library size correctly.
    Xtmp <- X
    sizeFactors(Xtmp) <- NULL
    out <- multiBatchNorm(Xtmp)
    expect_equal(sizeFactors(out[[1]]), scater::librarySizeFactors(Xtmp))
})

set.seed(200111)
test_that("multiBatchNorm size factor centering logic is correct", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Centers of new size factors are correct.
    out <- multiBatchNorm(X, X2, min.mean=0)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(sizeFactors(scater::centreSizeFactors(X)), sizeFactors(out[[1]]))

    ave1 <- scater::calcAverage(X)
    ave2 <- scater::calcAverage(X2)
    M <- median(ave2/ave1)
    expect_equal(mean(sizeFactors(out[[2]])), M)
    expect_equal(sizeFactors(scater::centreSizeFactors(X2)), sizeFactors(out[[2]]) / M)

    # RelogNormCountsd values have no composition biases.
    ave1 <- scater::calcAverage(out[[1]])
    ave2 <- scater::calcAverage(out[[2]]) / M # as calcAverage automatically centers out[[2]'s SFs.
    expect_equal(1, median(ave2/ave1))
})

set.seed(20012)
test_that("multiBatchNorm behaves correctly with gene filtering", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Creating a reference function using calcAverage() explicitly.
    library(scater)
    REFFUN <- function(..., min.mean=1) {
	    batches <- list(...)
	    nbatches <- length(batches)
        batches <- lapply(batches, centreSizeFactors)
	    collected.ave <- lapply(batches, calcAverage, use_size_factors=TRUE)

    	collected.ratios <- rep(1, nbatches)
	    first.ave <- collected.ave[[1]]
        for (second in 2:nbatches) {
	        second.ave <- collected.ave[[second]]
            keep <- calcAverage(cbind(first.ave, second.ave)) >= min.mean
            collected.ratios[second] <- median(second.ave[keep]/first.ave[keep]) 
		}
		collected.ratios
    }

    out <- multiBatchNorm(X, X2, min.mean=1)
	refA <- REFFUN(X, X2, min.mean=1)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refA[2])

    out <- multiBatchNorm(X, X2, min.mean=10)
	refB <- REFFUN(X, X2, min.mean=10)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refB[2])

    out <- multiBatchNorm(X, X2, min.mean=100)
	refC <- REFFUN(X, X2, min.mean=100)
    expect_equal(mean(sizeFactors(out[[2]]))/mean(sizeFactors(out[[1]])), refC[2])
    
    expect_false(isTRUE(all.equal(refA, refB)))
    expect_false(isTRUE(all.equal(refA, refC)))

    # Checking that gene subsetting works correctly.
    randoms <- sample(ngenes, 500)
    sub.out <- multiBatchNorm(X, X2, subset.row=randoms)
    ref.out <- multiBatchNorm(X[randoms,], X2[randoms,])
    expect_equal(logcounts(sub.out[[1]])[randoms,], logcounts(ref.out[[1]])) 
    expect_equal(logcounts(sub.out[[2]])[randoms,], logcounts(ref.out[[2]])) 
})
