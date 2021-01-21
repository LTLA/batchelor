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
    # Using fixed scaling for exact checks.
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

    # Checking it returns the same object as logNormCounts().
    solo.out <- multiBatchNorm(X3, X3)
    ref <- scuttle::logNormCounts(X3)
    expect_equal(solo.out[[1]], ref)
    expect_equal(solo.out[[2]], ref)

    # Reverts to the library size correctly.
    Xtmp <- X
    sizeFactors(Xtmp) <- NULL
    out <- multiBatchNorm(Xtmp, Xtmp)
    ref <- scuttle::librarySizeFactors(Xtmp)
    expect_equal(sizeFactors(out[[1]]), ref)
    expect_equal(sizeFactors(out[[2]]), ref)
})

set.seed(200111)
test_that("multiBatchNorm size factor centering logic is correct", {
    # Tests with even numbers of genes are unnecessarily fragile due to the
    # lack of symmetry for the median (i.e., med(x/y)!=1/med(y/x)). This causes
    # slight differences that break the conceptual tests below. 
    ngenes <- 999
    X <- X[seq_len(ngenes),]

    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Centers of new size factors are correct.
    out <- multiBatchNorm(X, X2, min.mean=0)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(sizeFactors(X)/mean(sizeFactors(X)), sizeFactors(out[[1]]))

    ave1 <- scuttle::calculateAverage(X)
    ave2 <- scuttle::calculateAverage(X2)
    M <- median(ave2/ave1)
    expect_equal(mean(sizeFactors(out[[2]])), M)
    expect_equal(sizeFactors(X2)/mean(sizeFactors(X2)), sizeFactors(out[[2]]) / M)

    # Normalized values have no composition biases.
    ave1 <- scuttle::calculateAverage(out[[1]])
    ave2 <- scuttle::calculateAverage(out[[2]]) / M # as calculateAverage automatically centers out[[2]]'s SFs.
    expect_equal(1, median(ave2/ave1))
})

set.seed(20012)
test_that("multiBatchNorm behaves correctly with mean filtering", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))

    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Creating a reference function using calculateAverage() explicitly. This
    # compares all other batches to the first batch, which is assumed to be the
    # lowest-coverage.
    library(scuttle)
    REFFUN <- function(..., min.mean=1) {
	    batches <- list(...)
	    nbatches <- length(batches)
	    collected.ave <- lapply(batches, calculateAverage)

    	collected.ratios <- rep(1, nbatches)
	    first.ave <- collected.ave[[1]]
        for (second in 2:nbatches) {
	        second.ave <- collected.ave[[second]]
            keep <- calculateAverage(cbind(first.ave, second.ave)) >= min.mean
            collected.ratios[second] <- median(first.ave[keep]/second.ave[keep]) 
		}

        # Again, we compute 1/med(first/second) rather than med(second/first),
        # as these are not exactly equal (especially for small numbers of genes).
		1/collected.ratios
    }

    out <- multiBatchNorm(X, X2, min.mean=1)
	refA <- REFFUN(X, X2, min.mean=1)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(mean(sizeFactors(out[[2]])), refA[2])

    out <- multiBatchNorm(X, X2, min.mean=10)
	refB <- REFFUN(X, X2, min.mean=10)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(mean(sizeFactors(out[[2]])), refB[2])

    out <- multiBatchNorm(X, X2, min.mean=100)
	refC <- REFFUN(X, X2, min.mean=100)
    expect_equal(mean(sizeFactors(out[[1]])), 1)
    expect_equal(mean(sizeFactors(out[[2]])), refC[2])
    
    expect_false(isTRUE(all.equal(refA, refB)))
    expect_false(isTRUE(all.equal(refA, refC)))
})

set.seed(20013)
test_that("multiBatchNorm behaves correctly with subsetting", {
    dummy2 <- matrix(rnbinom(ngenes*ncells, mu=means * 10, size=5), ncol=ncells, nrow=ngenes)
    rownames(dummy2) <- paste0("X", seq_len(ngenes))
    
    X2 <- SingleCellExperiment(list(counts=dummy2))
    sizeFactors(X2) <- runif(ncol(X2))

    # Checking that gene subsetting works correctly.
    randoms <- sample(ngenes, 500)
    ref.out <- multiBatchNorm(X[randoms,], X2[randoms,])

    sub.out <- multiBatchNorm(X, X2, subset.row=randoms)
    expect_equal(logcounts(sub.out[[1]]), logcounts(ref.out[[1]]))
    expect_equal(logcounts(sub.out[[2]]), logcounts(ref.out[[2]]))

    sub.out <- multiBatchNorm(X, X2, subset.row=randoms, normalize.all=TRUE)
    expect_equal(logcounts(sub.out[[1]])[randoms,], logcounts(ref.out[[1]]))
    expect_equal(logcounts(sub.out[[2]])[randoms,], logcounts(ref.out[[2]]))
})

set.seed(20013)
test_that("multiBatchNorm behaves correctly with a single batch", {
    X2 <- X
    counts(X2) <- counts(X2) * runif(nrow(X), 1.5, 2.5)
    X3 <- X
    counts(X3) <- counts(X3) * runif(nrow(X), 2.5, 3.5)

    ref <- multiBatchNorm(X, X2, X3)
    combined <- cbind(X, X2, X3)
    alt <- multiBatchNorm(combined, batch=rep(1:3, each=ncol(X)), preserve.single=FALSE)
    for (i in seq_along(ref)) {
        expect_equal(logcounts(ref[[i]]), logcounts(alt[[i]]))
        expect_equal(sizeFactors(ref[[i]]), sizeFactors(alt[[i]]))
    }
    expect_identical(names(alt), as.character(1:3))

    alt2 <- multiBatchNorm(combined, batch=rep(3:1, each=ncol(X)), preserve.single=FALSE)
    expect_equal(unname(alt2[3:1]), unname(alt))

    # Behaves correctly with subsetting and normalize.all=TRUE.
    reref <- multiBatchNorm(combined[1:100,], batch=rep(1:3, each=ncol(X)), preserve.single=FALSE)
    alt3 <- multiBatchNorm(combined, batch=rep(1:3, each=ncol(X)), subset.row=1:100, preserve.single=FALSE)
    expect_identical(reref, alt3)

    alt4 <- multiBatchNorm(combined, batch=rep(1:3, each=ncol(X)), subset.row=100:1,
        preserve.single=FALSE, normalize.all=TRUE)
    for (i in seq_along(ref)) {
        expect_identical(dim(alt4[[i]]), dim(ref[[i]]))
        expect_equal(logcounts(reref[[i]])[1:100,], logcounts(alt4[[i]])[1:100,])
        expect_equal(sizeFactors(reref[[i]]), sizeFactors(alt4[[i]]))
    }
})

set.seed(20014)
test_that("multiBatchNorm behaves correctly with a single batch and preserving single-ness", {
    X2 <- X
    counts(X2) <- counts(X2) * runif(nrow(X), 1.5, 2.5)
    X3 <- X
    counts(X3) <- counts(X3) * runif(nrow(X), 2.5, 3.5)

    combined <- cbind(X, X2, X3)
    alt <- multiBatchNorm(combined, batch=rep(1:3, each=ncol(X)), preserve.single=FALSE)
    alt2 <- multiBatchNorm(combined, batch=rep(1:3, each=ncol(X)))

    alt <- do.call(cbind, alt)
    int_metadata(alt) <- int_metadata(alt2) 
    expect_equal(alt2, alt)

    # Handles subsetting.
    alt3a <- multiBatchNorm(combined[1:100,], batch=rep(1:3, each=ncol(X)))
    alt3b <- multiBatchNorm(combined, batch=rep(3:1, each=ncol(X)), subset.row=1:100)
    expect_equal(alt3a[1:100,], alt3b)

    alt3c <- multiBatchNorm(combined, batch=rep(3:1, each=ncol(X)), subset.row=100:1, normalize.all=TRUE)
    expect_equal(alt3a[1:100,], alt3c[1:100,])
    expect_identical(dim(alt3c), dim(combined))

    # Handles reordering smoothly.
    o <- sample(ncol(combined))
    alt4 <- multiBatchNorm(combined[,o], batch=rep(1:3, each=ncol(X))[o])
    expect_equal(alt4, alt2[,o])
})

set.seed(20013)
test_that("multiBatchNorm behaves correctly with an input list", {
    X2 <- X
    counts(X2) <- counts(X2) * 2L
    X3 <- X
    counts(X3) <- counts(X3) * 3L

    ref <- multiBatchNorm(X, X2, X3)
    alt1 <- multiBatchNorm(X, list(X2, X3))
    expect_identical(ref, alt1)

    alt2 <- multiBatchNorm(list(X, X2, X3))
    expect_identical(ref, alt2)
})

set.seed(20014)
test_that("multiBatchNorm works correctly with alternative experiments", {
    X2 <- X
    counts(X2) <- counts(X2) * 2L
    X3 <- X
    counts(X3) <- counts(X3) * 3L
    ref <- multiBatchNorm(X, X2, X3)

    dummy1 <- dummy2 <- dummy3 <- X
    altExp(dummy1, "axel") <- X
    altExp(dummy2, "axel") <- X2
    altExp(dummy3, "axel") <- X3
    out <- multiBatchNorm(dummy1, dummy2, dummy3, as.altexp="axel")
    expect_identical(ref, lapply(out, altExp, e="axel"))

    # Same behavior for the combined method.
    ref <- multiBatchNorm(X, batch=rep(1:3, length.out=ncol(X)))
    alt <- multiBatchNorm(dummy1, batch=rep(1:3, length.out=ncol(X)), as.altexp="axel")
    expect_identical(ref, altExp(alt, "axel"))
    
    ref <- multiBatchNorm(X, batch=rep(1:3, length.out=ncol(X)), preserve.single=FALSE)
    alt <- multiBatchNorm(dummy1, batch=rep(1:3, length.out=ncol(X)), as.altexp="axel", preserve.single=FALSE)
    expect_identical(ref, lapply(alt, altExp, e="axel"))
})

