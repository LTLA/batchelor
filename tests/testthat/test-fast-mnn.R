# This tests the functions related to fastMNN.
# library(batchelor); library(testthat); source("test-fast-mnn.R")

set.seed(1200001)
test_that("averaging correction vectors works as expected", {
    test1 <- matrix(rnorm(1000), ncol=10)
    test2 <- matrix(rnorm(2000), ncol=10)
    mnn1 <- sample(nrow(test1), 250, replace=TRUE)
    mnn2 <- sample(nrow(test1), 250, replace=TRUE)

    # Slow reference calculation.
    correct <- test1[mnn1,] - test2[mnn2,]
    by.mnn <- split(seq_along(mnn2), mnn2)
    collected <- vector("list", length(by.mnn))
    for (idx in seq_along(by.mnn)) {
        collected[[idx]] <- colMeans(correct[by.mnn[[idx]],,drop=FALSE])
    }
    ref <- do.call(rbind, collected)

    # Comparing the implementation in batchelor.
    out <- batchelor:::.average_correction(test1, mnn1, test2, mnn2)  
    expect_equal(out$averaged, ref)
    expect_identical(out$second, sort(unique(mnn2)))
})

set.seed(1200002)
test_that("centering along a batch vector works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    batch <- rnorm(10) 
    centered <- batchelor:::.center_along_batch_vector(test, batch)
    new.locations <- centered %*% batch
    expect_true(mad(new.locations) < 1e-8)
})

set.seed(1200003)
test_that("tricube weighting works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    correction <- matrix(rnorm(500), ncol=10)
    involved <- sample(nrow(test), nrow(correction))

    # Setting up a reference function for comparison, operating truly row-by-row.
    FUN <- function(current, corvec, in.mnn, k=20, ndist=3) {
        cur.uniq <- current[in.mnn,,drop=FALSE]
        safe.k <- min(k, nrow(cur.uniq))
        closest <- BiocNeighbors::queryKNN(query=current, X=cur.uniq, k=safe.k)
        middle.k <- ceiling(safe.k/2L)

        for (x in seq_len(nrow(current))) {
            all.dists <- closest$distance[x,]
            all.index <- closest$index[x,]

            middist <- sort(all.dists)[middle.k] 
            weights <- (1 - pmin(1, all.dists/(middist*ndist))^3)^3
            weights <- weights/sum(weights)

            curcor <- colSums(corvec[all.index,] * weights) 
            current[x,] <- current[x,] + curcor
        }

        return(current)
    }

    out <- batchelor:::.tricube_weighted_correction(test, correction, involved, k=20, ndist=3)
    ref <- FUN(test, correction, involved, k=20, ndist=3)
    expect_equal(ref, out)

    out <- batchelor:::.tricube_weighted_correction(test, correction, involved, k=11, ndist=3)
    ref <- FUN(test, correction, involved, k=11, ndist=3)
    expect_equal(ref, out)

    out <- batchelor:::.tricube_weighted_correction(test, correction, involved, k=11, ndist=1)
    ref <- FUN(test, correction, involved, k=11, ndist=1)
    expect_equal(ref, out)
})

CHECK_PAIRINGS <- function(mnn.out) 
# Checks that the reported pairs belong to the same batch at each merge step,
# and are consistent with the reported merge order.
{
    origin <- as.vector(mnn.out$batch)
    pairings <- metadata(mnn.out)$pairs
    order <- metadata(mnn.out)$order

    sofar <- unique(origin[pairings[[1]]$first])
    expect_identical(sofar, order[1L])
    counter <- 1L

    for (p in pairings) {
        sbatch <- origin[p$second]
        usecond <- unique(sbatch)
        expect_identical(length(usecond), 1L)

        fbatch <- origin[p$first]
        expect_true(!any(fbatch == usecond))
   
        ufirst <- unique(fbatch) 
        expect_identical(sofar, sort(ufirst))
        sofar <- sort(c(ufirst, usecond))
        expect_identical(sofar, sort(order[seq_along(sofar)]))

        expect_identical(length(sofar), counter + 1L)
        counter <- counter + 1L
    }

    return(NULL)
}

set.seed(1200004)
test_that("fastMNN works as expected for two batches", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, d=50) # corrected values
    expect_identical(dim(out$corrected), c(ncol(B1) + ncol(B2), 50L))
    expect_identical(as.integer(out$batch), rep(1:2, c(ncol(B1), ncol(B2))))
    expect_identical(metadata(out)$order, 1:2)
    CHECK_PAIRINGS(out)
    
    # Dimension choice behaves correctly.
    out.10 <- fastMNN(B1, B2, d=10) 
    expect_identical(ncol(out.10$corrected), 10L)
    CHECK_PAIRINGS(out.10)

    # Handles names correctly.
    out.n <- fastMNN(X=B1, Y=B2, d=50) 
    expect_identical(out$corrected, out.n$corrected)
    expect_identical(as.character(out.n$batch), rep(c("X", "Y"), c(ncol(B1), ncol(B2))))
    CHECK_PAIRINGS(out.n)

    # Behaves if we turn off cosine-norm.
    nB1 <- t(t(B1)/ sqrt(colSums(B1^2)))
    nB2 <- t(t(B2)/ sqrt(colSums(B2^2)))
    out.ncos <- fastMNN(nB1, nB2, cos.norm=FALSE, d=50) 
    expect_equal(out.ncos, out) 

    # Subset.row behaves correctly.
    i <- sample(nrow(B1), 50)
    ref <- fastMNN(X=B1[i,], Y=B2[i,], d=50)
    out.s <- fastMNN(X=B1, Y=B2, d=50, subset.row=i)
    expect_identical(out.s, ref)

    # Behaves if we only use PCs.
    pcs <- multiBatchPCA(B1, B2, d=10)
    out.pre <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
    out.norm <- fastMNN(B1, B2, d=10, cos.norm=FALSE)

    expect_equal(metadata(pcs)$rotation, metadata(out.norm)$rotation)
    metadata(out.norm)$rotation <- NULL
    expect_equal(out.pre, out.norm)
})

set.seed(1200004)
test_that("variance loss calculations work as expected", {
    PC1 <- matrix(rnorm(10000), ncol=10) # Batch 1 
    PC2 <- matrix(rnorm(20000), ncol=10) # Batch 2

    out <- batchelor:::.compute_intra_var(PC1, PC2, list(PC1, PC2), c(1L, 2L))
    expect_identical(out[1], sum(DelayedMatrixStats::colVars(DelayedArray(PC1))))
    expect_identical(out[2], sum(DelayedMatrixStats::colVars(DelayedArray(PC2))))

    # Alternative ordering. 
    out2 <- batchelor:::.compute_intra_var(PC2, PC1, list(PC1, PC2), c(2L, 1L))
    expect_identical(out, rev(out2))

    # Multiple lengths.
    out3 <- batchelor:::.compute_intra_var(rbind(PC1, PC1), PC2, list(PC1, PC2), c(1L, 1L, 2L))
    expect_identical(out3, out[c(1,1,2)])

    # Checking that we compute something.
    mnn.out <- fastMNN(PC1, PC2, pc.input=TRUE, compute.variances=TRUE)
    expect_identical(length(metadata(mnn.out)$lost.var), 2L)
    mnn.out2 <- fastMNN(PC2, PC1, auto.order=2:1, pc.input=TRUE, compute.variances=TRUE)
    expect_identical(metadata(mnn.out)$lost.var, rev(metadata(mnn.out2)$lost.var)) # variance loss at first step should be symmetric.
})

set.seed(1200005)
test_that("fastMNN works as expected for three batches with re-ordering", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, B3, d=50) # corrected values
    expect_identical(dim(out$corrected), c(ncol(B1) + ncol(B2) + ncol(B3), 50L))
    expect_identical(as.integer(out$batch), rep(1:3, c(ncol(B1), ncol(B2), ncol(B3))))
    expect_identical(metadata(out)$order, 1:3)
    CHECK_PAIRINGS(out)
    
    # Testing the re-ordering algorithms.
    out.re <- fastMNN(B3=B3, B2=B2, B1=B1, auto.order=c(3,2,1))
    CHECK_PAIRINGS(out.re)

    back.to.original <- order(out.re$batch)
    expect_equal(out$corrected, out.re$corrected[back.to.original,])
    expect_identical(metadata(out.re)$order, c("B1", "B2", "B3"))

    old.pairs <- metadata(out)$pairs
    new.pairs <- metadata(out.re)$pairs
    for (i in seq_along(out$pairs)) {
        expect_identical(old.pairs[[i]]$first, match(new.pairs[[i]]$first, back.to.original))
        expect_identical(old.pairs[[i]]$second, match(new.pairs[[i]]$second, back.to.original))
    }

    expect_equal(metadata(out)$rotation, metadata(out.re)$rotation)
})

set.seed(1200006)
test_that("fastMNN works as expected for three batches with auto ordering", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000), nrow=100) # Batch 2
    ref <- fastMNN(B1, B2, B3, d=50) 

    # Testing the auto-ordering algorithms. 
    out.auto <- fastMNN(B1, B2, B3, d=50, auto.order=TRUE) 
    expect_identical(dim(ref$corrected), dim(out.auto$corrected))
    expect_identical(out.auto$batch, ref$batch)
    expect_identical(metadata(out.auto)$order, c(2L, 1L, 3L)) # 3 should be last, with the fewest cells => fewest MNNs.
    CHECK_PAIRINGS(out.auto)

    out.re.auto <- fastMNN(B1, B2, B3, auto.order=metadata(out.auto)$order)
    expect_equal(out.auto, out.re.auto)

    # Auto-ordering consistently handles the BiocNeighborIndex class.
    expect_error(out.approx <- fastMNN(B1, B2, B3, auto.order=TRUE, BNPARAM=BiocNeighbors::AnnoyParam()), NA)
    expect_identical(metadata(out.approx)$order, metadata(out.auto)$order)

    # Testing the internal auto-ordering functions.
    fmerge <- batchelor:::.define_first_merge(list(t(B1), t(B2), t(B3)), k=20)   
    expect_identical(fmerge$first, 2L) # 2 is arbitrarily 'first', and 1 is arbitrarily 'second'; but 3 should never show up.
    expect_identical(fmerge$second, 1L)
    expect_identical(fmerge$pairs, findMutualNN(t(B2), t(B1), k1=20, k2=20))

    expect_identical(BiocNeighbors::findKNN(precomputed=fmerge$precomputed[[1]], k=5), BiocNeighbors::findKNN(t(B1), k=5)) # checking that the precomputations are correct.
    expect_identical(BiocNeighbors::findKNN(precomputed=fmerge$precomputed[[2]], k=10), BiocNeighbors::findKNN(t(B2), k=10))
    expect_identical(BiocNeighbors::findKNN(precomputed=fmerge$precomputed[[3]], k=15), BiocNeighbors::findKNN(t(B3), k=15))

    nmerge <- batchelor:::.define_next_merge(t(B1), list(t(B1), t(B2), t(B3)), processed=1L, precomputed=fmerge$precomputed, k=20)   
    expect_identical(nmerge$other, 2L) # 1 should have more MNNs with 2 than with 3.
    expect_identical(nmerge$pairs, findMutualNN(t(B1), t(B2), k1=20, k2=20))
})

set.seed(12000051)
test_that("fastMNN works on SingleCellExperiment inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    ref <- fastMNN(sce1, sce2)
    out <- fastMNN(B1, B2)
    expect_equal(ref, out)

    # Behaves with spikes as input.
    isp <- rbinom(nrow(B1), 1, 0.1)==1L
    isSpike(sce1, "ERCC") <- isp
    isSpike(sce2, "ERCC") <- isp
    ref <- fastMNN(sce1, sce2, d=5)
    out <- fastMNN(B1[!isp,], B2[!isp,], d=5)
    expect_equal(ref, out)

    # Spikes and subsetting interact correctly
    i <- rbinom(nrow(B1), 1, 0.5)==1L
    ref <- fastMNN(sce1, sce2, d=2, subset.row=i)
    out <- fastMNN(sce1[i,], sce2[i,], d=2)
    expect_equal(ref, out)

    # Handles reduced dimensional inputs correctly.
    blah <- multiBatchPCA(sce1, sce2, get.spikes=TRUE)
    reducedDim(sce1) <- blah[[1]]
    reducedDim(sce2) <- blah[[2]]
    out <- fastMNN(sce1, sce2, use.dimred=1)

    ref <- fastMNN(B1, B2, cos.norm=FALSE)
    metadata(ref)$rotation <- NULL
    expect_equal(ref, out)
})

set.seed(12000052)
test_that("fastMNN works with within-object batches", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    sce3 <- SingleCellExperiment(list(logcounts=B3))
    combined <- cbind(sce1, sce2, sce3)
    batches <- rep(1:3, c(ncol(sce1), ncol(sce2), ncol(sce3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    ref <- fastMNN(B1, B2, B3)
    out <- fastMNN(combined, batch=batches)
    expect_equal(ref$corrected[shuffle,], out$corrected)
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # Checking the pairings.
    for (x in seq_along(out$pairs)) {
        curref <- metadata(ref)$pairs[[x]]
        curref[,1] <- match(curref[,1], shuffle)
        curref[,2] <- match(curref[,2], shuffle)
        curout <- out$pairs[[x]]
        expect_identical(
            curref[order(curref[,1], curref[,2]),],
            curout[order(curout[,1], curout[,2]),]
        )
    }
})

set.seed(1200006)
test_that("fastMNN fails on silly inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    # Throws errors properly with no genes or no cells.
    expect_error(fastMNN(), "at least two batches")
    expect_error(fastMNN(B1), "'batch' must be specified")
    expect_error(expect_warning(fastMNN(B1[0,], B2[0,]), "more requested"), "zero")
    expect_error(expect_warning(fastMNN(B1[,0], B2[,0]), "more requested"), "zero")

    # SCE vs matrix errors.    
    expect_error(multiBatchPCA(SingleCellExperiment(list(logcounts=B1)), B2), "cannot mix")

    # Throws errors upon row checks.
    expect_error(fastMNN(B1[1:10,], B2), "number of rows is not the same")
    xB1 <- B1
    xB2 <- B2
    rownames(xB1) <- sample(nrow(B1))
    rownames(xB2) <- sample(nrow(B2))
    expect_error(fastMNN(xB1, xB2), "row names are not the same")

    # Fails if auto.order is not in 1:nbatches.
    expect_error(fastMNN(B1, B2, auto.order=1:3), "permutation of")
})
