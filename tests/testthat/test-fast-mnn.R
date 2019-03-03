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

    # Comparing to the implementation in batchelor.
    out <- batchelor:::.average_correction(test1, mnn1, test2, mnn2)  
    expect_equal(out$averaged, ref)
    expect_identical(out$second, sort(unique(mnn2)))
    expect_identical(out$second, as.integer(names(by.mnn)))

    # Doesn't fail on dummy inputs.
    empty <- batchelor:::.average_correction(test1, integer(0), test2, integer(0))  
    expect_identical(dim(empty$averaged), c(0L, ncol(test2)))
    expect_identical(empty$second, integer(0))
})

set.seed(1200002)
test_that("centering along a batch vector works correctly", {
    test <- matrix(rnorm(1000), ncol=10)
    batch <- rnorm(10) 
    centered <- batchelor:::.center_along_batch_vector(test, batch)
    new.locations <- centered %*% batch
    expect_true(sd(new.locations) < 1e-8)

    # Support restriction.
    test <- matrix(rnorm(1000), ncol=10)
    test2 <- rbind(test, test[1:10,])
    batch <- rnorm(10) 

    original <- batchelor:::.center_along_batch_vector(test, batch)
    keep <- seq_len(nrow(test))
    current <- batchelor:::.center_along_batch_vector(test2, batch, restrict=keep)
    expect_identical(original, current[keep,])
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
    pairings <- metadata(mnn.out)$merge.info$pairs
    order <- metadata(mnn.out)$merge.order

    expect_identical(length(order), length(order))
    expect_true(all(order[1]==origin[pairings[[1]]$first]))

    for (counter in seq_along(pairings)) {
        p <- pairings[[counter]]
        sbatch <- origin[p$second]
        expect_true(all(sbatch == order[counter+1L]))

        fbatch <- origin[p$first]
        expect_true(all(fbatch %in% order[seq_len(counter)]))
    }

    return(NULL)
}

set.seed(1200004)
test_that("fastMNN works as expected for two batches", {
    # Note that all fastMNN checks will have batches that are offset 
    # to ensure that the correction is not skipped via min.batch.effect.
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, d=50) 
    expect_identical(dim(reducedDim(out)), c(ncol(B1) + ncol(B2), 50L))
    expect_identical(as.integer(out$batch), rep(1:2, c(ncol(B1), ncol(B2))))
    expect_identical(metadata(out)$merge.order, 1:2)
    CHECK_PAIRINGS(out)
    
    # Dimension choice behaves correctly.
    out.10 <- fastMNN(B1, B2, d=10) 
    expect_identical(ncol(reducedDim(out.10)), 10L)
    CHECK_PAIRINGS(out.10)

    # Behaves if we turn off cosine-norm.
    nB1 <- t(t(B1)/ sqrt(colSums(B1^2)))
    nB2 <- t(t(B2)/ sqrt(colSums(B2^2)))
    out.ncos <- fastMNN(nB1, nB2, cos.norm=FALSE, d=50) 
    expect_equal(out.ncos, out)

    # Subset.row behaves correctly.
    i <- sample(nrow(B1), 50)
    ref <- fastMNN(X=B1[i,], Y=B2[i,], d=50)
    out.s <- fastMNN(X=B1, Y=B2, d=50, subset.row=i)
    expect_identical(reducedDim(out.s), reducedDim(ref))
    expect_equal(out.s, ref)

    # Behaves if we only use PCs.
    pcs <- multiBatchPCA(B1, B2, d=10)
    out.pre <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
    out.norm <- fastMNN(B1, B2, d=10, cos.norm=FALSE)

    expect_equal(metadata(pcs)$rotation, rowData(out.norm)$rotation)
    expect_equal(out.pre$corrected, reducedDim(out.norm))
    expect_equal(out.pre$batch, out.norm$batch)
})

set.seed(12000041)
test_that("fastMNN handles names correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    out <- fastMNN(X=B1, Y=B2, d=50) 
    expect_identical(out$batch, rep(c("X", "Y"), c(ncol(B1), ncol(B2))))

    rownames(B1) <- rownames(B2) <- sprintf("GENE_%i", seq_len(nrow(B1)))
    out <- fastMNN(B1, B2, d=50) 
    expect_identical(rownames(out), rownames(B1))

    out <- fastMNN(B1, B2, d=50, subset.row=1:50) 
    expect_identical(rownames(out), rownames(B1)[1:50])

    colnames(B1) <- sprintf("CELL_%i", seq_len(ncol(B1)))
    colnames(B2) <- sprintf("CELL_%i", seq_len(ncol(B2)))
    out <- fastMNN(B1, B2, d=50) 
    expect_identical(colnames(out), c(colnames(B1), colnames(B2)))

    # Same for PCs.
    C1 <- matrix(rnorm(1000), ncol=10)
    C2 <- matrix(rnorm(2000), ncol=10)
    out <- fastMNN(Y=C1, X=C2, pc.input=TRUE) 
    expect_identical(out$batch, rep(c("Y", "X"), c(ncol(B1), ncol(B2))))

    rownames(C1) <- sprintf("CELL_%i", seq_len(nrow(C1)))
    rownames(C2) <- sprintf("CELL_%i", seq_len(nrow(C2)))
    out <- fastMNN(C1, C2, pc.input=TRUE)
    expect_identical(rownames(out), c(rownames(C1), rownames(C2)))
    expect_identical(rownames(out), rownames(out$corrected))

    colnames(C1) <- colnames(C2) <- sprintf("PC%i", seq_len(ncol(C1)))
    out <- fastMNN(C1, C2, pc.input=TRUE)
    expect_identical(colnames(out$corrected), colnames(C1))
})

set.seed(12000042)
test_that("variance loss calculations work as expected", {
    PC1 <- matrix(rnorm(10000, 0), ncol=10) # Batch 1 
    PC2 <- matrix(rnorm(20000, 1), ncol=10) # Batch 2

    store <- batchelor:::MNN_supplied_order(list(PC1, PC2))
    store <- batchelor:::.advance(store, k=20)
    out <- batchelor:::.compute_reference_var(PC1, store)
    expect_identical(out, sum(DelayedMatrixStats::colVars(PC1)))

    # Alternative ordering. 
    store <- batchelor:::MNN_supplied_order(list(PC1, PC2), ordering=2:1)
    store <- batchelor:::.advance(store, k=20)
    out2 <- batchelor:::.compute_reference_var(PC2, store)
    expect_identical(out2, sum(DelayedMatrixStats::colVars(PC2)))

    # Multiple lengths.
    PC3 <- matrix(rnorm(5000, 2), ncol=10) # Batch 3

    store <- batchelor:::MNN_supplied_order(list(PC1, PC2, PC3))
    store <- batchelor:::.advance(store, k=20)
    out3a <- batchelor:::.compute_reference_var(PC1, store)
    expect_identical(out3a, out)

    store <- batchelor:::.compile(store, corrected=PC2)
    store <- batchelor:::.advance(store, k=20)
    out3b <- batchelor:::.compute_reference_var(rbind(PC1, PC2), store)
    expect_identical(out3b, c(out, out2))

    # Multiple lengths and a different order.
    store <- batchelor:::MNN_supplied_order(list(PC1, PC2, PC3), ordering=c(2L,3L,1L))
    store <- batchelor:::.advance(store, k=20)
    store <- batchelor:::.compile(store, corrected=PC3)
    store <- batchelor:::.advance(store, k=20)
    out3c <- batchelor:::.compute_reference_var(rbind(PC2, PC3), store)
    expect_identical(out3c, c(out2, sum(DelayedMatrixStats::colVars(PC3))))

    # Checking that we compute something.
    mnn.out <- fastMNN(PC1, PC2, pc.input=TRUE)
    expect_identical(length(metadata(mnn.out)$lost.var), 2L)
    expect_identical(length(unique(metadata(mnn.out)$lost.var)), 2L)

    mnn.out2 <- fastMNN(PC2, PC1, auto.order=2:1, pc.input=TRUE)
    expect_identical(metadata(mnn.out)$lost.var, rev(metadata(mnn.out2)$lost.var)) # variance loss at first step should be symmetric.

    # Variance loss should be zero without a batch effect.
    PC1x <- matrix(rnorm(20000, 0), ncol=10)
    mnn.outx <- fastMNN(PC1, PC1x, pc.input=TRUE, min.batch.skip=TRUE)
    expect_identical(metadata(mnn.outx)$lost.var, numeric(2))
})

set.seed(12000043)
test_that("fastMNN changes the reference batch upon orthogonalization", {
    PC1 <- matrix(rnorm(10000, 0), ncol=10) # Batch 1 
    PC2 <- matrix(rnorm(20000, 1), ncol=10) # Batch 2

    mnn.out <- fastMNN(PC1, PC2, pc.input=TRUE)
    expect_false(isTRUE(all.equal(PC1, mnn.out$corrected[mnn.out$batch==1,])))
    expect_false(isTRUE(all.equal(PC2, mnn.out$corrected[mnn.out$batch==2,])))

    # ... except when there is no batch effect!
    PC1 <- matrix(rnorm(10000, 0), ncol=10)
    PC2 <- matrix(rnorm(20000, 0), ncol=10)

    mnn.out <- fastMNN(PC1, PC2, pc.input=TRUE, min.batch.skip=0.1)
    expect_true(isTRUE(all.equal(PC1, mnn.out$corrected[mnn.out$batch==1,])))
    expect_true(isTRUE(all.equal(PC2, mnn.out$corrected[mnn.out$batch==2,])))
})

set.seed(1200005)
test_that("fastMNN works as expected for three batches with re-ordering", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3

    out <- fastMNN(B1, B2, B3, d=50) 
    expect_identical(dim(reducedDim(out)), c(ncol(B1) + ncol(B2) + ncol(B3), 50L))
    expect_identical(as.integer(out$batch), rep(1:3, c(ncol(B1), ncol(B2), ncol(B3))))
    expect_identical(metadata(out)$merge.order, 1:3)
    CHECK_PAIRINGS(out)

    # Testing the re-ordering algorithms.
    out.re <- fastMNN(B3=B3, B2=B2, B1=B1, auto.order=c(3,2,1))
    CHECK_PAIRINGS(out.re)

    back.to.original <- order(out.re$batch)
    expect_equal(reducedDim(out), reducedDim(out.re)[back.to.original,])
    expect_identical(metadata(out.re)$merge.order, c("B1", "B2", "B3"))

    old.pairs <- metadata(out)$merge.info$pairs
    new.pairs <- metadata(out.re)$merge.info$pairs
    for (i in seq_along(old.pairs)) {
        expect_identical(old.pairs[[i]]$first, match(new.pairs[[i]]$first, back.to.original))
        expect_identical(old.pairs[[i]]$second, match(new.pairs[[i]]$second, back.to.original))
    }

    expect_equal(metadata(out)$rotation, metadata(out.re)$rotation)
})

set.seed(12000050)
test_that("fastMNN works as expected for three batches with auto ordering", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3
    ref <- fastMNN(B1, B2, B3, d=50) 

    # Testing the auto-ordering algorithms. 
    out.auto <- fastMNN(B1, B2, B3, d=50, auto.order=TRUE) 
    expect_identical(dim(reducedDim(ref)), dim(reducedDim(out.auto)))
    expect_identical(out.auto$batch, ref$batch)
    expect_identical(metadata(out.auto)$merge.order, c(2L, 1L, 3L)) # 3 should be last, with the fewest cells => fewest MNNs.
    CHECK_PAIRINGS(out.auto)

    out.re.auto <- fastMNN(B1, B2, B3, auto.order=metadata(out.auto)$merge.order)
    expect_equal(out.auto, out.re.auto)

    # Auto-ordering consistently handles the BiocNeighborIndex class.
    expect_error(out.approx <- fastMNN(B1, B2, B3, auto.order=TRUE, BNPARAM=BiocNeighbors::AnnoyParam()), NA)
    expect_identical(metadata(out.approx)$merge.order, metadata(out.auto)$merge.order)
})

set.seed(120000500)
test_that("fastMNN computes the batch size correctly and skips no-batch scenarios", {
    B1 <- matrix(rnorm(50000), nrow=100) 
    B2x <- matrix(rnorm(50000, 1), nrow=100) 
    B2y <- matrix(rnorm(50000), nrow=100) 

    out <- fastMNN(B1, B2x, d=50)
    expect_true(all(metadata(out)$merge.info$batch.size > 0.5))
    expect_false(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var > 0))

    out <- fastMNN(B1, B2y, d=50)
    expect_true(all(metadata(out)$merge.info$batch.size < 0.1))
    expect_false(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var > 0))
    CHECK_PAIRINGS(out)

    out <- fastMNN(B1, B2y, d=50, min.batch.skip=0.1)
    expect_true(all(metadata(out)$merge.info$batch.size < 0.1))
    expect_true(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var == 0))

    ref <- multiBatchPCA(cosineNorm(B1), cosineNorm(B2y))
    expect_identical(reducedDim(out), rbind(ref[[1]], ref[[2]]))
    CHECK_PAIRINGS(out)
})

set.seed(120000501)
test_that("fastMNN with three batches behaves in the absence of a batch effect", {
    # Just checking that min.batch.effect doesn't do weird things with auto.order.
    B1x <- matrix(rnorm(100000), nrow=100) # Batch 1 
    B2x <- matrix(rnorm(200000), nrow=100) # Batch 2
    B3x <- matrix(rnorm(50000), nrow=100) # Batch 3

    out3 <- fastMNN(B1x, B2x, B3x, d=50, min.batch.skip=0.1)
    expect_identical(metadata(out3)$merge.order, 1:3)
    expect_identical(metadata(out3)$lost.var, numeric(3))
    expect_identical(metadata(out3)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3)

    out3r <- fastMNN(B1x, B2x, B3x, d=50, auto.order=3:1, min.batch.skip=0.1)
    expect_identical(metadata(out3r)$merge.order, 3:1)
    expect_identical(metadata(out3r)$lost.var, numeric(3))
    expect_identical(metadata(out3r)$lost.var, numeric(3))
    expect_identical(metadata(out3r)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3r)

    out3a <- fastMNN(B1x, B2x, B3x, d=50, auto.order=TRUE, min.batch.skip=0.1)
    expect_identical(metadata(out3a)$lost.var, numeric(3))
    expect_identical(metadata(out3a)$merge.order, c(2L, 1L, 3L)) 
    expect_identical(metadata(out3a)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3a)
})

set.seed(120000502)
test_that("Orthogonalization is performed correctly across objects", {
    B1 <- matrix(rnorm(10000), ncol=10)
    B2 <- matrix(rnorm(20000, 1), ncol=10)
    B3 <- matrix(rnorm(5000, rep(0:1, each=5)), ncol=10) 
    B4 <- matrix(rnorm(8000, rep(0:1, 5)), ncol=10) 

    ref <- fastMNN(B1, B2, B3, B4, pc.input=TRUE)

    outA <- fastMNN(B1, B2, pc.input=TRUE)
    outB <- fastMNN(outA, B3)
    outC <- fastMNN(outB, B4)
    expect_equal(ref$corrected, outC$corrected, tol=0.001) # not theoretically equal, but should be close.
    expect_identical(nrow(metadata(outC)$orthogonalize), 2L)

    # Handles hierarhical merges.
    out.i <- fastMNN(B1, B2, pc.input=TRUE)
    out.ii <- fastMNN(B3, B4, pc.input=TRUE)
    out.iii <- fastMNN(out.i, out.ii)
    expect_identical(nrow(metadata(out.iii)$orthogonalize), 2L)

    # Handles restriction correctly.
    B1x <- rbind(B1, B1[1:10,])
    B2x <- rbind(B2, B2[1:10,])
    B3x <- rbind(B3, B3[1:10,])
    B4x <- rbind(B4, B4[1:10,])

    collected <- list(seq_len(nrow(B1)), seq_len(nrow(B2)))
    outAx <- fastMNN(B1x, B2x, restrict=collected, pc.input=TRUE)
    collected[[2]] <- collected[[2]] + nrow(B1x)
    collected[[3]] <- seq_len(nrow(B3))
    outBx <- fastMNN(outAx, B3x, restrict=list(unlist(collected[1:2]), collected[[3]]))
    collected[[3]] <- collected[[3]] + nrow(outAx)
    collected[[4]] <- seq_len(nrow(B4))
    outCx <- fastMNN(outBx, B4x, restrict=list(unlist(collected[1:3]), collected[[4]]))
    collected[[4]] <- collected[[4]] + nrow(outBx)
    
    expect_equal(outCx$corrected[unlist(collected),], outC$corrected)
    first.ten <- lapply(collected, head, 10)
    expect_equal(outCx$corrected[unlist(first.ten),], outCx$corrected[-unlist(collected),])

    # Handles skipping correctly.
    outAs <- fastMNN(B1, B1, pc.input=TRUE, min.batch.skip=0.1)
    outBs <- fastMNN(outAs, B2)
    expect_identical(nrow(metadata(outBs)$orthogonalize), 0L)

    outAs <- fastMNN(B1, B1, pc.input=TRUE, min.batch.skip=0)
    outBs <- fastMNN(outAs, B2)
    expect_identical(nrow(metadata(outBs)$orthogonalize), 1L)
})

set.seed(120000503)
test_that("orthogonalization is done correctly in exact toy examples", {
    core <- cbind(rep(1:10, each=10), rep(1:10, 10))
    b1 <- b2 <- core
    b1[,1] <- b1[,1] + 20
    b2[,2] <- b2[,2] + 20

    # Orthogonalized on x-axis.
    out1 <- fastMNN(core, b1, k=1, pc.input=TRUE)
    expect_equal(out1$corrected[,1], rep(5.5, nrow(out1)))
    expect_equal(out1$corrected[,2], c(core[,2], b1[,2]))

    # Orthogonalized on y-axis.
    out2 <- fastMNN(core, b1, b2, k=1, pc.input=TRUE)
    expect_equal(out2$corrected[,1], rep(5.5, nrow(out2)))
    expect_equal(out2$corrected[,2], rep(5.5, nrow(out2)))

    # Orthogonalized step by step.
    out3 <- fastMNN(out1, b2, k=1)
    expect_equal(out3$corrected[,1], rep(5.5, nrow(out3)))
    expect_equal(out3$corrected[,2], rep(5.5, nrow(out3)))

    # Uses the information in 'orthogonalize'.
    out4 <- fastMNN(out3, core + 10, k=1)
    expect_equal(out4$corrected[,1], rep(5.5, nrow(out4)))
    expect_equal(out4$corrected[,2], rep(5.5, nrow(out4)))

    # Hierarchical orthogonalization works.
    outX <- fastMNN(core, b1, pc.input=TRUE, k=1)
    outY <- fastMNN(core+10, b2+10, pc.input=TRUE, k=1)
    expect_equal(outY$corrected[,1], c(core[,1], b2[,1])+10)
    expect_equal(outY$corrected[,2], rep(15.5, nrow(outY)))

    outZ <- fastMNN(outX, outY, k=1)
    expect_equal(outZ$corrected[,1], rep(5.5, nrow(out4)))
    expect_equal(outZ$corrected[,2], rep(5.5, nrow(out4)))
})

set.seed(12000051)
test_that("fastMNN works on SingleCellExperiment inputs", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    ref <- fastMNN(sce1, sce2)
    out <- fastMNN(B1, B2)
    expect_equal(ref, out)

    # Behaves with spikes as input.
    isp <- rbinom(nrow(B1), 1, 0.1)==1L
    isSpike(sce1, "ERCC") <- isp
    isSpike(sce2, "ERCC") <- isp
    out <- fastMNN(sce1, sce2, d=5)
    ref <- fastMNN(B1[!isp,], B2[!isp,], d=5)
    expect_equal(ref, out)

    # Spikes and subsetting interact correctly
    i <- rbinom(nrow(B1), 1, 0.5)==1L
    out <- fastMNN(sce1, sce2, d=2, subset.row=i)
    ref <- fastMNN(sce1[i,], sce2[i,], d=2)
    expect_equal(ref, out)
    ref2 <- fastMNN(B1[!isp & i,], B2[!isp & i,], d=2)
    expect_equal(ref2, out)

    # Handles reduced dimensional inputs correctly.
    blah <- multiBatchPCA(sce1, sce2, get.spikes=TRUE)
    reducedDim(sce1) <- blah[[1]]
    reducedDim(sce2) <- blah[[2]]
    out <- fastMNN(sce1, sce2, use.dimred=1)

    ref <- fastMNN(B1, B2, cos.norm=FALSE)
    expect_identical(reducedDim(ref), out$corrected)
    expect_identical(ref$batch, out$batch)
    expect_identical(metadata(ref), metadata(out))
})

set.seed(12000052)
test_that("fastMNN works with within-object batches", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 2), nrow=100) # Batch 2

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
    expect_equal(reducedDim(ref)[shuffle,], reducedDim(out))
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))
    CHECK_PAIRINGS(out)

    # Checking the pairings in closer detail.
    pairings <- metadata(out)$pairs
    ref.pairings <- metadata(ref)$pairs
    for (x in seq_along(pairings)) {
        curref <- ref.pairings[[x]]
        curref[,1] <- match(curref[,1], shuffle)
        curref[,2] <- match(curref[,2], shuffle)
        curout <- pairings[[x]]
        expect_identical(
            curref[order(curref[,1], curref[,2]),],
            curout[order(curout[,1], curout[,2]),]
        )
    }
})

set.seed(120000521)
test_that("fastMNN works with within-object batches for PCs", {
    pcd <- list(
        matrix(rnorm(10000, 0), ncol=50), # Batch 1 
        matrix(rnorm(20000, 1), ncol=50), # Batch 2
        matrix(rnorm(15000, 2), ncol=50) # Batch 2
    )
    com.pcd <- do.call(rbind, pcd)
    batches <- rep(1:3, vapply(pcd, nrow, FUN.VALUE=0L))

    shuffle <- sample(nrow(com.pcd))
    com.pcd <- com.pcd[shuffle,]
    batches <- batches[shuffle]

    # Splits PC inputs properly.
    out <- fastMNN(com.pcd, batch=batches, pc.input=TRUE)
    ref <- fastMNN(pcd[[1]], pcd[[2]], pcd[[3]], pc.input=TRUE)
    expect_equal(ref$corrected[shuffle,], out$corrected)
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # ... and same for use.dimred=.
    sce1 <- SingleCellExperiment(matrix(0, 0, nrow(pcd[[1]])), reducedDims=list(PCA=pcd[[1]]))
    sce2 <- SingleCellExperiment(matrix(0, 0, nrow(pcd[[2]])), reducedDims=list(PCA=pcd[[2]]))
    sce3 <- SingleCellExperiment(matrix(0, 0, nrow(pcd[[3]])), reducedDims=list(PCA=pcd[[3]]))
    combined <- SingleCellExperiment(matrix(0, 0, nrow(com.pcd)), reducedDims=list(PCA=com.pcd)) 

    out <- fastMNN(combined, batch=batches, use.dimred="PCA")
    ref <- fastMNN(sce1, sce2, sce3, use.dimred="PCA")
    expect_equal(ref$corrected[shuffle,], out$corrected)
    expect_equal(as.character(ref$batch[shuffle]), out$batch)
})

set.seed(120000522)
test_that("fastMNN renames within-object batches correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 2), nrow=100) # Batch 2

    rownames(B1) <- rownames(B2) <- rownames(B3) <- sprintf("GENE_%i", sample(nrow(B1)))
    colnames(B1) <- sprintf("CELL_1_%i", seq_len(ncol(B1)))
    colnames(B2) <- sprintf("CELL_2_%i", seq_len(ncol(B2)))
    colnames(B3) <- sprintf("CELL_3_%i", seq_len(ncol(B3)))

    combined <- cbind(B1, B2, B3)
    batches <- rep(1:3, c(ncol(B1), ncol(B2), ncol(B3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    out <- fastMNN(combined, batch=batches)
    expect_identical(rownames(out), rownames(B1))
    expect_identical(colnames(out), colnames(combined))
})

set.seed(12000053)
test_that("fastMNN works correctly with restriction", {
    B1 <- matrix(rnorm(10000, 0), ncol=10) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), ncol=10) # Batch 2
    B3 <- matrix(rnorm(5000), ncol=10) # Batch 3

    # Restricted results are only directly comparable if we're not doing a PCA internally.
    ref <- fastMNN(B1, B2, B3, pc.input=TRUE) 

    i1 <- 100:50
    C1 <- rbind(B1, B1[i1,])
    i2 <- 1:20
    C2 <- rbind(B2, B2[i2,])
    i3 <- 50:100
    C3 <- rbind(B3, B3[i3,])

    keep1 <- seq_len(nrow(B1))
    keep2 <- seq_len(nrow(B2))
    keep3 <- seq_len(nrow(B3))
    out <- fastMNN(C1, C2, C3, pc.input=TRUE, restrict=list(keep1, keep2, keep3)) 
    CHECK_PAIRINGS(out)

    expect_identical(ref$corrected[ref$batch==1,], out$corrected[out$batch==1,][keep1,])
    expect_identical(ref$corrected[ref$batch==2,], out$corrected[out$batch==2,][keep2,])
    expect_identical(ref$corrected[ref$batch==3,], out$corrected[out$batch==3,][keep3,])

    expect_identical(ref$corrected[ref$batch==1,][i1,], out$corrected[out$batch==1,][-keep1,])
    expect_identical(ref$corrected[ref$batch==2,][i2,], out$corrected[out$batch==2,][-keep2,])
    expect_identical(ref$corrected[ref$batch==3,][i3,], out$corrected[out$batch==3,][-keep3,])

    # Also behaves with a single batch.
    DY <- rbind(C1, C2, C3)
    batch <- rep(1:3, c(nrow(C1), nrow(C2), nrow(C3)))
    shuffle <- sample(length(batch))

    out2 <- fastMNN(DY[shuffle,], batch=batch[shuffle], pc.input=TRUE, 
        restrict=list(shuffle %in% c(keep1, keep2 + nrow(C1), keep3 + nrow(C1) + nrow(C2))))
    CHECK_PAIRINGS(out2)
    expect_equal(out2$corrected, out$corrected[shuffle,])
    expect_identical(out2$batch, as.character(out$batch)[shuffle])
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
    expect_error(fastMNN(SingleCellExperiment(list(logcounts=B1)), B2), "cannot mix")

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
