# This tests the functions related to fastMNN.
# library(batchelor); library(testthat); source("test-fast-mnn.R")

library(BiocSingular)

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
    origin <- mnn.out$batch
    merge.info <- metadata(mnn.out)$merge.info

    expect_identical(length(unique(origin)) - 1L, nrow(merge.info))
    expect_identical(sort(unique(origin)), sort(unique(c(unlist(merge.info$left), unlist(merge.info$right)))))

    for (counter in seq_len(nrow(merge.info))) {
        left <- merge.info$left[[counter]]
        right <- merge.info$right[[counter]]
        expect_identical(length(intersect(left, right)), 0L)

        p <- merge.info$pairs[[counter]]

        allowed.left <- which(origin %in% left)
        expect_true(length(p$left) > 0)
        expect_true(all(p$left %in% allowed.left))
        
        allowed.right <- which(origin %in% right)
        expect_true(length(p$right) > 0)
        expect_true(all(p$right %in% allowed.right))
    }

    invisible(NULL)
}

set.seed(1200004)
test_that("fastMNN works as expected for two batches", {
    # Note that all fastMNN checks will have batches that are offset 
    # to ensure that the correction is not skipped via min.batch.effect.
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, d=50, BSPARAM=ExactParam()) 
    expect_identical(dim(reducedDim(out)), c(ncol(B1) + ncol(B2), 50L))
    expect_identical(as.integer(out$batch), rep(1:2, c(ncol(B1), ncol(B2))))
    expect_identical(metadata(out)$merge.info$left[[1]], 1L)
    expect_identical(metadata(out)$merge.info$right[[1]], 2L)
    CHECK_PAIRINGS(out)
    
    # Dimension choice behaves correctly.
    out.10 <- fastMNN(B1, B2, d=10, BSPARAM=ExactParam())
    expect_identical(ncol(reducedDim(out.10)), 10L)
    CHECK_PAIRINGS(out.10)

    # Behaves if we turn off cosine-norm.
    nB1 <- t(t(B1)/ sqrt(colSums(B1^2)))
    nB2 <- t(t(B2)/ sqrt(colSums(B2^2)))
    out.ncos <- fastMNN(nB1, nB2, cos.norm=FALSE, d=50, BSPARAM=ExactParam())
    expect_equal(out.ncos, out)
})

set.seed(120000401)
test_that("fastMNN subsets correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    # Subset.row behaves correctly.
    i <- sample(nrow(B1), 50)
    ref <- fastMNN(X=B1[i,], Y=B2[i,], d=50, BSPARAM=ExactParam())
    out.s <- fastMNN(X=B1, Y=B2, d=50, subset.row=i, BSPARAM=ExactParam())
    expect_identical(reducedDim(out.s), reducedDim(ref))
    expect_equal(out.s, ref)

    # Correct.all behaves correctly.
    i <- c(1:nrow(B1), 1:10)
    ref <- fastMNN(X=B1, Y=B2, d=50, BSPARAM=ExactParam())
    out <- fastMNN(X=B1[i,], Y=B2[i,], subset.row=1:nrow(B1), d=50, BSPARAM=ExactParam(), correct.all=TRUE)
    expect_identical(nrow(out), length(i))
    expect_identical(reducedDim(ref), reducedDim(out))
    expect_equal(as.matrix(assay(ref)[1:10,]), as.matrix(assay(out)[1:10+nrow(B1),]))
})

set.seed(12000041)
test_that("fastMNN handles names correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    out <- fastMNN(X=B1, Y=B2, BSPARAM=ExactParam()) 
    expect_identical(out$batch, rep(c("X", "Y"), c(ncol(B1), ncol(B2))))
    expect_identical(metadata(out)$merge.info$left[[1]], "X")
    expect_identical(metadata(out)$merge.info$right[[1]], "Y")

    # Handles row names.
    rownames(B1) <- rownames(B2) <- sprintf("GENE_%i", seq_len(nrow(B1)))
    out <- fastMNN(B1, B2, BSPARAM=ExactParam())
    expect_identical(rownames(out), rownames(B1))

    out <- fastMNN(B1, B2, subset.row=1:50, BSPARAM=ExactParam()) 
    expect_identical(rownames(out), rownames(B1)[1:50])

    out <- fastMNN(B1, B2, subset.row=1:50, BSPARAM=ExactParam(), correct.all=TRUE) 
    expect_identical(rownames(out), rownames(B1))

    # Handles column names.
    colnames(B1) <- sprintf("CELL_%i", seq_len(ncol(B1)))
    colnames(B2) <- sprintf("CELL_%i", seq_len(ncol(B2)))
    out <- fastMNN(B1, B2, d=50, BSPARAM=ExactParam())
    expect_identical(colnames(out), c(colnames(B1), colnames(B2)))
})

set.seed(12000042)
test_that("variance loss calculations work as expected", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    out <- fastMNN(B1, B2, BSPARAM=ExactParam()) 
    expect_true(all(metadata(out)$merge.info$lost.var > 0))
    expect_identical(length(unique(metadata(out)$merge.info$lost.var)), 2L)

    # Multiple lengths.
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3
    out2 <- fastMNN(B1, B2, B3, BSPARAM=ExactParam()) 
    expect_true(all(metadata(out2)$merge.info$lost.var[1,-3] > 0))
    expect_equal(metadata(out2)$merge.info$lost.var[1,3], 0)
    expect_true(all(metadata(out2)$merge.info$lost.var[2,] > 0))
    expect_identical(length(unique(metadata(out2)$merge.info$lost.var)), 6L)

    # Respects names.
    out3 <- fastMNN(X=B1, Y=B2, Z=B3, BSPARAM=ExactParam()) 
    expect_identical(colnames(metadata(out3)$merge.info$lost.var), c("X", "Y", "Z"))

    # Variance loss should be zero without a batch effect.
    mnn.outx <- fastMNN(B1, B1, min.batch.skip=TRUE, BSPARAM=ExactParam()) 
    expect_identical(metadata(mnn.outx)$merge.info$lost.var, matrix(0, 1, 2))

    mnn.outx <- fastMNN(B1, B1, B1, min.batch.skip=TRUE, BSPARAM=ExactParam()) 
    expect_identical(metadata(mnn.outx)$merge.info$lost.var, matrix(0, 2, 3))
})

set.seed(120000421)
test_that("fastMNN uses 'prop.k' correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(10000, 1), nrow=100) # Batch 2

    ref <- fastMNN(B1, B2, BSPARAM=ExactParam()) 
    out <- fastMNN(B1, B2, BSPARAM=ExactParam(), k=10, prop.k=20/ncol(B1)) 
    expect_identical(reducedDim(ref), reducedDim(out))

    # max() kicks in.
    ref <- fastMNN(B1, B2, BSPARAM=ExactParam()) 
    out <- fastMNN(B1, B2, BSPARAM=ExactParam(), prop.k=0)
    expect_identical(reducedDim(ref), reducedDim(out))

    # Actually causes a difference in results.
    B2a <- matrix(rnorm(20000, 1), nrow=100)
    ref <- fastMNN(B1, B2a, BSPARAM=ExactParam()) 
    out <- fastMNN(B1, B2a, BSPARAM=ExactParam(), prop.k=20/ncol(B1))
    expect_false(identical(reducedDim(ref), reducedDim(out)))
})

set.seed(12000043)
test_that("fastMNN changes the reference batch upon orthogonalization", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    PCA <- multiBatchPCA(B1, B2, BSPARAM=ExactParam())

    mnn.out <- fastMNN(B1, B2, BSPARAM=ExactParam(), cos.norm=FALSE)
    expect_false(isTRUE(all.equal(PCA[[1]], reducedDim(mnn.out)[mnn.out$batch==1,])))
    expect_false(isTRUE(all.equal(PCA[[2]], reducedDim(mnn.out)[mnn.out$batch==2,])))

    # ... except when there is no batch effect!
    B1 <- matrix(rnorm(10000, 0), nrow=100)
    B2 <- matrix(rnorm(20000, 0), nrow=100)
    PCA <- multiBatchPCA(B1, B2, BSPARAM=ExactParam())

    mnn.out <- fastMNN(B1, B2, min.batch.skip=0.2, BSPARAM=ExactParam(), cos.norm=FALSE) 
    expect_true(isTRUE(all.equal(PCA[[1]], reducedDim(mnn.out)[mnn.out$batch==1,])))
    expect_true(isTRUE(all.equal(PCA[[2]], reducedDim(mnn.out)[mnn.out$batch==2,])))
})

set.seed(1200005)
test_that("fastMNN works as expected for three batches with re-ordering", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3

    out <- fastMNN(B1, B2, B3, BSPARAM=ExactParam())
    expect_identical(dim(reducedDim(out)), c(ncol(B1) + ncol(B2) + ncol(B3), 50L))
    expect_identical(as.integer(out$batch), rep(1:3, c(ncol(B1), ncol(B2), ncol(B3))))
    CHECK_PAIRINGS(out)

    expect_identical(metadata(out)$merge.info$left[[1]], 1L)
    expect_identical(metadata(out)$merge.info$right[[1]], 2L)
    expect_identical(metadata(out)$merge.info$left[[2]], 1:2)
    expect_identical(metadata(out)$merge.info$right[[2]], 3L)

    # Testing the re-ordering algorithms.
    out.re <- fastMNN(B3=B3, B2=B2, B1=B1, merge.order=c(3,2,1), BSPARAM=ExactParam())
    CHECK_PAIRINGS(out.re)

    back.to.original <- order(out.re$batch)
    expect_equal(reducedDim(out), reducedDim(out.re)[back.to.original,])
    expect_identical(c("B1", "B2", "B3")[out$batch], out.re$batch[back.to.original])

    # Checking for coherent diagnostic data.
    expect_identical(metadata(out.re)$merge.info$left[[1]], "B1")
    expect_identical(metadata(out.re)$merge.info$right[[1]], "B2")
    expect_identical(metadata(out.re)$merge.info$left[[2]], c("B1", "B2"))
    expect_identical(metadata(out.re)$merge.info$right[[2]], "B3")

    old.pairs <- metadata(out)$merge.info$pairs
    new.pairs <- metadata(out.re)$merge.info$pairs
    for (i in seq_along(old.pairs)) {
        expect_identical(old.pairs[[i]]$left, match(new.pairs[[i]]$left, back.to.original))
        expect_identical(old.pairs[[i]]$right, match(new.pairs[[i]]$right, back.to.original))
    }

    expect_equivalent(metadata(out)$merge.info$lost.var, 
        metadata(out.re)$merge.info$lost.var[,c("B1", "B2", "B3")])

    # Works when specified by name.
    out.re2 <- fastMNN(B3=B3, B2=B2, B1=B1, merge.order=c("B1","B2","B3"), BSPARAM=ExactParam())
    expect_identical(reducedDim(out.re), reducedDim(out.re2))
})

set.seed(12000050)
test_that("fastMNN works as expected for three batches with auto ordering", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3
    ref <- fastMNN(B1, B2, B3, BSPARAM=ExactParam())

    # Testing the auto-ordering algorithms. 
    out.auto <- fastMNN(B1, B2, B3, auto.merge=TRUE, BSPARAM=ExactParam())
    expect_identical(dim(reducedDim(ref)), dim(reducedDim(out.auto)))
    expect_identical(out.auto$batch, ref$batch)
    CHECK_PAIRINGS(out.auto)

    # Inspecting the merge order.
    # 3 should be last, with the fewest cells => fewest MNNs.
    expect_identical(metadata(out.auto)$merge.info$left[[1]], 2L)
    expect_identical(metadata(out.auto)$merge.info$right[[1]], 1L)
    expect_identical(metadata(out.auto)$merge.info$left[[2]], c(2L, 1L))
    expect_identical(metadata(out.auto)$merge.info$right[[2]], 3L)

    # Checking out that I get the same result if I pass in the merge order exactly.
    out.re.auto <- fastMNN(B1, B2, B3, merge.order=c(2L, 1L, 3L), BSPARAM=ExactParam())
    expect_equal(out.auto, out.re.auto)
})

set.seed(12000050)
test_that("fastMNN works as expected for four batches with trees", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3
    B4 <- matrix(rnorm(3000, 3), nrow=100) # Batch 4

    # Checking consistency with linear merges.
    out.tree <- fastMNN(B1, B2, B3, B4, merge.order=list(list(list(1, 2), 3), 4), BSPARAM=ExactParam())
    out.default <- fastMNN(B1, B2, B3, B4, BSPARAM=ExactParam())
    expect_equal(out.default, out.tree)

    out.tree2 <- fastMNN(B1, B2, B3, B4, merge.order=list(list(list(4, 3), 2), 1), BSPARAM=ExactParam())
    out.reorder <- fastMNN(B1, B2, B3, B4, merge.order=4:1, BSPARAM=ExactParam())
    expect_equal(out.reorder, out.tree2)

    # Trying an actual hierarchical merge.
    out.treeX <- fastMNN(B1, B2, B3, B4, merge.order=list(list(4, 3), list(1, 2)), BSPARAM=ExactParam())

    CHECK_PAIRINGS(out.treeX)
    expect_identical(dim(reducedDim(out.treeX)), dim(reducedDim(out.tree)))
    expect_identical(out.treeX$batch, out.tree$batch)
    
    expect_false(isTRUE(all.equal(reducedDim(out.treeX), reducedDim(out.tree))))
    expect_false(isTRUE(all.equal(reducedDim(out.treeX), reducedDim(out.tree2))))

    # Checking that the each of the subtrees were merged in isolation,
    # and thus have the same diagnostics as if they were merged at the start of a linear order.
    expect_equal(metadata(out.reorder)$merge.info[1,], metadata(out.treeX)$merge.info[2,])
    expect_equal(metadata(out.default)$merge.info[1,], metadata(out.treeX)$merge.info[1,])

    # Works for strings and factors as well.
    out.treeY <- fastMNN(X=B1, A=B2, Z=B3, B=B4, merge.order=list(list("B", "Z"), list("X", "A")), BSPARAM=ExactParam())
    expect_true(isTRUE(all.equal(reducedDim(out.treeX), reducedDim(out.treeY))))

    f <- factor(paste0(letters[1:4], "_"))
    out.treeZ <- fastMNN(a_=B1, b_=B2, c_=B3, d_=B4, merge.order=list(list(f[4], f[3]), list(f[1], f[2])), BSPARAM=ExactParam())
    expect_identical(reducedDim(out.treeX), reducedDim(out.treeZ))
})

set.seed(120000051)
test_that("fastMNN works as expected for many batches with auto ordering", {
    # Creating several batches with different composition.
    collected <- list()
    for (i in 1:10) {
        n <- round(runif(1, 5, 20))*10
        stuff <- runif(100)
        collected[[i]] <- matrix(rnorm(n*length(stuff), stuff), nrow=length(stuff)) 
    }

    ref <- do.call(fastMNN, c(collected, list(BSPARAM=ExactParam(), auto.merge=TRUE)))

    # Checking that all our batches are present in the output.
    by.batch <- rle(ref$batch)
    expect_identical(by.batch$value, seq_along(collected))
    expect_identical(by.batch$length, vapply(collected, ncol, 0L)) 

    # The algorithm should yield the same results after reordering if it is correct.
    # Note that the reordering is done in a manner that preserves the chosen reference at each step
    # (the later batch is chosen as the reference, but the reference is reported first, hence the rev()).
    # Otherwise the results will not be identical due to the asymmetry of each merge step.
    last <- metadata(ref)$merge.info[length(collected)-1,]
    s <- rev(c(last$left[[1]], last$right[[1]]))
    expect_identical(sort(s), seq_along(s))

    alt <- do.call(fastMNN, c(collected[s], list(BSPARAM=ExactParam(), auto.merge=TRUE)))
    expect_identical(metadata(alt)$merge.info$left[[length(collected)-1L]], 10:2)

    o <- order(s[alt$batch])
    expect_equal(reducedDim(ref), reducedDim(alt)[o,])
})

set.seed(120000500)
test_that("fastMNN computes the batch size correctly and skips no-batch scenarios", {
    B1 <- matrix(rnorm(50000), nrow=100) 
    B2x <- matrix(rnorm(50000, 1), nrow=100) 
    B2y <- matrix(rnorm(50000), nrow=100) 

    out <- fastMNN(B1, B2x, d=50, BSPARAM=ExactParam())
    expect_true(all(metadata(out)$merge.info$batch.size > 0.5))
    expect_false(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var > 0))

    out <- fastMNN(B1, B2y, d=50, BSPARAM=ExactParam())
    expect_true(all(metadata(out)$merge.info$batch.size < 0.1))
    expect_false(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var > 0))
    CHECK_PAIRINGS(out)

    out <- fastMNN(B1, B2y, d=50, min.batch.skip=0.1, BSPARAM=ExactParam())
    expect_true(all(metadata(out)$merge.info$batch.size < 0.1))
    expect_true(metadata(out)$merge.info$skipped)
    expect_true(all(metadata(out)$merge.info$lost.var == 0))

    ref <- multiBatchPCA(cosineNorm(B1), cosineNorm(B2y), BSPARAM=ExactParam())
    expect_identical(reducedDim(out), rbind(ref[[1]], ref[[2]]))
    CHECK_PAIRINGS(out)
})

set.seed(120000501)
test_that("fastMNN with three batches behaves in the absence of a batch effect", {
    # Just checking that min.batch.effect doesn't do weird things with auto.order.
    B1x <- matrix(rnorm(100000), nrow=100) # Batch 1 
    B2x <- matrix(rnorm(200000), nrow=100) # Batch 2
    B3x <- matrix(rnorm(50000), nrow=100) # Batch 3

    out3 <- fastMNN(B1x, B2x, B3x, d=50, min.batch.skip=0.1, BSPARAM=ExactParam())
    expect_identical(metadata(out3)$merge.info$lost.var, matrix(0, 2, 3))
    expect_identical(metadata(out3)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3)

    out3r <- fastMNN(B1x, B2x, B3x, d=50, merge.order=3:1, min.batch.skip=0.1, BSPARAM=ExactParam())
    expect_identical(metadata(out3r)$merge.info$lost.var, matrix(0, 2, 3))
    expect_identical(metadata(out3r)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3r)

    out3a <- fastMNN(B1x, B2x, B3x, d=50, auto.merge=TRUE, min.batch.skip=0.1, BSPARAM=ExactParam())
    expect_identical(metadata(out3a)$merge.info$lost.var, matrix(0, 2, 3))
    expect_identical(metadata(out3a)$merge.info$skipped, !logical(2))
    CHECK_PAIRINGS(out3a)
})

set.seed(12000051)
test_that("fastMNN works on SingleCellExperiment inputs", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    ref <- fastMNN(sce1, sce2, BSPARAM=ExactParam())
    out <- fastMNN(B1, B2, BSPARAM=ExactParam())
    expect_equal(ref, out)
})

set.seed(12000052)
test_that("fastMNN works with within-object batches", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 2), nrow=100) # Batch 3

    sce1 <- SingleCellExperiment(list(logcounts=B1))
    sce2 <- SingleCellExperiment(list(logcounts=B2))
    sce3 <- SingleCellExperiment(list(logcounts=B3))
    combined <- cbind(sce1, sce2, sce3)
    batches <- rep(1:3, c(ncol(sce1), ncol(sce2), ncol(sce3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    ref <- fastMNN(B1, B2, B3, BSPARAM=ExactParam())
    out <- fastMNN(combined, batch=batches, BSPARAM=ExactParam())
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

    # Works for character vectors and factors as well.
    alt <- fastMNN(combined, batch=paste0("BATCH", batches), BSPARAM=ExactParam())
    expect_identical(reducedDim(out), reducedDim(alt))
    expect_identical(alt$batch, paste0("BATCH", batches))

    alt <- fastMNN(combined, batch=factor(batches), BSPARAM=ExactParam())
    expect_identical(reducedDim(out), reducedDim(alt))
    expect_identical(alt$batch, as.character(batches))

    alt <- fastMNN(combined, batch=batches, merge.order=list(factor(1), factor(2), factor(3)), BSPARAM=ExactParam())
    expect_identical(reducedDim(out), reducedDim(alt))

    # Handles non-obvious level orderings.
    f <- factor(4 - batches, 3:1) 
    alt <- fastMNN(combined, batch=f, BSPARAM=ExactParam())
    expect_identical(reducedDim(out), reducedDim(alt))

    f <- factor(batches, 3:1)
    alt <- fastMNN(combined, batch=f, BSPARAM=ExactParam())
    rev <- fastMNN(combined, batch=batches, merge.order=3:1, BSPARAM=ExactParam())
    expect_equal(reducedDim(rev), reducedDim(alt))
    expect_identical(alt$batch, as.character(batches))

    f <- factor(batches, 3:1)
    alt <- fastMNN(combined, batch=f, merge.order=factor(1:3), BSPARAM=ExactParam())
    expect_equal(reducedDim(out), reducedDim(alt))
})

set.seed(120000521)
test_that("fastMNN works with within-object batches and subsetting", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 2), nrow=100) # Batch 3
    combined <- cbind(B1, B2, B3)
    batches <- rep(1:3, c(ncol(B1), ncol(B2), ncol(B3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    # Subset.row behaves correctly.
    i <- sample(nrow(B1), 50)
    ref <- fastMNN(combined[i,], batch=batches, d=50, BSPARAM=ExactParam())
    out.s <- fastMNN(combined, batch=batches, d=50, subset.row=i, BSPARAM=ExactParam())
    expect_identical(reducedDim(out.s), reducedDim(ref))
    expect_equal(out.s, ref)

    # Correct.all behaves correctly.
    i <- c(1:nrow(B1), 1:10)
    ref <- fastMNN(combined, batch=batches, d=50, BSPARAM=ExactParam())
    out <- fastMNN(combined[i,], batch=batches, subset.row=1:nrow(B1), d=50, BSPARAM=ExactParam(), correct.all=TRUE)
    expect_identical(nrow(out), length(i))
    expect_identical(reducedDim(ref), reducedDim(out))
    expect_equal(as.matrix(assay(ref)[1:10,]), as.matrix(assay(out)[1:10+nrow(B1),]))
})

set.seed(120000522)
test_that("fastMNN renames within-object batches correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 2), nrow=100) # Batch 3

    rownames(B1) <- rownames(B2) <- rownames(B3) <- sprintf("GENE_%i", sample(nrow(B1)))
    colnames(B1) <- sprintf("CELL_1_%i", seq_len(ncol(B1)))
    colnames(B2) <- sprintf("CELL_2_%i", seq_len(ncol(B2)))
    colnames(B3) <- sprintf("CELL_3_%i", seq_len(ncol(B3)))

    combined <- cbind(B1, B2, B3)
    batches <- rep(1:3, c(ncol(B1), ncol(B2), ncol(B3)))

    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]

    out <- fastMNN(combined, batch=batches, BSPARAM=ExactParam())
    expect_identical(rownames(out), rownames(B1))
    expect_identical(colnames(out), colnames(combined))
})

set.seed(12000053)
test_that("fastMNN works correctly with restriction", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    B3 <- matrix(rnorm(5000, 2), nrow=100) # Batch 3
    B4 <- matrix(rnorm(8000, 2), nrow=100) # Batch 4

    # Restricted results are only directly comparable if we're not doing a PCA internally.
    # So here, we'll just check that the pairs are only formed between allowable instances. 
    restricted <- list(1:80, 1:100, 1:40, 1:50)
    ref <- fastMNN(B1, B2, B3, B4, restrict=restricted)
    CHECK_PAIRINGS(ref)

    RESTRICTED_CHECK <- function(mnn.out, restricted) {
        origin <- mnn.out$batch
        merge.info <- metadata(mnn.out)$merge.info
        for (counter in seq_len(nrow(merge.info))) {
            left <- merge.info$left[[counter]]
            right <- merge.info$right[[counter]]
            p <- merge.info$pairs[[counter]]

            allowed.left <- unlist(lapply(left, function(x) which(x==origin)[restricted[[x]]]))
            expect_true(length(p$left) > 0)
            expect_true(all(p$left %in% allowed.left))
            
            allowed.right <- unlist(lapply(right, function(x) which(x==origin)[restricted[[x]]]))
            expect_true(length(p$right) > 0)
            expect_true(all(p$right %in% allowed.right))
        }
    }

    RESTRICTED_CHECK(ref, restricted)

    # Checking that restriction respects alternative orderings.
    out2 <- fastMNN(B1, B2, B3, B4, merge.order=4:1, restrict=restricted)
    RESTRICTED_CHECK(out2, restricted)

    out2 <- fastMNN(B1, B2, B3, B4, merge.order=list(list(4,1), list(2,3)), restrict=restricted)
    RESTRICTED_CHECK(out2, restricted)

    # Also behaves with a single batch.
    DY <- cbind(B1, B2, B3, B4)
    batch <- rep(1:4, c(ncol(B1), ncol(B2), ncol(B3), ncol(B4)))
    shuffle <- sample(length(batch))

    accumulated <- 0L
    single.res <- restricted
    for (i in seq_along(single.res)) {
        single.res[[i]] <- single.res[[i]] + accumulated
        accumulated <- accumulated + sum(batch==i)
    }

    out2 <- fastMNN(DY[,shuffle], batch=batch[shuffle], 
        restrict=list(shuffle %in% unlist(single.res)))

    CHECK_PAIRINGS(out2)
    expect_equal(out2$corrected, ref$corrected[shuffle,])
    expect_identical(out2$batch, as.character(ref$batch)[shuffle])
})

set.seed(12000054)
test_that("fastMNN works correctly with weighting of PCs", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    pcs <- multiBatchPCA(B1, B2, d=10, weights=c(5, 1), BSPARAM=ExactParam())
    out.pre <- reducedMNN(pcs[[1]], pcs[[2]])
    out.norm <- fastMNN(B1, B2, d=10, weights=c(5, 1), cos.norm=FALSE, BSPARAM=ExactParam())

    expect_equal(metadata(pcs)$rotation, rowData(out.norm)$rotation)
    expect_equal(out.pre$corrected, reducedDim(out.norm))
    expect_equal(out.pre$batch, out.norm$batch)

    # Also trying with a single batch.
    DY <- cbind(B1, B2)
    batch <- rep(LETTERS[1:2], c(ncol(B1), ncol(B2)))
    pcs <- multiBatchPCA(DY, batch=batch, d=10, weights=c(A=5, B=1), BSPARAM=ExactParam())
    out.pre <- reducedMNN(A=pcs[[1]], B=pcs[[2]])
    out.norm <- fastMNN(DY, batch=batch, d=10, weights=c(A=5, B=1), cos.norm=FALSE, BSPARAM=ExactParam())

    expect_equal(metadata(pcs)$rotation, rowData(out.norm)$rotation)
    expect_equal(out.pre$corrected, reducedDim(out.norm))
    expect_equal(out.pre$batch, out.norm$batch)
})

set.seed(120000541)
test_that("fastMNN works correctly with no PCA", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    no.pc <- fastMNN(B1, B2, d=NA, cos.norm=FALSE)
    center <- (rowMeans(B1) + rowMeans(B2))/2
    ref <- reducedMNN(t(B1 - center), t(B2 - center))
    expect_identical(ref$corrected, reducedDim(no.pc))

    # Same result when we combine them.
    no.pc2 <- fastMNN(cbind(B1, B2), batch=rep(1:2, c(ncol(B1), ncol(B2))), d=NA, cos.norm=FALSE)
    expect_identical(reducedDim(no.pc), reducedDim(no.pc2))

    # Subsetting behaves as expected.
    out <- fastMNN(B1, B2, d=NA, subset.row=20:5)
    ref <- fastMNN(B1[20:5,], B2[20:5,], d=NA)
    expect_identical(out, ref)

    # Correct.all behaves as expected.
    out <- fastMNN(B1, B2, d=NA, subset.row=20:5, correct.all=TRUE)
    expect_identical(reducedDim(out), reducedDim(ref))
    expect_identical(as.matrix(assay(out)[20:5,]), as.matrix(assay(ref)))
    expect_true(all(assay(out)[-(20:5),]==0))
})

set.seed(120000542)
test_that("fastMNN reports PC statistics correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2

    out <- fastMNN(B1, B2, get.variance=TRUE)
    expect_true(!is.null(metadata(out)$pca.info$var.total))
    expect_true(!is.null(metadata(out)$pca.info$var.explained))
})

set.seed(1200006)
test_that("fastMNN fails on silly inputs", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    # Throws errors properly with no genes or no cells.
    expect_error(fastMNN(), "at least two batches")
    expect_error(fastMNN(B1), "'batch' must be specified")
    expect_error(expect_warning(fastMNN(B1[0,], B2[0,], BSPARAM=ExactParam()), "more requested"), "zero")
    expect_error(expect_warning(fastMNN(B1[,0], B2[,0], BSPARAM=ExactParam()), "more requested"), "zero")

    # Throws errors upon row checks.
    expect_error(fastMNN(B1[1:10,], B2), "number of rows is not the same")
    xB1 <- B1
    xB2 <- B2
    rownames(xB1) <- sample(nrow(B1))
    rownames(xB2) <- sample(nrow(B2))
    expect_error(fastMNN(xB1, xB2), "row names are not the same")

    # Throws errors upon restrict name checks.
    expect_error(fastMNN(B=B1, A=B2, restrict=list(A=1, B=2)), "same names")
})
