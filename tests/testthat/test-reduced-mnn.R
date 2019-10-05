# This tests the reducedMNN function.
# library(testthat); library(batchelor); source("test-reduced-mnn.R")

test_that("reducedMNN works with two batches", {
    # Note that all fastMNN checks will have batches that are offset 
    # to ensure that the correction is not skipped via min.batch.effect.
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 1), nrow=100) # Batch 2
    
    # Behaves if we only use PCs.
    pcs <- multiBatchPCA(B1, B2, d=10)
    out.pre <- reducedMNN(pcs[[1]], pcs[[2]])
    out.norm <- fastMNN(B1, B2, d=10, cos.norm=FALSE, BSPARAM=BiocSingular::ExactParam())

    expect_equal(metadata(pcs)$rotation, rowData(out.norm)$rotation)
    expect_equal(out.pre$corrected, reducedDim(out.norm))
    expect_equal(out.pre$batch, out.norm$batch)
})

test_that("reducedMNN works with names", {    
    C1 <- matrix(rnorm(1000), ncol=10)
    C2 <- matrix(rnorm(2000), ncol=10)
    out <- reducedMNN(Y=C1, X=C2)
    expect_identical(out$batch, rep(c("Y", "X"), c(nrow(C1), nrow(C2))))

    rownames(C1) <- sprintf("CELL_%i", seq_len(nrow(C1)))
    rownames(C2) <- sprintf("CELL_%i", seq_len(nrow(C2)))
    out <- reducedMNN(C1, C2)
    expect_identical(rownames(out), c(rownames(C1), rownames(C2)))
    expect_identical(rownames(out), rownames(out$corrected))

    colnames(C1) <- colnames(C2) <- sprintf("PC%i", seq_len(ncol(C1)))
    out <- reducedMNN(C1, C2)
    expect_identical(colnames(out$corrected), colnames(C1))
})

set.seed(12000052)
test_that("reducedMNN uses 'prop.k' correctly", {
    B1 <- matrix(rnorm(10000, 0), ncol=100) # Batch 1
    B2 <- matrix(rnorm(10000, 1), ncol=100) # Batch 2

    ref <- reducedMNN(B1, B2)
    out <- reducedMNN(B1, B2, k=10, prop.k=20/nrow(B1))
    expect_identical(ref, out)

    # max() kicks in.
    ref <- reducedMNN(B1, B2)
    out <- reducedMNN(B1, B2, prop.k=0)
    expect_identical(ref, out)

    # Actually causes a difference in results.
    B2a <- matrix(rnorm(20000, 1), ncol=100)
    ref <- reducedMNN(B1, B2a)
    out <- reducedMNN(B1, B2a, prop.k=20/nrow(B1))
    expect_false(identical(ref, out))
})

set.seed(120000521)
test_that("reducedMNN works with within-object batches", {
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

    out <- reducedMNN(com.pcd, batch=batches)
    ref <- reducedMNN(pcd[[1]], pcd[[2]], pcd[[3]])
    expect_equal(ref$corrected[shuffle,], out$corrected)
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))
})

set.seed(120000503)
test_that("orthogonalization is done correctly in exact toy examples", {
    core <- cbind(rep(1:10, each=10), rep(1:10, 10))
    b1 <- b2 <- core
    b1[,1] <- b1[,1] + 20
    b2[,2] <- b2[,2] + 20

    # Orthogonalized on x-axis.
    out1 <- reducedMNN(core, b1, k=1)
    expect_equal(out1$corrected[,1], rep(5.5, nrow(out1)))
    expect_equal(out1$corrected[,2], c(core[,2], b1[,2]))

    # Orthogonalized on y-axis.
    out2 <- reducedMNN(core, b1, b2, k=1)
    expect_equal(out2$corrected[,1], rep(5.5, nrow(out2)))
    expect_equal(out2$corrected[,2], rep(5.5, nrow(out2)))

    # Hierarchical orthogonalization works.
    outY <- reducedMNN(core+10, b2+10, k=1)
    expect_equal(outY$corrected[,1], c(core[,1], b2[,1])+10)
    expect_equal(outY$corrected[,2], rep(15.5, nrow(outY)))

    outZ <- reducedMNN(core, b1, core+10, b2+10, merge.order=list(list(1,2), list(3,4)), k=1)
    expect_equal(outZ$corrected[,1], rep(5.5, nrow(outZ)))
    expect_equal(outZ$corrected[,2], rep(5.5, nrow(outZ)))
})

set.seed(12000053)
test_that("reducedMNN works correctly with restriction", {
    B1 <- matrix(rnorm(10000, 0), ncol=10) # Batch 1
    B2 <- matrix(rnorm(20000, 1), ncol=10) # Batch 2
    B3 <- matrix(rnorm(5000, 2), ncol=10) # Batch 3

    # Restricted results are only directly comparable if we're not doing a PCA internally.
    ref <- reducedMNN(B1, B2, B3)

    i1 <- 100:50
    C1 <- rbind(B1, B1[i1,])
    i2 <- 1:20
    C2 <- rbind(B2, B2[i2,])
    i3 <- 50:100
    C3 <- rbind(B3, B3[i3,])

    keep1 <- seq_len(nrow(B1))
    keep2 <- seq_len(nrow(B2))
    keep3 <- seq_len(nrow(B3))
    out <- reducedMNN(C1, C2, C3, restrict=list(keep1, keep2, keep3))

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

    out2 <- reducedMNN(DY[shuffle,], batch=batch[shuffle],
        restrict=list(shuffle %in% c(keep1, keep2 + nrow(C1), keep3 + nrow(C1) + nrow(C2))))
    expect_equal(out2$corrected, out$corrected[shuffle,])
    expect_identical(out2$batch, as.character(out$batch)[shuffle])
})
