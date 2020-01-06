# Tests the clusterMNN() function.
# library(testthat); library(batchelor); source("setup.R"); source("test-cluster-mnn.R")

# Mocking up some data for multiple batches:
set.seed(10000001)
means <- matrix(rnorm(3000), ncol=3)
colnames(means) <- LETTERS[1:3]

B1 <- means[,sample(LETTERS[1:3], 500, replace=TRUE)]
B1 <- B1 + rnorm(length(B1))

B2 <- means[,sample(LETTERS[1:3], 500, replace=TRUE)]
B2 <- B2 + rnorm(length(B2)) + rnorm(nrow(B2)) # batch effect.

cluster1 <- kmeans(t(B1), centers=10)$cluster
cluster2 <- kmeans(t(B2), centers=10)$cluster

test_that("clusterMNN behaves like fastMNN on pseudo-bulk samples", {
    output <- clusterMNN(B1, B2, clusters=list(cluster1, cluster2))

    library(scater)
    norm1 <- sumCountsAcrossCells(cosineNorm(B1), cluster1, average=TRUE)
    norm2 <- sumCountsAcrossCells(cosineNorm(B2), cluster2, average=TRUE)
    ref <- fastMNN(norm1, norm2, cos.norm=FALSE, k=1, 
        BSPARAM=BiocSingular::ExactParam())

    expect_identical(metadata(output)$merge.info, metadata(ref)$merge.info)
    expect_identical(metadata(output)$cluster$cluster, c(colnames(norm1), colnames(norm2)))
    expect_identical(metadata(output)$cluster$batch, rep(1:2, c(ncol(norm1), ncol(norm2))))
})

test_that("clusterMNN's full-rank PCA preserves relative distances effect", {
    stuff1 <- matrix(rnorm(1000), nrow=20)
    stuff2 <- matrix(rnorm(500), nrow=20)
    stuff3 <- matrix(rnorm(2000), nrow=20)
    
    pcas <- batchelor:::.full_rank_pca(list(stuff1, stuff2, stuff3), correct.all=FALSE, subset.row=NULL,
        BSPARAM=BiocSingular::ExactParam(), BPPARAM=BiocParallel::SerialParam())
    out <- dist(do.call(rbind, as.list(pcas)))

    ref <- dist(t(cbind(stuff1, stuff2, stuff3)))
    expect_equal(as.matrix(out), as.matrix(ref))
})

test_that("clusterMNN can correct beyond the subset", {
    ref <- clusterMNN(B1, B2, clusters=list(cluster1, cluster2))

    expanded <- c(1:10, seq_len(nrow(B1)))
    subset <- 10+seq_len(nrow(B1))

    # Basic subsetting works correctly.
    output <- clusterMNN(B1[expanded,], B2[expanded,], clusters=list(cluster1, cluster2), subset.row=subset)
    expect_identical(output, ref)

    # Extrapolated correction works correctly.
    output <- clusterMNN(B1[expanded,], B2[expanded,], clusters=list(cluster1, cluster2), subset.row=subset, correct.all=TRUE)
    expect_identical(output[-(1:10),], ref)
    expect_equal(output[1:10,], ref[1:10,])
})

test_that("clusterMNN restriction works correctly", {
    ref <- clusterMNN(B1, B2, clusters=list(cluster1, cluster2))

    expanded1 <- c(1:10, seq_len(ncol(B1)))
    expanded2 <- c(1:10, seq_len(ncol(B2)))

    output <- clusterMNN(
        B1[,expanded1], B2[,expanded2], 
        clusters=list(cluster1[expanded1], cluster2[expanded2]),
        restrict=list(-(1:10), -(1:10))
    )

    # Same results as if those cells weren't included.
    keep <- c(10+seq_len(ncol(B1)), 10 + ncol(B1) + 10 + seq_len(ncol(B2)))
    expect_equal(ref, output[,keep])

    expect_equal(
        output[,c(1:10, 10 + ncol(B1) + 1:10)],
        output[,c(1:10 + 10, 10 + ncol(B1) + 10 + 1:10)]
    )
})

test_that("clusterMNN single-batch mode works correctly", {
    ref <- clusterMNN(A=B1, X=B2, clusters=list(cluster1, cluster2))

    single <- clusterMNN(cbind(B1, B2), 
        batch=rep(c("A", "X"), c(ncol(B1), ncol(B2))),
        clusters=list(c(cluster1, cluster2)))

    expect_identical(ref, single)

    # Also works for reordered structures.
    reorder <- sample(ncol(ref))

    single2 <- clusterMNN(cbind(B1, B2)[,reorder], 
        batch=rep(c("A", "X"), c(ncol(B1), ncol(B2)))[reorder],
        clusters=list(c(cluster1, cluster2)[reorder]))

    expect_equal(single[, reorder], single2)
})

test_that("clusterMNN single-batch works correctly with SCEs", {
    sce1 <- SingleCellExperiment(B1)
    sce2 <- SingleCellExperiment(B2)

    ref <- clusterMNN(B1, B2, clusters=list(cluster1, cluster2))
    out <- clusterMNN(sce1, sce2, clusters=list(cluster1, cluster2), assay.type=1)

    expect_identical(ref, out)
})
