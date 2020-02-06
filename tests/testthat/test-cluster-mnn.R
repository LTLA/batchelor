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
    ref <- fastMNN(assay(norm1), assay(norm2), cos.norm=FALSE, k=1, 
        BSPARAM=BiocSingular::ExactParam())

    expect_identical(metadata(output)$merge.info, metadata(ref)$merge.info)
    expect_identical(metadata(output)$cluster$cluster, c(colnames(norm1), colnames(norm2)))
    expect_identical(metadata(output)$cluster$batch, rep(1:2, c(ncol(norm1), ncol(norm2))))
})

test_that("clusterMNN's full-rank PCA preserves relative distances", {
    stuff1 <- matrix(rnorm(1000), nrow=20)
    stuff2 <- matrix(rnorm(500), nrow=20)
    stuff3 <- matrix(rnorm(2000), nrow=20)
    
    pcas <- batchelor:::.full_rank_pca(list(stuff1, stuff2, stuff3), correct.all=FALSE, subset.row=NULL,
        BSPARAM=BiocSingular::ExactParam(), BPPARAM=BiocParallel::SerialParam())
    out <- dist(do.call(rbind, as.list(pcas)))

    ref <- dist(t(cbind(stuff1, stuff2, stuff3)))
    expect_equal(as.matrix(out), as.matrix(ref))
})

test_that("clusterMNN performs the careful gaussian smoothing correctly", {
    pcs <- matrix(rnorm(1000), ncol=20)
    centers <- matrix(rnorm(200), ncol=20)
    delta <- matrix(rnorm(200), ncol=20) - centers

    output <- batchelor:::.smooth_gaussian_from_centroids(pcs, 
        centers=centers, sigma=0.5, delta=delta)

    # Computing it the naive way:
    distances <- matrix(0, nrow(pcs), nrow(centers))
    for (i in seq_len(ncol(distances))) {
        distances[,i] <- colSums((t(pcs) - centers[i,])^2)
    }

    w <- exp(-distances/0.5^2)
    w <- w/rowSums(w)
    ref <- pcs + w %*% delta

    expect_equal(ref, output)
})

set.seed(10000002)
test_that("clusterMNN's propagation of the correction vectors works correctly", {
    cluster <- rep(1:3, each=100)
    y <- matrix(cluster, ncol=300, nrow=50, byrow=TRUE)

    # This set of tests relies on the sigma being small enough that 
    # each cell is really only affected by its closest cluster.
    y1 <- y + runif(length(y), -0.01, 0.01)
    y2 <- y + runif(length(y), -0.01, 0.01)

    # Introduce an orthogonal batch effect to avoid edge case failure from fastMNN.
    y1 <- rbind(y1, 0)
    y2 <- rbind(y2, 1000)

    # Turn off cosine normalization to avoid weird effects from the 1000 above.
    out <- clusterMNN(y1, y2, cos.norm=FALSE, clusters=list(cluster, cluster))

    # Check that the cluster means are the same between batches after
    # correction, indicating the propagation of the per-center batch correction
    # to to per-cell corrected coordinates was successful.
    rd <- reducedDim(out)
    for (i in unique(cluster)) {
        LEFT <- rd[out$batch==1 & out$cluster==i,]
        RIGHT <- rd[out$batch==2 & out$cluster==i,]
        expect_equal(colMeans(LEFT), colMeans(RIGHT))
    }
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
