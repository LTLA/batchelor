# This tests the various quickCorrect functions.
# library(testthat); library(batchelor); source("test-quick-correct.R")

set.seed(100)
d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
sce1 <- SingleCellExperiment(list(counts=d1))
sizeFactors(sce1) <- runif(ncol(d1))
rownames(sce1) <- paste0("GENE", 1:500)

d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
sce2 <- SingleCellExperiment(list(counts=d2))
sizeFactors(sce2) <- runif(ncol(d2))
rownames(sce2) <- paste0("GENE", 201:700)

universe <- intersect(rownames(sce1), rownames(sce2))

set.seed(1000)
output <- quickCorrect(sce1, sce2)

test_that("quickCorrect works as expected", {
    expect_identical(rownames(output$corrected), universe)

    # Same results.
    set.seed(1000)
    normed <- multiBatchNorm(sce1[universe,], sce2[universe,])
    dec1 <- scran::modelGeneVar(normed[[1]])
    dec2 <- scran::modelGeneVar(normed[[2]])

    pre <- quickCorrect(sce1, sce2, precomputed=list(dec1, dec2))
    expect_equal(reducedDim(output$corrected), reducedDim(pre$corrected))

    expect_error(quickCorrect(sce1, sce2, precomputed=list(dec2)), "not the same")
})

test_that("quickCorrect works with a single object", {
    com <- cbind(sce1[universe,], sce2[universe,])
    b <- rep(1:2, c(ncol(sce1), ncol(sce2)))

    set.seed(1000)
    single <- quickCorrect(com, batch=b)
    expect_equal(reducedDim(output$corrected), reducedDim(single$corrected))

    # Same results.
    normed <- multiBatchNorm(com, batch=b)
    dec <- scran::modelGeneVar(normed, block=b)

    set.seed(1000)
    pre <- quickCorrect(com, batch=b, precomputed=list(dec))
    expect_equal(reducedDim(output$corrected), reducedDim(pre$corrected))

    expect_error(quickCorrect(com, batch=b, precomputed=list(dec, dec)), "not the same")
})

test_that("quickCorrect actually uses its HVGs", {
    BPPARAM <- FastMnnParam(d=10) # avoid warnings

    set.seed(0)
    pre1 <- quickCorrect(sce1, sce2, PARAM=BPPARAM, hvg.args=list(n=50))
    set.seed(0)
    pre2 <- quickCorrect(sce1, sce2, PARAM=BPPARAM, hvg.args=list(n=100))

    # Actually has an effect.
    expect_false(identical(pre1$corrected, pre2$corrected))
    expect_identical(nrow(pre1$corrected), nrow(pre2$corrected)) # negative control
q})
