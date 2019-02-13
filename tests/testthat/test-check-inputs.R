# Unit tests the input checking utilities.
# library(batchelor); library(testthat); source("test-check-inputs.R")

set.seed(1000001)
test_that("checkBatchConsistency works correctly", {
    # Testing row checks.
    A1 <- matrix(runif(1000), nrow=10)
    A2 <- matrix(runif(2000), nrow=10)
    expect_identical(checkBatchConsistency(list(A1, A2)), list(NULL, list(NULL, NULL)))
    expect_identical(checkBatchConsistency(list(A1, A2, cbind(A1, A2))), list(NULL, list(NULL, NULL, NULL)))

    rownames(A1) <- rownames(A2) <- sample(nrow(A1))
    expect_identical(checkBatchConsistency(list(A1, A2)), list(rownames(A1), list(NULL, NULL)))
    
    rownames(A2) <- NULL
    expect_error(checkBatchConsistency(list(A1, A2)), "row names are not the same")
    expect_identical(checkBatchConsistency(list(A1, A2), ignore.null=TRUE), list(rownames(A1), list(NULL, NULL)))

    rownames(A2) <- sample(nrow(A1))
    expect_error(checkBatchConsistency(list(A1, A2)), "row names are not the same")

    # Testing column checks.
    B1 <- matrix(runif(1000), ncol=10)
    B2 <- matrix(runif(2000), ncol=10)
    expect_identical(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), list(list(NULL, NULL), NULL))
    expect_identical(checkBatchConsistency(list(B1, B2, rbind(B1, B2)), cells.in.columns=FALSE), list(list(NULL, NULL, NULL), NULL))

    colnames(B1) <- colnames(B2) <- seq_len(ncol(B1))
    expect_identical(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), list(list(NULL, NULL), colnames(B1)))
    
    colnames(B2) <- NULL
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), "column names are not the same")
    expect_identical(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE, ignore.null=TRUE), list(list(NULL, NULL), colnames(B1)))

    colnames(B2) <- rev(seq_len(ncol(B1)))
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), "column names are not the same")
})

set.seed(10000011)
test_that("checkBatchConsistency fills in NULL character names", {
    # For columns:
    A <- B <- matrix(runif(1000), nrow=10)
    rownames(A) <- rownames(B) <- sprintf("GENE_%i", seq_len(nrow(A)))
    colnames(B) <- seq_len(ncol(A))

    expected <- list(rownames(A), list(character(ncol(A)), colnames(B)))
    expect_identical(checkBatchConsistency(list(A, B)), expected)
    names(expected[[2]]) <- c("X", "Y")
    expect_identical(checkBatchConsistency(list(X=A, Y=B)), expected)

    # For rows:
    C <- D <- matrix(runif(2000), nrow=10)
    colnames(C) <- colnames(D) <- sprintf("CELL_%i", seq_len(ncol(C)))
    rownames(D) <- seq_len(nrow(D))

    expected <- list(list(character(nrow(C)), rownames(D)), colnames(C))
    expect_identical(checkBatchConsistency(list(C, D), cells.in.columns=FALSE), expected)
    names(expected[[1]]) <- c("X", "Y")
    expect_identical(checkBatchConsistency(list(X=C, Y=D), cells.in.columns=FALSE), expected)
})

set.seed(10000012)
test_that("checkBatchConsistency fills in character names", {
    A <- matrix(runif(1000), nrow=10)
    B <- matrix(runif(2000), nrow=10)

    # Getting coverage of extreme cases.
    expect_identical(checkBatchConsistency(list(A)), list(NULL, list(NULL)))
    expect_identical(checkBatchConsistency(list(B), cells.in.columns=FALSE), list(list(NULL), NULL))

    expect_identical(checkBatchConsistency(list()), list(NULL, list()))
    expect_identical(checkBatchConsistency(list(), cells.in.columns=FALSE), list(list(), NULL))

    expect_identical(checkBatchConsistency(list(A[0,], A[0,])), list(NULL, list(NULL, NULL)))
    expect_identical(checkBatchConsistency(list(B[,0], B[,0]), cells.in.columns=FALSE), list(list(NULL, NULL), NULL))
})

set.seed(1000002)
test_that("checkSpikeConsistency works correctly", {
    sce1 <- SingleCellExperiment(list(logcounts=matrix(runif(5000), nrow=10)))
    sce2 <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))
    expect_error(checkSpikeConsistency(list(sce1, sce2)), NA)
    expect_error(checkSpikeConsistency(list(sce1, sce2, cbind(sce1, sce2))), NA)

    isSpike(sce1, "ERCC") <- isSpike(sce2, "ERCC") <- 1:5
    expect_error(checkSpikeConsistency(list(sce1, sce2)), NA)

    isSpike(sce2, "ERCC") <- NULL
    expect_error(checkSpikeConsistency(list(sce1, sce2)), "spike-in sets differ")

    isSpike(sce2, "SIRV") <- isSpike(sce1, "ERCC")
    expect_error(checkSpikeConsistency(list(sce1, sce2)), "spike-in sets differ")
    
    expect_error(checkSpikeConsistency(list(sce1)), NA)
    expect_error(checkSpikeConsistency(list()), NA)
})

set.seed(1000003)
test_that("checkIfSCE works correctly", {
    sce1 <- SingleCellExperiment(list(logcounts=matrix(runif(5000), nrow=10)))
    sce2 <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))

    expect_true(checkIfSCE(list(sce1, sce2)))
    expect_false(checkIfSCE(list(assay(sce1), assay(sce2))))

    expect_error(checkIfSCE(list(assay(sce1), sce2)), "cannot mix")
    expect_error(checkIfSCE(list(sce1, assay(sce2))), "cannot mix")

    expect_true(checkIfSCE(list(sce1)))
    expect_false(checkIfSCE(list(assay(sce1))))
    expect_false(checkIfSCE(list()))
})

set.seed(10000051)
test_that("checkRestrictions works correctly", {
    # Testing row checks.
    A1 <- matrix(runif(1000), nrow=10)
    A2 <- matrix(runif(2000), nrow=10)

    expect_identical(checkRestrictions(list(A1, A2), NULL), NULL)
    keep1 <- 1:10
    keep2 <- rbinom(ncol(A2), 1, 0.5)==1
    expect_identical(checkRestrictions(list(A1, A2), list(keep1, keep2)), list(keep1, which(keep2)))
    expect_identical(checkRestrictions(list(A1, A2), list(NULL, keep2)), list(NULL, which(keep2)))

    colnames(A2) <- sprintf("CELL_%i", seq_len(ncol(A2)))
    keep3 <- sample(colnames(A2))
    expect_identical(checkRestrictions(list(A1, A2), list(NULL, keep3)), list(NULL, match(keep3, colnames(A2))))

    # Testing column checks.
    B1 <- matrix(runif(1000), ncol=10)
    B2 <- matrix(runif(2000), ncol=10)

    expect_identical(checkRestrictions(list(B1, B2), NULL, cells.in.columns=FALSE), NULL)
    keep1 <- 1:10
    keep2 <- rbinom(nrow(B2), 1, 0.5)==1
    expect_identical(checkRestrictions(list(B1, B2), cells.in.columns=FALSE, list(keep1, keep2)), list(keep1, which(keep2)))
    expect_identical(checkRestrictions(list(B1, B2), cells.in.columns=FALSE, list(NULL, keep2)), list(NULL, which(keep2)))

    rownames(B2) <- sprintf("CELL_%i", seq_len(nrow(B2)))
    keep3 <- sample(rownames(B2), 10)
    expect_identical(checkRestrictions(list(B1, B2), cells.in.columns=FALSE, list(NULL, keep3)), list(NULL, match(keep3, rownames(B2))))

    # Throwing upon errors.
    expect_error(checkRestrictions(list(A1, A2), list()), "number of batches")
    expect_error(checkRestrictions(list(A=A1, B=A2), list(3, 4)), "same names")
    expect_error(checkRestrictions(list(A=A1, B=A2), list(B=3, A=4)), "same names")
    expect_error(checkRestrictions(list(A1, A2), list(3, integer(0))), "no cells")

    expect_identical(checkRestrictions(list(A=A1, B=A2), list(A=3, B=4)), list(A=3L, B=4L)) # handles names.
})
