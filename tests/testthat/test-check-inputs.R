# Unit tests the input checking utilities.
# library(batchelor); library(testthat); source("test-check-inputs.R")

set.seed(1000001)
test_that("checkBatchConsistency works correctly", {
    # Testing row checks.
    A1 <- matrix(runif(1000), nrow=10)
    A2 <- matrix(runif(2000), nrow=10)
    expect_error(checkBatchConsistency(list(A1, A2)), NA)
    expect_error(checkBatchConsistency(list(A1, A2, cbind(A1, A2))), NA)

    rownames(A1) <- rownames(A2) <- sample(nrow(A1))
    expect_error(checkBatchConsistency(list(A1, A2)), NA)
    
    rownames(A2) <- NULL
    expect_error(checkBatchConsistency(list(A1, A2)), "row names are not the same")

    rownames(A2) <- sample(nrow(A1))
    expect_error(checkBatchConsistency(list(A1, A2)), "row names are not the same")

    # Testing column checks.
    B1 <- matrix(runif(1000), ncol=10)
    B2 <- matrix(runif(2000), ncol=10)
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), NA)
    expect_error(checkBatchConsistency(list(B1, B2, rbind(B1, B2)), cells.in.columns=FALSE), NA)

    colnames(B1) <- colnames(B2) <- seq_len(ncol(B1))
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), NA)
    
    colnames(B2) <- NULL
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), "column names are not the same")

    colnames(B2) <- rev(seq_len(ncol(B1)))
    expect_error(checkBatchConsistency(list(B1, B2), cells.in.columns=FALSE), "column names are not the same")
})

set.seed(10000012)
test_that("checkBatchConsistency handles corner cases", {
    A <- matrix(runif(1000), nrow=10)
    B <- matrix(runif(2000), nrow=10)

    # Getting coverage of extreme cases.
    expect_error(checkBatchConsistency(list(A)), NA)
    expect_error(checkBatchConsistency(list(B), cells.in.columns=FALSE), NA)

    expect_error(checkBatchConsistency(list()), NA)
    expect_error(checkBatchConsistency(list(), cells.in.columns=FALSE), NA)

    expect_error(checkBatchConsistency(list(A[0,], A[0,])), NA)
    expect_error(checkBatchConsistency(list(B[,0], B[,0]), cells.in.columns=FALSE), NA)
})

set.seed(1000003)
test_that("checkIfSCE works correctly", {
    sce1 <- SingleCellExperiment(list(logcounts=matrix(runif(5000), nrow=10)))
    sce2 <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))

    expect_identical(checkIfSCE(list(sce1, sce2)), !logical(2))
    expect_identical(checkIfSCE(list(assay(sce1), assay(sce2))), logical(2))

    expect_identical(checkIfSCE(list(assay(sce1), sce2)), c(FALSE, TRUE))
    expect_identical(checkIfSCE(list(sce1, assay(sce2))), c(TRUE, FALSE))

    expect_true(checkIfSCE(list(sce1)))
    expect_false(checkIfSCE(list(assay(sce1))))
    expect_identical(checkIfSCE(list()), logical(0))
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
