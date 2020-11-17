# Tests the noCorrect() function.
# library(testthat); library(batchelor); source('test-no-correct.R')

test_that("noCorrect works as expected", {
    test1 <- matrix(rnorm(1000), nrow=100)
    test2 <- matrix(rnorm(2000), nrow=100)

    out <- noCorrect(test1, test2)
    expect_equivalent(assay(out), cbind(test1, test2))
    expect_equivalent(out$batch, rep(1:2, c(ncol(test1), ncol(test2))))

    out <- noCorrect(X=test1, Y=test2)
    expect_equivalent(assay(out), cbind(test1, test2))
    expect_equivalent(out$batch, rep(c("X", "Y"), c(ncol(test1), ncol(test2))))

    # Handles three batches.
    test3 <- matrix(rnorm(3000), nrow=100)
    out <- noCorrect(X=test1, Y=test2, Z=test3)
    expect_equivalent(assay(out), cbind(test1, test2, test3))
    expect_equivalent(out$batch, rep(c("X", "Y", "Z"), c(ncol(test1), ncol(test2), ncol(test3))))

    # Handles row names.
    rownames(test1) <- rownames(test2) <- rownames(test3) <- sprintf("GENE_%i", seq_len(nrow(test1)))
    out <- noCorrect(test1, test2, test3)
    expect_equivalent(rownames(out), rownames(test1))
})

test_that("noCorrect responds to subsetting", {
    test1 <- matrix(rnorm(1000), nrow=100)
    test2 <- matrix(rnorm(2000), nrow=100)

    subset.row <- 1:10
    out <- noCorrect(X=test1, Y=test2, subset.row=subset.row)
    ref <- noCorrect(X=test1[subset.row,], Y=test2[subset.row,])
    expect_equal(out, ref)

    # Handles subsetted row names.
    rownames(test1) <- rownames(test2) <- sprintf("GENE_%i", seq_len(nrow(test1)))
    out <- noCorrect(X=test1, Y=test2, subset.row=subset.row)
    ref <- noCorrect(X=test1[subset.row,], Y=test2[subset.row,])
    expect_equal(out, ref)

    # Handles correct.all=TRUE.
    out <- noCorrect(X=test1, Y=test2, subset.row=subset.row, correct.all=TRUE)
    ref <- noCorrect(X=test1, Y=test2)
    expect_equal(out, ref)
})

test_that("noCorrect handles a single object", {
    test <- matrix(rnorm(10000), nrow=100)
    batch <- sample(LETTERS[1:5], ncol(test), replace=TRUE)

    out <- noCorrect(test, batch=batch)
    expect_equivalent(assay(out), test)
    expect_identical(batch, out$batch)

    # Interacts with subsetting.
    subset.row <- sample(nrow(test), 20)
    out <- noCorrect(test, batch=batch, subset.row=subset.row)
    expect_equivalent(assay(out), test[subset.row,])
    expect_identical(batch, out$batch)
})

test_that("noCorrect handles SingleCellExperiment inputs", {
    test1 <- SingleCellExperiment(list(logcounts=matrix(rnorm(1000), nrow=100)))
    test2 <- SingleCellExperiment(list(logcounts=matrix(rnorm(2000), nrow=100)))

    out <- noCorrect(test1, test2)
    ref <- noCorrect(assay(test1), assay(test2))
    expect_equal(out, ref)

    assayNames(test1) <- assayNames(test2) <- "whee"
    out2 <- noCorrect(test1, test2, assay.type="whee")
    expect_equal(out2, ref)
})

test_that("noCorrect respects names properly", {
    test <- matrix(rnorm(10000), nrow=100)
    colnames(test) <- seq_len(nrow(test))
    batch <- sample(LETTERS[1:5], ncol(test), replace=TRUE)
    names(batch) <- paste0("whee", seq_len(nrow(test)))

    out <- noCorrect(test, batch=batch)
    expect_identical(colnames(out), colnames(test))
    expect_identical(out$batch , batch)
})
