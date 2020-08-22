# This tests that correctExperiments() does its job.
# library(testthat); library(batchelor); source("test-correct-exps.R")

library(scuttle)
sce1 <- mockSCE(ngenes=100)
sce2 <- mockSCE(ngenes=100)
sce1 <- logNormCounts(sce1)
sce2 <- logNormCounts(sce2)

test_that("correctExperiments works properly", {
    merged <- correctExperiments(sce1, sce2, PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_s4_class(assay(merged), "LowRankMatrix")
    expect_equivalent(counts(merged), cbind(counts(sce1), counts(sce2)))
    expect_equivalent(logcounts(merged), cbind(logcounts(sce1), logcounts(sce2)))

    expect_identical(length(merged$batch), ncol(merged))
    expect_identical(colData(merged)[,-1], rbind(colData(sce1), colData(sce2)))

    # Also works with a single batch... though this is less meaningful.
    b <- rep(1:2, c(ncol(sce1), ncol(sce2)))
    merged2 <- correctExperiments(cbind(sce1, sce2), batch=b, add.single=FALSE,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(reducedDim(merged), reducedDim(merged2))

    # Also works with packed batches.
    merged3 <- correctExperiments(list(sce1, sce2),
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(reducedDim(merged), reducedDim(merged3))
})

test_that("correctExperiments responds to combining arguments for assays", {
    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()),
        combine.assays="counts")

    expect_true('counts' %in% assayNames(merged))
    expect_false('logcounts' %in% assayNames(merged))

    # Smart enough to automatically restrict itself:
    assays(sce1) <- assays(sce1)["logcounts"]
    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_true('logcounts' %in% assayNames(merged))
    expect_false('counts' %in% assayNames(merged))

    # Throws a warning if the names overlap.
    assayNames(sce1)[1] <- assayNames(sce2)[1] <- "reconstructed"
    expect_warning(merged <- correctExperiments(sce1, sce2, assay.type="reconstructed",
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), 
        "ignoring 'assays'")
    expect_identical(assayNames(merged), "reconstructed")
    expect_s4_class(assay(merged), "LowRankMatrix")
})

test_that("correctExperiments responds to combining arguments for coldata", {
    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()),
        combine.coldata=c("Mutation_Status", "Treatment"))

    expect_true('Mutation_Status' %in% colnames(colData((merged))))
    expect_true('Treatment' %in% colnames(colData((merged))))
    expect_false('Cell_Cycle' %in% colnames(colData((merged))))

    # Smart enough to restrict itself.
    sce1$Treatment <- NULL

    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_true('Mutation_Status' %in% colnames(colData((merged))))
    expect_false('Treatment' %in% colnames(colData((merged))))
    expect_true('Cell_Cycle' %in% colnames(colData((merged))))

    # Throws a warning if the names overlap.
    sce1$batch <- sce2$batch <- FALSE
    expect_warning(merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), 
        "ignoring 'colData'")
    expect_type(merged$batch, "integer")
})

test_that("correctExperiments responds to including rowdata", {
    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_true(!is.null(rowData(merged)$rotation))
    rowData(merged)$rotation <- NULL
    expect_identical(rowRanges(merged), rowRanges(sce1))

    # Handles rowRanges correctly.
    dummy <- GRanges("chrA", IRanges(1:nrow(sce1), width=1))
    names(dummy) <- rownames(sce1)
    rowRanges(sce1) <- dummy

    expect_warning(merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), "ignoring non-identical 'rowRanges'")
    expect_true(!is.null(rowData(merged)$rotation))
    rowData(merged)$rotation <- NULL
    expect_identical(rowRanges(merged), rowRanges(sce2))

    rowRanges(sce2) <- rowRanges(sce1)
    merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_true(!is.null(rowData(merged)$rotation))
    rowData(merged)$rotation <- NULL
    expect_true(all(rowRanges(merged)==rowRanges(sce1)))

    # Handles rowData correctly.
    rowData(sce1)$blah <- runif(nrow(sce1))
    rowData(sce1)$stuff <- sample(LETTERS, nrow(sce1), replace=TRUE)
    rowData(sce2)$blah <- runif(nrow(sce2))

    expect_warning(merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), "ignoring non-identical 'blah'")

    expect_identical(rowRanges(merged)$stuff, rowData(sce1)$stuff)
    expect_true(is.null(rowRanges(merged)$blah))

    rowData(sce2)$blah <- rowData(sce1)$blah
    expect_warning(merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), NA)

    expect_identical(rowRanges(merged)$blah, rowData(sce1)$blah)

    # Handles conflicts with the batchCorrect output.
    rowData(sce2)$rotation <- rowData(sce1)$rotation <- 5
    expect_warning(merged <- correctExperiments(sce1, sce2, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())), 
        "ignoring 'rowData'")

    expect_true(!is.null(dim(rowData(merged)$rotation)))
    rowData(merged)$rotation <- NULL
    expect_true(all(rowRanges(merged)==rowRanges(sce1)))

    # Handles subset.row correctly.
    expect_warning(merged <- correctExperiments(sce1, sce2, subset.row=10:1,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())))
    expect_identical(rowData(merged)$stuff, rowData(sce1)$stuff[10:1])

    expect_warning(merged <- correctExperiments(sce1, sce2, subset.row=10:1, correct.all=TRUE,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())))
    expect_identical(rowData(merged)$stuff, rowData(sce1)$stuff)
})

test_that("correctExperiments works properly with add.single=TRUE", {
    b <- rep(1:2, c(ncol(sce1), ncol(sce2)))
    combined <- cbind(sce1, sce2)
    reducedDim(combined, "PCA") <- matrix(0, ncol(combined), 2)

    ref <- correctExperiments(combined, batch=b, add.single=FALSE,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    out <- correctExperiments(combined, batch=b, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_identical(reducedDim(ref), reducedDim(out))
    expect_identical(assays(ref), assays(out))
    expect_identical(colData(ref), colData(out))
    expect_identical(rowData(ref), rowData(out))

    # Points of difference from 'ref'.
    expect_identical(altExps(out), altExps(combined)) 
    expect_identical(reducedDim(combined, "PCA"), reducedDim(out, "PCA"))

    # Defend against overlaps.
    combined2 <- combined
    assayNames(combined2)[1] <- "merged"
    expect_warning(out2 <- correctExperiments(combined2, batch=b, assay.type="merged", PARAM=NoCorrectParam()), 
        "ignoring 'assays'")

    combined2 <- combined
    colnames(colData(combined2))[1] <- "batch"
    expect_warning(out2 <- correctExperiments(combined2, batch=b, PARAM=NoCorrectParam()), "ignoring 'colData'")

    combined2 <- combined
    reducedDimNames(combined2)[1] <- "corrected"
    expect_warning(out2 <- correctExperiments(combined2, batch=b, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam())),
        "ignoring 'reducedDims'")

    # Handles correct.all= correctly.
    sub <- correctExperiments(combined, batch=b, subset.row=1:10, 
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(rownames(sub), rownames(combined)[1:10])

    sub2 <- correctExperiments(combined, batch=b, subset.row=1:10, correct.all=TRUE,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(rownames(sub2), rownames(combined))
    expect_identical(reducedDim(sub), reducedDim(sub2))
})

test_that("correctExperiments respects other arguments", {
    # Respects subsetting.
    merged <- correctExperiments(sce1, sce2, subset.row=1:10,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(nrow(merged), 10L)

    # Respects correct.all=TRUE.
    merged2 <- correctExperiments(sce1, sce2, subset.row=1:10, correct.all=TRUE,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))
    expect_identical(nrow(merged2), nrow(sce1))

    expect_identical(reducedDim(merged), reducedDim(merged2))
})
