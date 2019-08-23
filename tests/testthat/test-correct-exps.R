# This tests that correctExperiments() does its job.
# library(testthat); library(batchelor); source("test-correct-exps.R")

library(scater)
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
    merged2 <- correctExperiments(cbind(sce1, sce2), batch=b,
        PARAM=FastMnnParam(BSPARAM=BiocSingular::ExactParam()))

    expect_identical(reducedDim(merged), reducedDim(merged2))
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
        "ignoring assays")
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
        "overlapping")
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
