# This tests the functionality around applyMultiSCE.
# library(testthat); library(batchelor); source("test-apply-multi.R")

# Setting up some objects with alternative Experiments.
d1 <- matrix(rnbinom(5000, mu=10, size=1), ncol=100)
sce1 <- SingleCellExperiment(list(counts=d1))
sizeFactors(sce1) <- runif(ncol(d1))

d1 <- matrix(rnbinom(5000, mu=10, size=1), ncol=100)
altExp(sce1, "Spike") <- SingleCellExperiment(list(counts=d1))
d1 <- matrix(rnbinom(5000, mu=10, size=1), ncol=100)
altExp(sce1, "Protein") <- SingleCellExperiment(list(counts=d1))

d2 <- matrix(rnbinom(2000, mu=50, size=1), ncol=40)
sce2 <- SingleCellExperiment(list(counts=d2))
sizeFactors(sce2) <- runif(ncol(d2))

d2 <- matrix(rnbinom(2000, mu=50, size=1), ncol=40)
altExp(sce2, "Spike") <- SingleCellExperiment(list(counts=d2))
d2 <- matrix(rnbinom(2000, mu=50, size=1), ncol=40)
altExp(sce2, "Protein") <- SingleCellExperiment(list(counts=d2))

TESTFUN_ALPHA <- function(..., n=6) {
    lapply(list(...), FUN=head, n=n)
}

test_that("applyMultiSCE works as expected for multiple objects", {
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA)
    ref <- TESTFUN_ALPHA(sce1, sce2)
    expect_identical(
        lapply(normed, removeAltExps),
        lapply(ref, removeAltExps)
    )
    expect_identical(
        lapply(normed, altExp),
        TESTFUN_ALPHA(altExp(sce1), altExp(sce2))
    )
    expect_identical(
        lapply(normed, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2))
    )

    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=c("Spike", "Protein"))
    expect_identical(normed, normed2)

    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=1:2)
    expect_identical(normed, normed2)

    # Responds to various settings to enable/disable .
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=character(0))
    ref <- TESTFUN_ALPHA(sce1, sce2)
    expect_identical(normed, ref)

    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=2)
    expect_identical(
        lapply(normed, removeAltExps),
        lapply(ref, removeAltExps)
    )
    expect_identical(
        lapply(normed, altExp),
        list(altExp(sce1), altExp(sce2))
    )
    expect_identical(
        lapply(normed, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2))
    )

    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH="Protein")
    expect_identical(normed, normed2)

    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, MAIN.ARGS=NULL)
    expect_identical(vapply(normed, nrow, 0L), integer(2))
    expect_identical(
        lapply(normed, altExp),
        TESTFUN_ALPHA(altExp(sce1), altExp(sce2))
    )
    expect_identical(
        lapply(normed, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2))
    )
})

test_that("applyMultiSCE handles arguments correctly", {
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, n=10)
    ref <- TESTFUN_ALPHA(sce1, sce2, n=10)
    expect_identical(
        lapply(normed, removeAltExps),
        lapply(ref, removeAltExps)
    )
    expect_identical(
        lapply(normed, altExp),
        TESTFUN_ALPHA(altExp(sce1), altExp(sce2), n=10)
    )
    expect_identical(
        lapply(normed, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2), n=10)
    )

    # Works the same in the common arguments. 
    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, COMMON.ARGS=list(n=10))
    expect_identical(normed, normed2)

    # Main only.
    ref <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA)
    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, MAIN.ARGS=list(n=10))
    expect_identical(
        lapply(normed, removeAltExps),
        lapply(normed2, removeAltExps)
    )
    expect_identical(lapply(ref, altExps), lapply(normed2, altExps))

    ref <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, n=10)
    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, n=10, MAIN.ARGS=list(n=5)) # overriding.
    expect_identical(
        lapply(TESTFUN_ALPHA(sce1, sce2, n=5), removeAltExps),
        lapply(normed2, removeAltExps),
    )
    expect_identical(lapply(ref, altExps), lapply(normed2, altExps))

    # One altExp only with overriding behavior.
    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, n=5, ALT.ARGS=list(Spike=list(n=10)))
    ref <- TESTFUN_ALPHA(sce1, sce2, n=5)
    expect_identical(
        lapply(normed2, removeAltExps),
        lapply(ref, removeAltExps)
    )
    expect_identical(
        lapply(normed2, altExp),
        TESTFUN_ALPHA(altExp(sce1), altExp(sce2), n=10)
    )
    expect_identical(
        lapply(normed2, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2), n=5)
    )

    # Multiple altExps.
    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, ALT.ARGS=list(Spike=list(n=10), Protein=list(n=2)))
    ref <- TESTFUN_ALPHA(sce1, sce2)
    expect_identical(
        lapply(normed2, removeAltExps),
        lapply(ref, removeAltExps)
    )
    expect_identical(
        lapply(normed2, altExp),
        TESTFUN_ALPHA(altExp(sce1), altExp(sce2), n=10)
    )
    expect_identical(
        lapply(normed2, altExp, e=2),
        TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2), n=2)
    )
})

test_that("applyMultiSCE respects a lack of simplification", {
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, SIMPLIFY=FALSE)
    expect_identical(names(normed), c("", altExpNames(sce1)))
    expect_identical(normed[[1]], TESTFUN_ALPHA(sce1, sce2))
    expect_identical(normed[["Spike"]], TESTFUN_ALPHA(altExp(sce1), altExp(sce2)))
    expect_identical(normed[["Protein"]], TESTFUN_ALPHA(altExp(sce1, 2), altExp(sce2, 2)))

    # Automatically declines to simplify.
    anames <- applyMultiSCE(sce1, sce2, FUN=function(...) lapply(list(...), assayNames))
    cnames <- list("counts", "counts")
    expect_identical(anames, list(cnames, Spike=cnames, Protein=cnames))

    # Emits the right warnings.
    set.seed(100)
    expect_warning(sub <- applyMultiSCE(sce1, sce2, FUN=function(...) { lapply(list(...), function(y) { y[,seq_len(sample(ncol(y), 1))] }) }),
        "could not simplify")
    expect_identical(length(sub), 3L)
    expect_warning(sub <- applyMultiSCE(sce1, sce2, FUN=function(...) { list(...)[seq_len(sample(0:2, 1))] }),
        "variable")
    expect_identical(length(sub), 3L)
})

TESTFUN_BETA <- function(..., n=6) {
    strip <- lapply(list(...), removeAltExps)
    out <- lapply(strip, head, n=n)
    do.call(cbind, out)
}

test_that("applyMultiSCE handles simplification to a single object", {
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_BETA, n=10)
    ref <- TESTFUN_BETA(sce1, sce2, n=10)
    expect_identical(
        removeAltExps(normed),
        removeAltExps(ref)
    )
    expect_identical(
        altExp(normed),
        TESTFUN_BETA(altExp(sce1), altExp(sce2), n=10)
    )
    expect_identical(
        altExp(normed, 2),
        TESTFUN_BETA(altExp(sce1, 2), altExp(sce2, 2), n=10)
    )
})

test_that("applyMultiSCE handles variable altExps", {
    # Only keeps the ones present in all of them.
    altExpNames(sce1)[1] <- c("WHAT")
    normed <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA)
    ref <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH="Protein")
    expect_identical(normed, ref)
    expect_identical(altExp(sce1), altExp(normed[[1]]))
    expect_identical(altExp(sce2), altExp(normed[[2]]))

    normed2 <- applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=2)
    expect_identical(normed2, ref)

    expect_error(applyMultiSCE(sce1, sce2, FUN=TESTFUN_ALPHA, WHICH=1))
})

