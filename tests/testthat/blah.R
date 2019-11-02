    # Check with sparse matrices.
    library(Matrix)
    sp1 <- rsparsematrix(50, 25, density=0.1)
    sp2 <- rsparsematrix(100, 25, density=0.1)
    mnn1 <- 1:10
    mnn2 <- 30:21


# Constructing a reference function.
REF <- function(data1, data2, mnn1, mnn2, s2) {
    d <- as.matrix(dist(data2))
    w <- exp(-d^2/s2)

    mnn.dens <- rowSums(w[,unique(mnn2)])
    N <- tabulate(mnn2, nbins=nrow(data2))
    kernel <- t(w/(N*mnn.dens))[,mnn2]

    kernel <- kernel/rowSums(kernel)
    vect <- data1[mnn1,] - data2[mnn2,]
    out <- kernel %*% vect
    dimnames(out) <- NULL
    return(out)
}

    xx <- batchelor:::.compute_correction_vectors(sp1, sp2, mnn1, mnn2, t(sp2), 0.5)
    ref <- REF(as.matrix(sp1), as.matrix(sp2), mnn1, mnn2, 0.5)
    testthat:::expect_equal(xx, ref)
