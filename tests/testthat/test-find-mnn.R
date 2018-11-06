# Checks the mnn discovery in findMutualNN.
# require(batchelor); require(testthat); source("test-find-mnn.R")

set.seed(10000)
test_that("Mutual NN detection is correct", {
    # Reference NNs.
    library(Matrix)
    REF <- function(d1, d2, k1, k2) {
        n1 <- nrow(d1)
        n2 <- nrow(d2)
        n.total <- n1 + n2
   
        W21 <- BiocNeighbors::queryKNN(d2, query=d1, k=k2)
        W12 <- BiocNeighbors::queryKNN(d1, query=d2, k=k1)
        W <- sparseMatrix(i=c(rep(seq_len(n1), k2), rep(n1 + seq_len(n2), k1)),
                          j=c(n1 + W21$index, W12$index),
                          x=rep(1, n1*k2 + n2*k1), dims=c(n.total, n.total))

        W <- W * t(W) # elementwise multiplication to keep mutual nns only
        A <- which(W>0, arr.ind=TRUE) # row/col indices of mutual NNs

        A1 <- A[,1]
        A1 <- A1[A1 <= n1]
        A2 <- A[,2] - n1
        A2 <- A2[A2 > 0]
        return(list(first=A1, second=A2))
    }
    
    # Check that the values are identical.
    comparator <- function(x, y) {
        ox <- order(x$first, x$second)
        oy <- order(y$first, y$second)
        expect_identical(x$first[ox], y$first[oy])
        expect_identical(x$second[ox], y$second[oy]) 
    }
    
    # Compare to actual run.
    A <- matrix(rnorm(10000), ncol=50)
    B <- matrix(rnorm(20000), ncol=50)
    comparator(REF(A, B, 10, 10), findMutualNN(A, B, 10, 10))
    comparator(REF(A, B, 5, 20), findMutualNN(A, B, 5, 20))
    comparator(REF(A, B, 20, 5), findMutualNN(A, B, 20, 5))

    A <- matrix(rpois(25000, lambda=20), ncol=100)
    B <- matrix(rpois(15000, lambda=50), ncol=100)
    comparator(REF(A, B, 10, 10), findMutualNN(A, B, 10, 10))
    comparator(REF(A, B, 5, 20), findMutualNN(A, B, 5, 20))
    comparator(REF(A, B, 20, 5), findMutualNN(A, B, 20, 5))
})
