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
   
        nn21 <- BiocNeighbors::queryKNN(d2, query=d1, k=k2)
        nn12 <- BiocNeighbors::queryKNN(d1, query=d2, k=k1)

        all.pairs.1 <- all.pairs.2 <- vector("list", n1)

        for (i in seq_len(n1)) {
            neighbors <- nn21$index[i,]
            converse <- nn12$index[neighbors,,drop=FALSE]
            mutual <- rowSums(converse==i) > 0L
            all.pairs.1[[i]] <- rep(i, sum(mutual))
            all.pairs.2[[i]] <- neighbors[mutual]
        }

        list(first=unlist(all.pairs.1), second=unlist(all.pairs.2))
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
    comparator(REF(A, B, 1, 1), findMutualNN(A, B, 1, 1))

    A <- matrix(rpois(25000, lambda=20), ncol=100)
    B <- matrix(rpois(15000, lambda=50), ncol=100)
    comparator(REF(A, B, 10, 10), findMutualNN(A, B, 10, 10))
    comparator(REF(A, B, 5, 20), findMutualNN(A, B, 5, 20))
    comparator(REF(A, B, 20, 5), findMutualNN(A, B, 20, 5))
    comparator(REF(A, B, 1, 1), findMutualNN(A, B, 1, 1))
})
