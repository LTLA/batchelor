# Unit tests utilities that are indirectly covered by the other tests.
# library(batchelor); library(testthat); source("test-utils.R")

set.seed(1000004)
test_that(".rename_output works correctly", {
    # For columns:
    A <- B <- matrix(runif(1000), nrow=10)
    rownames(A) <- rownames(B) <- sprintf("GENE_%i", seq_len(nrow(A)))
    colnames(B) <- seq_len(ncol(B))

    output <- matrix(0, nrow=10, ncol=2*ncol(A))
    output2 <- batchelor:::.rename_output(output, list(A, B))
    expect_equivalent(output2, output)   
    expect_identical(rownames(output2), rownames(A))
    expect_identical(colnames(output2), c(character(ncol(A)), colnames(B)))

    output3 <- batchelor:::.rename_output(output[1:5,], list(A, B), subset.row=1:5)
    expect_identical(output3, output2[1:5,])
    i <- sample(rownames(A), 5)
    output4 <- batchelor:::.rename_output(output[1:5,], list(A, B), subset.row=i)
    expect_identical(output4, output2[i,])

    # For rows:
    A <- B <- matrix(runif(1000), ncol=10)
    colnames(A) <- colnames(B) <- sprintf("PC_%i", seq_len(ncol(A)))
    rownames(A) <- seq_len(nrow(A))

    output <- matrix(0, nrow=2*nrow(A), ncol=10)
    output2 <- batchelor:::.rename_output(output, list(A, B), cells.in.columns=FALSE)
    expect_equivalent(output2, output) 
    expect_identical(rownames(output2), c(rownames(A), character(nrow(B))))
    expect_identical(colnames(output2), colnames(A))

    output3 <- batchelor:::.rename_output(output, list(A, B), cells.in.columns=FALSE, subset.row=1:5)
    expect_identical(output3, output2) # no effect!
})

set.seed(1000006)
test_that(".row_subset_to_index works correctly", {
    A <- matrix(runif(1000), ncol=1)
    rownames(A) <- paste0("GENE_", seq_along(A))

    i <- sample(nrow(A), nrow(A)/2)
    expect_identical(i, batchelor:::.row_subset_to_index(A, i))

    l <- rbinom(nrow(A), 1, 0.1)==1
    expect_identical(which(l), batchelor:::.row_subset_to_index(A, l))

    s <- sample(rownames(A), nrow(A)/2)
    expect_identical(match(s, rownames(A)), batchelor:::.row_subset_to_index(A, s))

    expect_identical(batchelor:::.row_subset_to_index(A, NULL), seq_along(A))

    # Testing silly inputs.
    expect_identical(batchelor:::.row_subset_to_index(A, integer(0)), integer(0))
    expect_identical(batchelor:::.row_subset_to_index(A, logical(0)), integer(0))
    expect_identical(batchelor:::.row_subset_to_index(A, character(0)), integer(0))
})

set.seed(10000061)
test_that(".col_subset_to_index works correctly", {
    A <- matrix(runif(1000), ncol=1000)
    colnames(A) <- paste0("CELL_", seq_along(A))

    i <- sample(ncol(A), ncol(A)/2)
    expect_identical(i, batchelor:::.col_subset_to_index(A, i))

    l <- rbinom(ncol(A), 1, 0.1)==1
    expect_identical(which(l), batchelor:::.col_subset_to_index(A, l))

    s <- sample(colnames(A), ncol(A)/2)
    expect_identical(match(s, colnames(A)), batchelor:::.col_subset_to_index(A, s))

    expect_identical(batchelor:::.col_subset_to_index(A, NULL), seq_along(A))

    # Testing silly inputs.
    expect_identical(batchelor:::.col_subset_to_index(A, integer(0)), integer(0))
    expect_identical(batchelor:::.col_subset_to_index(A, logical(0)), integer(0))
    expect_identical(batchelor:::.col_subset_to_index(A, character(0)), integer(0))
})

set.seed(1000009)
test_that("tricube calculations work correctly", {
    A <- matrix(runif(1000), ncol=20)
    self <- BiocNeighbors::findKNN(A, k=10)

    out <- batchelor:::.compute_tricube_average(A, indices=self$index, distances=self$distance)
    expect_identical(dim(out), dim(A)) 

    # Logic checks.
    out <- batchelor:::.compute_tricube_average(A, indices=self$index[,1,drop=FALSE], distances=self$distance[,1,drop=FALSE])
    expect_identical(out, A[self$index[,1],])

    uni.dist <- self$distance
    uni.dist[] <- 1
    out <- batchelor:::.compute_tricube_average(A, indices=self$index, distances=uni.dist)
    assembly <- A
    for (i in seq_len(nrow(assembly))) {
        assembly[i,] <- colMeans(A[self$index[i,],])
    }
    expect_equal(out, assembly)

    uni.index <- self$index
    uni.index[] <- seq_len(nrow(A))
    out <- batchelor:::.compute_tricube_average(A, indices=uni.index, distances=self$distance)
    expect_equal(out, A)

    # Silly input checks.
    out <- batchelor:::.compute_tricube_average(A, indices=self$index[,0], distances=self$distance[,0])
    expect_identical(dim(out), dim(A))
    expect_true(all(out==0))

    out <- batchelor:::.compute_tricube_average(A[0,], indices=self$index[0,], distances=self$distance[0,])
    expect_identical(dim(out), c(0L, ncol(A)))
})

set.seed(1000010)
test_that(".restore_original_order works correctly", {
    out <- batchelor:::.restore_original_order(c(2,1,3), c(10L, 20L, 30L))
    expect_identical(out, c(21:30, 1:20, 31:60))

    original <- list(runif(35), runif(13), runif(23), runif(2), runif(42))
    s <- c(5,3,1,4,2)
    shuffled <- original[s]
    out <- batchelor:::.restore_original_order(s, lengths(original))
    expect_equal(unlist(shuffled)[out], unlist(original))

    # Testing silly inputs.
    expect_error(out <- batchelor:::.restore_original_order(integer(0), 1), "not equal")
    out <- batchelor:::.restore_original_order(integer(0), integer(0))
    expect_identical(out, integer(0))
})

set.seed(1000011)
test_that(".reindex_pairings works correctly", {
    S <- sample(40)
    pairings <- list(
        data.frame(left=sample(10, 20, replace=TRUE), right=11:30), 
        data.frame(left=30:1, right=sample(33:40, 30, replace=TRUE))
    )

    out <- batchelor:::.reindex_pairings(pairings, S)
    expect_identical(S[out[[1]]$left], pairings[[1]]$left)
    expect_identical(S[out[[1]]$right], pairings[[1]]$right)
    expect_identical(S[out[[2]]$left], pairings[[2]]$left)
    expect_identical(S[out[[2]]$right], pairings[[2]]$right)

    # Works on empty inputs.
    empty <- lapply(pairings, "[", i=0,, drop=FALSE)
    out <- batchelor:::.reindex_pairings(empty, integer(0))
    expect_identical(empty, out)
})
