# Unit tests utilities that are indirectly covered by the other tests.
# library(batchelor); library(testthat); source("test-utils.R")

set.seed(1000001)
test_that(".check_batch_consistency works correctly", {
    # Testing row checks.
    A1 <- matrix(runif(1000), nrow=10)
    A2 <- matrix(runif(2000), nrow=10)
    expect_identical(batchelor:::.check_batch_consistency(list(A1, A2)), list(NULL, list(NULL, NULL)))
    expect_identical(batchelor:::.check_batch_consistency(list(A1, A2, cbind(A1, A2))), list(NULL, list(NULL, NULL, NULL)))

    rownames(A1) <- rownames(A2) <- sample(nrow(A1))
    expect_identical(batchelor:::.check_batch_consistency(list(A1, A2)), list(rownames(A1), list(NULL, NULL)))
    
    rownames(A2) <- NULL
    expect_error(batchelor:::.check_batch_consistency(list(A1, A2)), "row names are not the same")
    expect_identical(batchelor:::.check_batch_consistency(list(A1, A2), ignore.null=TRUE), list(rownames(A1), list(NULL, NULL)))

    rownames(A2) <- sample(nrow(A1))
    expect_error(batchelor:::.check_batch_consistency(list(A1, A2)), "row names are not the same")

    # Testing column checks.
    B1 <- matrix(runif(1000), ncol=10)
    B2 <- matrix(runif(2000), ncol=10)
    expect_identical(batchelor:::.check_batch_consistency(list(B1, B2), byrow=FALSE), list(list(NULL, NULL), NULL))
    expect_identical(batchelor:::.check_batch_consistency(list(B1, B2, rbind(B1, B2)), byrow=FALSE), list(list(NULL, NULL, NULL), NULL))

    colnames(B1) <- colnames(B2) <- seq_len(ncol(B1))
    expect_identical(batchelor:::.check_batch_consistency(list(B1, B2), byrow=FALSE), list(list(NULL, NULL), colnames(B1)))
    
    colnames(B2) <- NULL
    expect_error(batchelor:::.check_batch_consistency(list(B1, B2), byrow=FALSE), "column names are not the same")
    expect_identical(batchelor:::.check_batch_consistency(list(B1, B2), byrow=FALSE, ignore.null=TRUE), list(list(NULL, NULL), colnames(B1)))

    colnames(B2) <- rev(seq_len(ncol(B1)))
    expect_error(batchelor:::.check_batch_consistency(list(B1, B2), byrow=FALSE), "column names are not the same")

    # Flipping around the dimension names.
    C1 <- A1
    colnames(C1) <- seq_len(ncol(C1))
    expected <- list(rownames(A1), list(character(ncol(A1)), colnames(C1)))
    expect_identical(batchelor:::.check_batch_consistency(list(A1, C1)), expected)
    names(expected[[2]]) <- c("X", "Y")
    expect_identical(batchelor:::.check_batch_consistency(list(X=A1, Y=C1)), expected)

    D1 <- B1
    rownames(D1) <- seq_len(nrow(D1))
    expected <- list(list(character(nrow(B1)), rownames(D1)), colnames(B1))
    expect_identical(batchelor:::.check_batch_consistency(list(B1, D1), byrow=FALSE), expected)
    names(expected[[1]]) <- c("X", "Y")
    expect_identical(batchelor:::.check_batch_consistency(list(X=B1, Y=D1), byrow=FALSE), expected)

    # Getting coverage of extreme cases.
    expect_identical(batchelor:::.check_batch_consistency(list(A1)), list(rownames(A1), list(NULL)))
    expect_identical(batchelor:::.check_batch_consistency(list(B1), byrow=FALSE), list(list(NULL), colnames(B1)))

    expect_identical(batchelor:::.check_batch_consistency(list()), list(NULL, list()))
    expect_identical(batchelor:::.check_batch_consistency(list(), byrow=FALSE), list(list(), NULL))

    expect_identical(batchelor:::.check_batch_consistency(list(A1[0,], A2[0,])), list(NULL, list(NULL, NULL)))
    expect_identical(batchelor:::.check_batch_consistency(list(B1[,0], B2[,0]), byrow=FALSE), list(list(NULL, NULL), NULL))
})

set.seed(1000002)
test_that(".check_spike_consistency works correctly", {
    sce1 <- SingleCellExperiment(list(logcounts=matrix(runif(5000), nrow=10)))
    sce2 <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))
    expect_error(batchelor:::.check_spike_consistency(list(sce1, sce2)), NA)
    expect_error(batchelor:::.check_spike_consistency(list(sce1, sce2, cbind(sce1, sce2))), NA)

    isSpike(sce1, "ERCC") <- isSpike(sce2, "ERCC") <- 1:5
    expect_error(batchelor:::.check_spike_consistency(list(sce1, sce2)), NA)

    isSpike(sce2, "ERCC") <- NULL
    expect_error(batchelor:::.check_spike_consistency(list(sce1, sce2)), "spike-in sets differ")

    isSpike(sce2, "SIRV") <- isSpike(sce1, "ERCC")
    expect_error(batchelor:::.check_spike_consistency(list(sce1, sce2)), "spike-in sets differ")
    
    expect_error(batchelor:::.check_spike_consistency(list(sce1)), NA)
    expect_error(batchelor:::.check_spike_consistency(list()), NA)
})

set.seed(1000003)
test_that(".check_if_SCEs works correctly", {
    sce1 <- SingleCellExperiment(list(logcounts=matrix(runif(5000), nrow=10)))
    sce2 <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))

    expect_true(batchelor:::.check_if_SCEs(list(sce1, sce2)))
    expect_false(batchelor:::.check_if_SCEs(list(assay(sce1), assay(sce2))))

    expect_error(batchelor:::.check_if_SCEs(list(assay(sce1), sce2)), "cannot mix")
    expect_error(batchelor:::.check_if_SCEs(list(sce1, assay(sce2))), "cannot mix")

    expect_true(batchelor:::.check_if_SCEs(list(sce1)))
    expect_false(batchelor:::.check_if_SCEs(list(assay(sce1))))
    expect_false(batchelor:::.check_if_SCEs(list()))
})

set.seed(1000005)
test_that(".create_batch_names works correctly", {
    out <- batchelor:::.create_batch_names(LETTERS[1:3], c(10, 20, 30))
    expect_identical(out$labels, LETTERS[1:3])
    expect_identical(out$ids, rep(out$labels, c(10, 20, 30)))
    
    out <- batchelor:::.create_batch_names(NULL, c(10, 20, 30))
    expect_identical(out$ids, rep(out$labels, c(10, 20, 30)))
})

set.seed(10000051)
test_that(".check_restrictions works correctly", {
    expect_error(batchelor:::.check_restrictions(list(1, 2), NULL), NA)
    expect_error(batchelor:::.check_restrictions(list(1, 2), list()), "number of batches")
    expect_error(batchelor:::.check_restrictions(list(1, 2), list(3, 4)), NA)
    expect_error(batchelor:::.check_restrictions(list(A=1, B=2), list(3, 4)), "same names")
    expect_error(batchelor:::.check_restrictions(list(A=1, B=2), list(B=3, A=4)), "same names")
    expect_error(batchelor:::.check_restrictions(list(A=1, B=2), list(A=3, B=4)), NA)
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

set.seed(1000007)
test_that(".spike_subset works correctly", {
    sce <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))

    # Behaves without any spike-ins.
    expect_identical(batchelor:::.spike_subset(sce, TRUE), NULL)
    expect_identical(batchelor:::.spike_subset(sce, FALSE), NULL)

    # Behaves once you throw spike-ins in.
    isSpike(sce, "ERCC") <- 1:5
    expect_identical(batchelor:::.spike_subset(sce, TRUE), NULL)
    expect_identical(batchelor:::.spike_subset(sce, FALSE), !isSpike(sce))

    # Testing silly inputs.
    expect_identical(batchelor:::.spike_subset(sce[0,], TRUE), NULL)
    expect_identical(batchelor:::.spike_subset(sce[0,], FALSE), NULL)
})

set.seed(1000008)
test_that(".SCE_subset_genes works correctly", {
    sce <- SingleCellExperiment(list(logcounts=matrix(runif(2000), nrow=10)))

    # Behaves without any spike-ins.
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce, TRUE), NULL)
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce, FALSE), NULL)

    i <- 1:5
    expect_identical(batchelor:::.SCE_subset_genes(i, sce, TRUE), i)
    expect_identical(batchelor:::.SCE_subset_genes(i, sce, FALSE), i)

    l <- seq_len(nrow(sce)) %% 2 == 0
    expect_identical(batchelor:::.SCE_subset_genes(l, sce, TRUE), which(l))
    expect_identical(batchelor:::.SCE_subset_genes(l, sce, FALSE), which(l))

    # Behaves once you throw spike-ins in.
    isSpike(sce, "ERCC") <- 1:5
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce, TRUE), NULL)
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce, FALSE), !isSpike(sce))

    i <- 3:10
    expect_identical(batchelor:::.SCE_subset_genes(i, sce, TRUE), i)
    expect_identical(batchelor:::.SCE_subset_genes(i, sce, FALSE), 6:10)

    l <- seq_len(nrow(sce)) %% 2 == 0
    expect_identical(batchelor:::.SCE_subset_genes(l, sce, TRUE), which(l))
    expect_identical(batchelor:::.SCE_subset_genes(l, sce, FALSE), setdiff(which(l), 1:5))

    # Testing silly inputs.
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce[0,], TRUE), NULL)
    expect_identical(batchelor:::.SCE_subset_genes(NULL, sce[0,], FALSE), NULL)

    expect_identical(batchelor:::.SCE_subset_genes(integer(0), sce[0,], TRUE), integer(0))
    expect_identical(batchelor:::.SCE_subset_genes(integer(0), sce[0,], FALSE), integer(0))
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
    pairings <- list(data.frame(first=sample(10, 20, replace=TRUE), second=11:30), data.frame(first=30:1, second=sample(33:40, 30, replace=TRUE)))

    out <- batchelor:::.reindex_pairings(pairings, S)
    expect_identical(S[out[[1]]$first], pairings[[1]]$first)
    expect_identical(S[out[[1]]$second], pairings[[1]]$second)
    expect_identical(S[out[[2]]$first], pairings[[2]]$first)
    expect_identical(S[out[[2]]$second], pairings[[2]]$second)

    # Works on empty inputs.
    empty <- lapply(pairings, "[", i=0,, drop=FALSE)
    out <- batchelor:::.reindex_pairings(empty, integer(0))
    expect_identical(empty, out)
})
