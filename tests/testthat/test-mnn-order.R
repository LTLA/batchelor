# Tests the MNN re-ordering.
# library(batchelor); library(testthat); source("test-mnn-order.R")

set.seed(123001)
test_that("use of a pre-supplied store works as expected", {
    batches <- list(matrix(rnorm(10000), ncol=10), matrix(rnorm(20000), ncol=10), matrix(rnorm(5000), ncol=10))

    for (it in 1:3) {
        if (it==1L) {
            store <- batchelor:::MNN_supplied_order(batches)
            ordering <- 1:3
        } else {
            if (it==2L) {
                ordering <- c(2L, 3L, 1L)
            } else {
                ordering <- 3:1
            }
            store <- batchelor:::MNN_supplied_order(batches, ordering=ordering)
        }
  
        expect_identical(batchelor:::.get_batches(store), batches)
        expect_identical(batchelor:::.get_restrict(store), NULL)
        running.ref <- batches[[ordering[1]]]

        for (j in 2:length(ordering)) {
            k <- 20 + j

            store <- batchelor:::.advance(store, k=k)
            expect_identical(batchelor:::.get_reference(store), running.ref)
            expect_identical(batchelor:::.get_reference_restrict(store), NULL)
            expect_identical(batchelor:::.get_reference_indices(store), head(ordering, j-1))

            current <- ordering[j]
            expect_identical(batchelor:::.get_current_index(store), current)
            expect_identical(batchelor:::.get_mnn_result(store), findMutualNN(running.ref, batches[[current]], k1=k, k2=k))

            store <- batchelor:::.compile(store, batches[[current]])
            running.ref <- rbind(running.ref, batches[[current]])
            expect_identical(batchelor:::.get_reference(store), running.ref)
            expect_identical(batchelor:::.get_reference_indices(store), head(ordering, j))
            expect_identical(batchelor:::.get_current_index(store), integer(0))
            expect_identical(batchelor:::.get_mnn_result(store), list())
        }
    }
})

set.seed(123002)
test_that("use of a pre-supplied store works with restriction", {
    batches <- list(matrix(rnorm(10000), ncol=10), matrix(rnorm(20000), ncol=10), matrix(rnorm(5000), ncol=10))
    restrict <- list(NULL, sample(nrow(batches[[2]]), nrow(batches[[2]])/4), sample(nrow(batches[[3]]), nrow(batches[[3]])/2))
    processed <- mapply(FUN=batchelor:::.row_subset_to_index, x=batches, index=restrict)

    for (it in 1:3) {
        if (it==1L) {
            store <- batchelor:::MNN_supplied_order(batches, restrict)
            ordering <- 1:3
        } else {
            # test the case where the NULL is not at the start.
            if (it==2L) {
                ordering <- c(2L, 3L, 1L)
            } else {
                ordering <- 3:1
            }
            store <- batchelor:::MNN_supplied_order(batches, restrict, ordering=ordering)
        }
  
        expect_identical(batchelor:::.get_batches(store), batches)
        expect_identical(batchelor:::.get_restrict(store), restrict)
        running.ref <- batches[[ordering[1]]]
        running.restrict <- processed[[ordering[1]]]

        for (j in 2:length(ordering)) {
            k <- 20 + j

            store <- batchelor:::.advance(store, k=k)
            expect_identical(batchelor:::.get_reference(store), running.ref)
            expect_identical(batchelor:::.get_reference_restrict(store), running.restrict)
            expect_identical(batchelor:::.get_reference_indices(store), head(ordering, j-1))
            
            current <- ordering[j]
            expect_identical(batchelor:::.get_current_index(store), current)

            cur.restrict <- processed[[current]]
            subset.res <- findMutualNN(running.ref[running.restrict,,drop=FALSE], 
                batches[[current]][cur.restrict,,drop=FALSE], k1=k, k2=k)
            subset.res$first <- running.restrict[subset.res$first]
            subset.res$second <- cur.restrict[subset.res$second]
            expect_identical(batchelor:::.get_mnn_result(store), subset.res)

            store <- batchelor:::.compile(store, batches[[current]])
            expect_identical(batchelor:::.get_reference_indices(store), head(ordering, j))
            expect_identical(batchelor:::.get_current_index(store), integer(0))
            expect_identical(batchelor:::.get_mnn_result(store), list())

            running.restrict <- c(running.restrict, nrow(running.ref) + cur.restrict) # before modifying running.ref!
            expect_identical(batchelor:::.get_reference_restrict(store), running.restrict)
            running.ref <- rbind(running.ref, batches[[current]])
            expect_identical(batchelor:::.get_reference(store), running.ref)
        }
    }
})

UNRESTRICT <- function(store) {
    store@reference.restrict <- NULL
    store@restrict <- NULL
    store
}

set.seed(123002)
test_that("comparing a pre-supplied run with and without restriction", {
    batches <- list(matrix(rnorm(10000), ncol=10), matrix(rnorm(20000), ncol=10), matrix(rnorm(5000), ncol=10))
    restrict <- lapply(batches, FUN=function(b) seq_len(nrow(b)))

    for (it in 1:2) {
        if (it==1L) {
            ref <- batchelor:::MNN_supplied_order(batches)
            store <- batchelor:::MNN_supplied_order(batches, restrict)
            ordering <- 1:3
        } else {
            ordering <- 3:1
            ref <- batchelor:::MNN_supplied_order(batches, ordering=ordering)
            store <- batchelor:::MNN_supplied_order(batches, restrict, ordering=ordering)
        }
  
        for (j in 2:length(ordering)) {
            k <- 20 + j

            store <- batchelor:::.advance(store, k=k)
            ref <- batchelor:::.advance(ref, k=k)
            expect_identical(UNRESTRICT(store), ref)

            curbatch <- batches[[ordering[j]]]
            store <- batchelor:::.compile(store, curbatch)
            ref <- batchelor:::.compile(ref, curbatch)
            expect_identical(UNRESTRICT(store), ref)
        }
    }
})
