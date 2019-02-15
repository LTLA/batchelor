# Tests the MNN re-ordering.
# library(batchelor); library(testthat); source("test-mnn-order.R")

GENERATOR <- function(ncells, ndim) {
    lapply(ncells, FUN=function(n) matrix(rnorm(n*ndim), ncol=ndim))
}

#####################################################

set.seed(123001)
test_that("use of a pre-supplied store works as expected", {
    batches <- GENERATOR(c(1000, 2000, 500), 10)

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
    ncells <- c(1000, 2000, 500)
    for (restrict.it in 1:4) {
        if (restrict.it==1L) {
            restrict <- lapply(ncells, FUN=seq_len)
        } else if (restrict.it==2L) {
            restrict <- lapply(ncells, FUN=function(n) NULL)
        } else if (restrict.it==3L) {
            restrict <- lapply(ncells, FUN=function(n) sample(n, n/2))
        } else {
            restrict <- lapply(ncells, FUN=function(n) sample(n, n/2))
            restrict[sample(length(restrict), 1)] <- list(NULL) # mix of integer and NULLs.
        }

        batches <- GENERATOR(ncells, 10)
        processed <- mapply(FUN=batchelor:::.row_subset_to_index, x=batches, index=restrict)

        ordering <- sample(3) # random order, for some variety.
        store <- batchelor:::MNN_supplied_order(batches, restrict, ordering=ordering)
  
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
    batches <- GENERATOR(c(1000, 2000, 500), 10)
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

#####################################################

compare_common <- function(x, y) {
    common <- intersect(slotNames(x), slotNames(y))
    for (i in common) expect_identical(slot(x, i), slot(y, i))
    TRUE
}

set.seed(123003)
test_that("use of an auto-ordered store works as expected", {
    for (it in 1:3) {
        if (it==1L) {
            ncells <- c(100, 1000, 2000, 500)
            ordering <- c(3L, 2L, 4L, 1L) # ncells in decreasing order.
        } else if (it==2L) {
            ncells <- c(100, 1000, 2000, 5000)
            ordering <- c(4L, 3L, 2L, 1L) 
        } else {
            ncells <- c(5000, 2000, 1000, 100)
            ordering <- c(2L, 1L, 3L, 4L) # known swap at start.
        }

        # This tests compares to the output of the true pre-supplied ordering.
        batches <- GENERATOR(ncells, 10)
        store <- batchelor:::MNN_auto_order(batches)
        ref <- batchelor:::MNN_supplied_order(batches, ordering=ordering)
        compare_common(store, ref)

        for (j in 2:length(batches)) {
            k <- 20 + j

            store <- batchelor:::.advance(store, k=k)
            ref <- batchelor:::.advance(ref, k=k)
            compare_common(store, ref)
            
            expect_identical(head(ordering, j-1), batchelor:::.get_reference_indices(store))
            current <- batchelor:::.get_current_index(store)
            expect_identical(ordering[j], current)

            store <- batchelor:::.compile(store, batches[[current]])
            ref <- batchelor:::.compile(ref, batches[[current]])
            compare_common(store, ref)

            expect_identical(head(ordering, j), batchelor:::.get_reference_indices(store))
        }
    }
})

set.seed(123004)
test_that("use of an auto-ordered store works with restriction", {
    ncells <- c(100, 1000, 2000, 500)
    ordering <- c(3L, 2L, 4L, 1L) # ncells in decreasing order.

    for (restrict.it in 1:4) {
        if (restrict.it==1L) {
            restrict <- lapply(ncells, FUN=seq_len)
        } else if (restrict.it==2L) {
            restrict <- lapply(ncells, FUN=function(n) NULL)
        } else if (restrict.it==3L) {
            restrict <- lapply(ncells, FUN=function(n) sample(n, n/2))
        } else {
            restrict <- lapply(ncells, FUN=function(n) sample(n, n/2))
            restrict[sample(length(restrict), 1)] <- list(NULL) # mix of integer and NULLs.
        }

        batches <- GENERATOR(ncells, 10)
        store <- batchelor:::MNN_auto_order(batches, restrict=restrict)
        ref <- batchelor:::MNN_supplied_order(batches, restrict=restrict, ordering=ordering)
        compare_common(store, ref)

        for (j in 2:length(batches)) {
            k <- 20 + j

            store <- batchelor:::.advance(store, k=k)
            ref <- batchelor:::.advance(ref, k=k)
            compare_common(store, ref)
            
            expect_identical(head(ordering, j-1), batchelor:::.get_reference_indices(store))
            current <- batchelor:::.get_current_index(store)
            expect_identical(ordering[j], current)

            store <- batchelor:::.compile(store, batches[[current]])
            ref <- batchelor:::.compile(ref, batches[[current]])
            compare_common(store, ref)

            expect_identical(head(ordering, j), batchelor:::.get_reference_indices(store))
        }
    }
})

#####################################################

set.seed(123005)
test_that("replacement of the reference works in MNN store", {
    batches <- GENERATOR(c(1000, 2000, 500), 10)

    store <- batchelor:::MNN_supplied_order(batches)
    store <- batchelor:::.advance(store, 10)

    running.ref <- batches[[1]] 
    running.ref <- running.ref + rnorm(length(running.ref))
    store <- batchelor:::.compile(store, batches[[2]], new.reference=running.ref)
    running.ref <- rbind(running.ref, batches[[2]])
    expect_identical(batchelor:::.get_reference(store), running.ref)
    
    store <- batchelor:::.advance(store, 10)

    running.ref <- batchelor:::.get_reference(store)
    running.ref <- running.ref + rnorm(length(running.ref))
    store <- batchelor:::.compile(store, batches[[3]], new.reference=running.ref)
    running.ref <- rbind(running.ref, batches[[3]])
    expect_identical(batchelor:::.get_reference(store), running.ref)
})

set.seed(123006)
test_that("error conditions trigger successfully", {
    batches <- GENERATOR(c(1000, 2000, 500), 10)
    expect_error(batchelor:::MNN_supplied_order(list(batches[[1]][,1:5], batches[[2]])), "number of columns must be the same")

    ref <- batchelor:::MNN_supplied_order(batches)
    expect_error(batchelor:::.compile(ref), "not available yet")
    ref <- batchelor:::.advance(ref, 10)
    expect_error(batchelor:::.advance(ref, 20), "not been assimilated yet")
    
    expect_error(ref <- batchelor:::.compile(ref, batches[[2]][,1:5]), "invalid dimensions for the corrected batch")
    expect_error(ref <- batchelor:::.compile(ref, batches[[3]]), "invalid dimensions for the corrected batch")
    expect_error(ref <- batchelor:::.compile(ref, batches[[2]], new.reference=batches[[3]]), "invalid dimensions for the updated reference")
    ref <- batchelor:::.compile(ref, batches[[2]])
    expect_error(ref <- batchelor:::.compile(ref, batches[[2]]), "not available yet")
})
