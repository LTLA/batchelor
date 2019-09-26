# Tree binarization.
# library(batchelor); library(testthat); source("test-tree.R")

test_that("tree binarization works correctly", {
    expect_identical(
        batchelor:::.binarize_tree(list(1,2,3)),
        list(list(1,2),3)
    )

    expect_identical(
        batchelor:::.binarize_tree(list(1,2,3,4,5)),
        list(list(list(list(1,2),3),4),5)
    )

    # Handles useless internal nodes.
    expect_identical(
        batchelor:::.binarize_tree(list(list(1,2,3))),
        list(list(1,2),3)
    )

    expect_identical(
        batchelor:::.binarize_tree(list(list(1),list(2))),
        list(1,2)
    )

    # A more complex example.
    expect_identical(
        batchelor:::.binarize_tree(list(list(1,2,3), list(4,5,6))),
        list(list(list(1,2),3), list(list(4,5),6))
    )

    # No change if it's already binary.
    ref <- list(list(list(1,2),list(3,4)), list(list(5,6), list(7,8)))
    expect_identical(ref, batchelor:::.binarize_tree(ref))

    expect_error(batchelor:::.binarize_tree(list(list(), list(1,2,3), list(4,5,6))), "node with no children")
})

test_that("asdasd", {
    B1 <- matrix(1)
    B2 <- matrix(2)
    B3 <- matrix(3)
    B4 <- matrix(4)

    out <- batchelor:::.create_tree_predefined(list(B1, B2, B3), restrict=1:3*10L, merge.order=1:3)
    expect_identical(B1, out[[1]][[1]]@data)
    expect_identical(10L, out[[1]][[1]]@restrict)
    expect_identical(B2, out[[1]][[2]]@data)
    expect_identical(20L, out[[1]][[2]]@restrict)
    expect_identical(B3, out[[2]]@data)
    expect_identical(30L, out[[2]]@restrict)

    out <- batchelor:::.create_tree_predefined(list(B1, B2, B3), restrict=1:3*10L, merge.order=3:1)
    expect_identical(B3, out[[1]][[1]]@data)
    expect_identical(30L, out[[1]][[1]]@restrict)
    expect_identical(B2, out[[1]][[2]]@data)
    expect_identical(20L, out[[1]][[2]]@restrict)
    expect_identical(B1, out[[2]]@data)
    expect_identical(10L, out[[2]]@restrict)

    out <- batchelor:::.create_tree_predefined(list(B1, B2, B3, B4), restrict=NULL, merge.order=list(list(1,4), list(3,2)))
    expect_identical(B1, out[[1]][[1]]@data)
    expect_identical(NULL, out[[1]][[1]]@restrict)
    expect_identical(B4, out[[1]][[2]]@data)
    expect_identical(NULL, out[[1]][[2]]@restrict)
    expect_identical(B3, out[[2]][[1]]@data)
    expect_identical(NULL, out[[2]][[1]]@restrict)
    expect_identical(B2, out[[2]][[2]]@data)
    expect_identical(NULL, out[[2]][[2]]@restrict)

    # Converts from character properly.
    expect_identical(
        batchelor:::.create_tree_predefined(list(B1, B2, B3), restrict=1:3*10L, merge.order=1:3),
        batchelor:::.create_tree_predefined(list(A=B1, B=B2, C=B3), restrict=1:3*10L, merge.order=LETTERS[1:3])
    )

    # Various failure modes for the merge.order specification.
    expect_error(batchelor:::.create_tree_predefined(list(B1, B2), merge.order=1:3), "invalid leaf nodes")
    expect_error(batchelor:::.create_tree_predefined(list(B1, B2), merge.order=c(1,1)), "invalid leaf nodes")
    expect_error(batchelor:::.create_tree_predefined(list(A=B1, B=B2), merge.order=c("A","C")), "invalid leaf nodes")
})
