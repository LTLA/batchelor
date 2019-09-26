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
