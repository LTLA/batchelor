#' No correction
#'
#' Provides a no-correction method that has the same interface as the correction functions.
#' This allows users to easily swap function calls to examine the effect of correction.
#'
#' @inheritParams fastMNN
#'
#' @return
#' A \linkS4class{SingleCellExperiment} is returned where each row is a gene and each column is a cell. 
#' This contains:
#' \itemize{
#' \item A \code{merged} matrix in the \code{assays} slot, containing the merged expression values from all elements of \code{...}.
#' \item A \code{batch} column in the \code{colData} slot, containing the batch of origin for each row (i.e., cell) in \code{corrected}.
#' }
#'
#' @details
#' This function is effectively equivalent to \code{cbind}ing the matrices together without any correction.
#' The aim is to provide a consistent interface that allows users to simply combine batches without additional operations. 
#' This is often desirable as a negative control to see if the transformation is actually beneficial.
#' It also allows for convenient downstream analyses that are based on the uncorrected data, e.g., differential expression.
#'
#' @author Aaron Lun
#'
#' @examples
#' B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
#' out <- noCorrect(B1, B2)
#'
#' # Same as combining the expression values.
#' stopifnot(all(assay(out)==cbind(B1, B2)))
#'
#' # Specifies which cell came from which batch:
#' str(out$batch)
#' @export
#' @importFrom BiocGenerics cbind
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
noCorrect <- function(..., batch=NULL, subset.row=NULL, assay.type="logcounts") {
    batches <- .unpack_batches(...)
    checkBatchConsistency(batches, cells.in.columns=TRUE)

    # Extracting information from SCEs.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }

    if (!is.null(subset.row)) {
        for (i in seq_along(batches)) {
            batches[[i]] <- batches[[i]][subset.row,,drop=FALSE]
        }
    }

    if (length(batches)==1L) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        output <- batches[[1]]
    } else {
        output <- do.call(cbind, batches)
        batch <- rep(seq_along(batches), vapply(batches, ncol, 0L))
        if (!is.null(names(batches))) {
            if (anyDuplicated(names(batches))) {
                stop("names of batches should be unique")
            }
            batch <- names(batches)[batch]
        }
    }

    SingleCellExperiment(list(merged=output), colData=DataFrame(batch=batch))
}
