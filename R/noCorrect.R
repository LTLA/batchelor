#' No correction
#'
#' Provides a no-correction method that has the same interface as the correction functions.
#' This allows users to easily swap function calls to examine the effect of correction.
#' 
#' @param ... One or more log-expression matrices where genes correspond to rows and cells correspond to columns, if \code{pc.input=FALSE}.
#' Each matrix should contain the same number of rows, corresponding to the same genes in the same order.
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a log-expression matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for gene expression matrix inputs.
#'
#' In all cases, each object contains cells from a single batch; multiple objects represent separate batches of cells.
#' Objects of different types can be mixed together. 
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param subset.row A vector specifying which features to retain.
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values.
#' Only used for SingleCellExperiment inputs. 
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
    batches <- list(...)
    is.sce <- checkIfSCE(batches)

    # Extracting information from SCEs.
    if (any(is.sce)) {
        sce.batches <- batches[is.sce]
        # NOTE: do NOT set withDimnames=FALSE, this will break the consistency check.
        sce.batches <- lapply(sce.batches, assay, i=assay.type)
        batches[is.sce] <- sce.batches
    }

    # Batch consistency checks.
    checkBatchConsistency(batches, cells.in.columns=TRUE)
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
