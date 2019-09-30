#' Regress out batch effects
#'
#' Fit a linear model to regress out uninteresting factors of variation.
#'
#' @param ... Two or more log-expression matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a count matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for gene expression matrix inputs.
#'
#' If multiple objects are supplied, each object is assumed to contain all and only cells from a single batch.
#' Objects of different types can be mixed together. 
#' If a single object is supplied, \code{batch} should also be specified.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param restrict A list of length equal to the number of objects in \code{...}.
#' Each entry of the list corresponds to one batch and specifies the cells to use when computing the correction.
#' @param subset.row A vector specifying which features to use for correction.
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values, if SingleCellExperiment objects are present in \code{...}.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} object containing the \code{corrected} assay.
#' This contains corrected log-expression values for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
#'
#' Cells in the output object are always ordered in the same manner as supplied in \code{...}.
#' For a single input object, cells will be reported in the same order as they are arranged in that object.
#' In cases with multiple input objects, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#' 
#' @details
#' This function fits a linear model to the log-expression values for each gene and returns the residuals.
#' The model is parameterized as a one-way layout with the batch of origin,
#' so the residuals represent the expression values after correcting for the batch effect.
#'
#' The novelty of this function is that it returns a \linkS4class{ResidualMatrix} in as the \code{"corrected"} assay.
#' This avoids explicitly computing the residuals, which would result in a loss of sparsity or similar problems.
#' Rather, the residuals are either computed as needed or are never explicitly computed as all (e.g., during matrix multiplication).
#'
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs, though this is purely for convenience as each gene is processed independently of other genes.
#'
#' See \code{?"\link{batchelor-restrict}"} for a description of the \code{restrict} argument.
#' Specifically, this function will compute the model coefficients using only the specified subset of cells.
#' The regression will then be applied to all cells in each batch.
#'
#' @author Aaron Lun
#'
#' @examples
#' means <- 2^rgamma(1000, 2, 1)
#' A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
#' A2 <- matrix(rpois(10000, lambda=means*runif(1000, 0, 2)), ncol=50) # Batch 2
#'
#' B1 <- log2(A1 + 1)
#' B2 <- log2(A2 + 1)
#' out <- regressBatches(B1, B2) 
#'
#' @seealso
#' \code{\link{rescaleBatches}}, for another approach to regressing out the batch effect.
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics cbind
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix t 
#' @importFrom BiocSingular ResidualMatrix
#' @importFrom stats model.matrix
#' @importFrom S4Vectors DataFrame
#' @importFrom utils head
#' @importFrom DelayedArray seed seed<-
regressBatches <- function(..., batch=NULL, restrict=NULL, subset.row=NULL, assay.type="logcounts") {
    batches <- list(...)
    checkBatchConsistency(batches)
    restrict <- checkRestrictions(batches, restrict)

    # Pulling out information from the SCE objects.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        sce.batches <- batches[is.sce]
        checkSpikeConsistency(sce.batches)
        subset.row <- .SCE_subset_genes(subset.row, sce.batches[[1]], FALSE)
        batches[is.sce] <- lapply(sce.batches, assay, i=assay.type)
    }

    # Creating a single matrix object.
    if (length(batches) > 1L) {
        combined <- do.call(cbind, batches)
        ncells <- vapply(batches, ncol, 0L)
        batch <- rep(seq_along(batches), ncells)
        if (!is.null(nm <- names(batches))) {
            if (anyDuplicated(nm)) {
                stop("names of batches should be unique")
            }
            batch <- nm[batch]
        }
        if (!is.null(restrict)) {
            restrict <- unlist(mapply("+", restrict, c(0L, head(ncells, -1L)), SIMPLIFY=FALSE))
        }
    } else if (length(batches)==1L) {
        combined <- batches[[1]]
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        restrict <- restrict[[1]]
    } else {
        stop("at least two batches must be specified")
    }

    if (!is.null(subset.row)) {
        combined <- combined[subset.row,,drop=FALSE]
    }

    design <- model.matrix(~0 + factor(batch))
    combined <- t(combined)
    corrected <- ResidualMatrix(combined, design)

    if (!is.null(restrict)) {
        subdesign <- design[restrict,,drop=FALSE]
        QR <- qr(subdesign)
        Q <- as.matrix(qr.Q(QR))
        Qty <- as.matrix(crossprod(Q, combined[restrict,,drop=FALSE]))

        # TODO: modify ResidualMatrix() to accommodate non-standard left %*% right.
        coefs <- backsolve(qr.R(QR), Qty)
        seed(corrected)@Qty <- coefs
        seed(corrected)@Q <- as.matrix(design)
        seed(corrected)@centered <- FALSE # no guarantee of centering if you add extra cells.
    }

    SingleCellExperiment(list(corrected=t(corrected)), colData=DataFrame(batch=batch)) 
}
