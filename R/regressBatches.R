#' Regress out batch effects
#'
#' Fit a linear model to regress out uninteresting factors of variation.
#'
#' @inheritParams fastMNN
#' @param design A numeric design matrix with number of rows equal to the total number of cells,
#' specifying the experimental factors to remove.
#' Each row corresponds to a cell in the order supplied in \code{...}.
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
#' By default, the model is parameterized as a one-way layout with the batch of origin,
#' so the residuals represent the expression values after correcting for the batch effect.
#' The novelty of this function is that it returns a \linkS4class{ResidualMatrix} in as the \code{"corrected"} assay.
#' This avoids explicitly computing the residuals, which would result in a loss of sparsity or similar problems.
#' Rather, residuals are either computed as needed or are never explicitly computed at all (e.g., during matrix multiplication).
#' This means that \code{regressBatches} is faster and lighter than naive regression or even \code{\link{rescaleBatches}}.
#' 
#' More complex designs should be explicitly specified with the \code{design} argument, e.g., to regress out a covariate.
#' This can be any full-column-rank matrix that is typically constructed with \code{\link{model.matrix}}.
#' If \code{design} is specified with a single object in \code{...}, \code{batch} is ignored.
#' If \code{design} is specified with multiple objects, regression is applied to the matrix obtained by \code{cbind}ing all of those objects together; this means that the first few rows of \code{design} correspond to the cells from the first object, then the next rows correspond to the second object and so on.
#' 
#' Like \code{\link{rescaleBatches}}, this function assumes that the uninteresting factors described in \code{design} are orthogonal to the interesting factors of variation.
#' For example, each batch is assumed to have the same composition of cell types.
#' In the continuous case (e.g., regression of cell cycle), the assumption is that all cell types are similarly distributed across cell cycle phases.
#' If this is not true, the correction will not only be incomplete but can introduce spurious differences.
#' 
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs, though this is purely for convenience as each gene is processed independently of other genes.
#' Indeed, setting \code{correct.all=TRUE} is equivalent to forcing \code{subset.row=NULL},
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
regressBatches <- function(..., batch=NULL, design=NULL, restrict=NULL, 
    subset.row=NULL, correct.all=FALSE, assay.type="logcounts") 
{
    batches <- .unpack_batches(...)
    checkBatchConsistency(batches)
    restrict <- checkRestrictions(batches, restrict)

    # Pulling out information from the SCE objects.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
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

        # We just need a placeholder for the colData setting below.
        if (!is.null(design)) {
            batch <- matrix(0L, ncol(combined), 0) 
        }

        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        restrict <- restrict[[1]]
    } else {
        stop("at least two batches must be specified")
    }

    if (!correct.all && !is.null(subset.row)) {
        combined <- combined[subset.row,,drop=FALSE]
    }

    if (is.null(design)) {
        design <- model.matrix(~0 + factor(batch))
    } else if (nrow(design)!=ncol(combined)) {
        stop("'nrow(design)' should be equal to the total number of cells")
    }

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
