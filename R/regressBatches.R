#' Regress out batch effects
#'
#' Fit a linear model to each gene regress out uninteresting factors of variation, returning a matrix of residuals.
#'
#' @inheritParams fastMNN
#' @param design A numeric design matrix with number of rows equal to the total number of cells,
#' specifying the experimental factors to remove.
#' Each row corresponds to a cell in the order supplied in \code{...}.
#' @param keep Integer vector specifying the coefficients of \code{design} to \emph{not} regress out,
#' see the \code{\link{ResidualMatrix}} constructor for more details.
#' @param d Numeric scalar specifying the number of dimensions to use for PCA via \code{\link{multiBatchPCA}}.
#' If \code{NA}, no PCA is performed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA should be parallelized.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} object containing the \code{corrected} assay.
#' This contains the computed residuals for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
#'
#' Cells in the output object are always ordered in the same manner as supplied in \code{...}.
#' For a single input object, cells will be reported in the same order as they are arranged in that object.
#' In cases with multiple input objects, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#'
#' If \code{d} is not \code{NA}, a PCA is performed on the residual matrix via \code{\link{multiBatchPCA}},
#' and an additional \code{corrected} field is present in the \code{\link{reducedDims}} of the output object.
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
#' Like \code{\link{rescaleBatches}}, this function assumes that the batch effect is orthogonal to the interesting factors of variation.
#' For example, each batch is assumed to have the same composition of cell types.
#' The same reasoning applies to any uninteresting factors specified in \code{design}, including continuous variables.
#' For example, if one were to use this function to regress out cell cycle, the assumption is that all cell types are similarly distributed across cell cycle phases.
#' If this is not true, the correction will not only be incomplete but can introduce spurious differences.
#' 
#' See \code{?"\link{batchelor-restrict}"} for a description of the \code{restrict} argument.
#' Specifically, this function will compute the model coefficients using only the specified subset of cells.
#' The regression will then be applied to all cells in each batch.
#'
#' If set, the \code{d} option will perform a PCA via \code{\link{multiBatchPCA}}.
#' This is provided for convenience as efficiently executing a PCA on a \linkS4class{ResidualMatrix} is not always intuitive.
#' (Specifically, \linkS4class{BiocSingularParam} objects must be set up with \code{deferred=TRUE} for best performance.)
#' The arguments \code{BSPARAM}, \code{deferred} and \code{BPPARAM} only have an effect when \code{d} is set to a non-\code{NA} value.
#'
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' If a subset of genes is specified, residuals are only returned for that subset.
#' Similarly, if \code{d} is set, only the genes in the subset are used to perform the PCA.
#' If additionally \code{correct.all=TRUE}, residuals are returned for all genes but only the subset is used for the PCA.
#'
#' If \code{as.altexp} is specified, the corresponding \code{\link{altExp}} is used for all calculations.
#' The result is the same as if the alternative Experiments were directly passed to this function; subsetting, PCA and so on are applied to data from the alternative Experiments.
#' For studies with multiple feature sets, it may be desirable to run this function repeatedly with different \code{as.altexp} and then collate the objects into another SingleCellExperiment.
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
#' The \linkS4class{ResidualMatrix} class, for the class of the residual matrix.
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics cbind
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim<- altExp
#' @importFrom Matrix t 
#' @importFrom ResidualMatrix ResidualMatrix
#' @importFrom stats model.matrix
#' @importFrom S4Vectors DataFrame
#' @importFrom utils head
#' @importFrom DelayedArray seed seed<-
#' @importFrom scuttle .unpackLists
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular IrlbaParam 
regressBatches <- function(..., batch=NULL, design=NULL, keep=NULL, restrict=NULL, 
    subset.row=NULL, correct.all=FALSE, d=NA, deferred=TRUE, 
    assay.type="logcounts", as.altexp=NULL, 
    BSPARAM=IrlbaParam(), BPPARAM=SerialParam())
{
    batches <- .unpackLists(...)
    if (!is.null(as.altexp)) {
        batches <- lapply(batches, altExp, e=as.altexp)
    }
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

        if (is.null(batch)) {
            if (!is.null(design)) {
                batch <- rep(1L, ncol(combined))
            } else {
                stop("'batch' must be specified if '...' has only one object")
            }
        }

        restrict <- restrict[[1]]
    } else {
        stop("at least two batches must be specified")
    }

    if (!correct.all && !is.null(subset.row)) {
        combined <- combined[subset.row,,drop=FALSE]
        subset.row <- NULL
    }

    if (is.null(design)) {
        design <- model.matrix(~0 + factor(batch))
    } else if (nrow(design)!=ncol(combined)) {
        stop("'nrow(design)' should be equal to the total number of cells")
    }

    combined <- t(combined)
    corrected <- ResidualMatrix(combined, design, keep=keep, restrict=restrict)
    sce <- SingleCellExperiment(list(corrected=t(corrected)), colData=.create_unnamed_coldata(batch)) 

    if (!is.na(d)) {
        BSPARAM <- .set_deferred(BSPARAM, deferred)
        pc.out <- .multi_pca_single(assay(sce), batch, subset.row=subset.row, d=d, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
        reducedDim(sce, "corrected") <- pc.out[[1]]
    }

    sce
}
