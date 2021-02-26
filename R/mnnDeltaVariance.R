#' Computes the variance of the paired MNN deltas
#'
#' @inheritParams fastMNN
#' @param pairs A \linkS4class{DataFrame} or list of \linkS4class{DataFrame}s containing MNN pairing information.
#' Each row of each DataFrame specifies an MNN pair; each DataFrame should have two columns containing the column indices of paired cells.
#' This is typically produced by \code{\link{fastMNN}}, see that documentation for more information.
#' 
#' @details
#' The \dQuote{MNN delta} is defined as the difference in the \emph{uncorrected} expression values of MNN-paired cells.
#' We compute the variance of these values across all pairs for each gene; genes with highly variable deltas are strongly affected by the correction beyond a simple shift.
#' By looking through the top genes, we may conclude that the correction did something unusual if, e.g., we find important markers in the top set.
#' Conversely, if key genes do not have unusually high variances in their delta, we may gain some confidence in the reliability of the correction -
#' at least with respect to the behavior of those genes.
#'
#' We compute the variance of the deltas rather than using their magnitude, 
#' as a shift in expression space is not particularly concerning (and is, in fact, par for the course for batch correction).
#' Rather, we are interested in identifying non-linear changes that might be indicative of correction errors, e.g., inappropriate merging of two different subpopulations.
#' This would manifest as high variances for the relevant marker genes.
#' Of course, whether or not this is actually an error can be a subjective decision - loss of some biological heterogeneity is often an acceptable cost for flexible correction.
#'
#' @return
#' A numeric vector of variances for each row in the input objects (or as specified by \code{subset.row}).
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom DelayedArray blockApply rowAutoGrid
mnnDeltaVariance <- function(..., pairs, subset.row=NULL, BPPARAM=SerialParam()) {
    batches <- .unpackLists(...)
    checkBatchConsistency(batches, cells.in.columns=TRUE)

    # Extracting information from SCEs.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }
    
    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    if (!is(pairs, "DataFrame")) {
        first <- unlist(lapply(pairs, function(x) x[,1]))
        second <- unlist(lapply(pairs, function(x) x[,2]))
    } else {
        first <- pairs[,1]
        second <- pairs[,2]
    }

    X <- blockApply(x, FUN=.compute_mnn_variance, first=first, second=second, grid=rowAutoGrid(x), BPPARAM=BPPARAM)
    unlist(X)    
}

#' @importFrom methods is as 
#' @importFrom DelayedMatrixStats rowVars
#' @importClassesFrom DelayedArray SparseArraySeed
#' @importClassesFrom Matrix sparseMatrix
.compute_mnn_variance <- function(block, pairs) {
    if (is(block, "SparseArraySeed")) {
        block <- as(block, "sparseMatrix")
    }
    b1 <- block[,pairs$left,drop=FALSE]
    b2 <- block[,pairs$right,drop=FALSE]
    rowVars(b1 - b2)
}
