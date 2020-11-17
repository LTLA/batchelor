#' Scale counts across batches
#'
#' Scale counts so that the average count within each batch is the same for each gene.
#'
#' @inheritParams fastMNN
#' @param log.base A numeric scalar specifying the base of the log-transformation.
#' @param pseudo.count A numeric scalar specifying the pseudo-count used for the log-transformation.
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
#' This function assumes that the log-expression values were computed by a log-transformation of normalized count data, plus a pseudo-count.
#' It reverses the log-transformation and scales the underlying counts in each batch so that the average (normalized) count is equal across batches.
#' The assumption here is that each batch contains the same population composition.
#' Thus, any scaling difference between batches is technical and must be removed.
#'
#' This function is approximately equivalent to centering in log-expression space, the simplest application of linear regression methods for batch correction.
#' However, by scaling the raw counts, it avoids loss of sparsity that would otherwise result from centering.
#' It also mitigates issues with artificial differences in variance due to log-transformation.
#' This is done by always downscaling to the lowest average expression for each gene such that differences in variance are dampened by the addition of the pseudo-count.
#'
#' Use of \code{rescaleBatches} assumes that the uninteresting factors described in \code{design} are orthogonal to the interesting factors of variation.
#' For example, each batch is assumed to have the same composition of cell types.
#' If this is not true, the correction will not only be incomplete but may introduce spurious differences.
#'
#' The output values are always re-log-transformed with the same \code{log.base} and \code{pseudo.count}.
#' These can be used directly in place of the input values for downstream operations.
#'
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs, though this is purely for convenience as each gene is processed independently of other genes.
#'
#' See \code{?"\link{batchelor-restrict}"} for a description of the \code{restrict} argument.
#' Specifically, the function will compute the scaling differences using only the specified subset of cells, and then apply the re-scaling to all cells in each batch.
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
#' out <- rescaleBatches(B1, B2) 
#' 
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom scuttle .unpackLists
rescaleBatches <- function(..., batch=NULL, restrict=NULL, log.base=2, pseudo.count=1, subset.row=NULL, correct.all=FALSE, assay.type="logcounts") {
    originals <- batches <- .unpackLists(...)
    checkBatchConsistency(batches)
    restrict <- checkRestrictions(batches, restrict)

    # Pulling out information from the SCE objects.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type, withDimnames=FALSE)
    }

    # Subsetting by 'batch'. 
    do.split <- length(batches)==1L
    if (do.split) {
        divided <- divideIntoBatches(batches[[1]], batch=batch, byrow=FALSE, restrict=restrict[[1]])
        batches <- divided$batches
        restrict <- divided$restrict
    } 

    if (correct.all) {
        subset.row <- NULL
    }
    output <- do.call(.rescale_batches, c(batches, list(log.base=log.base, pseudo.count=pseudo.count, subset.row=subset.row, restrict=restrict)))

    # Reordering the output for correctness.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
    }

    # Adding dimension names.
    .rename_output(output, originals, subset.row=subset.row)
}

############################################
# Internal main function, to separate the input handling from the actual calculations.

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BiocGenerics cbind 
#' @importFrom Matrix rowMeans
.rescale_batches <- function(..., log.base=2, pseudo.count=1, subset.row=NULL, restrict=NULL) {
    batches <- list(...)
    nbatches <- length(batches)
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }

    # Computing the unlogged means for each matrix, using only the restricted subset of cells.
    for (b in seq_along(batches)) {
        if (!is.null(subset.row)) {
            batches[[b]] <- batches[[b]][subset.row,,drop=FALSE]
        }
        batches[[b]] <- .unlog(batches[[b]], log.base=log.base, pseudo.count=pseudo.count)
    }

    averages <- vector("list", nbatches)
    for (b in seq_along(batches)) {
        curbatch <- batches[[b]]
        if (!is.null(currestrict <- restrict[[b]])) {
            curbatch <- curbatch[,currestrict,drop=FALSE]
        }
        averages[[b]] <- rowMeans(curbatch)
    }

    # Defining the reference.
    reference <- do.call(pmin, averages)
    for (b in seq_along(batches)) {
        rescale <- reference / averages[[b]] 
        rescale[!is.finite(rescale)] <- 0
        batches[[b]] <- .relog(batches[[b]] * rescale, log.base=log.base, pseudo.count=pseudo.count)
    }

    batch.labels <- names(batches)
    if (is.null(batch.labels)) {
        batch.labels <- seq_along(batches)
    } else if (anyDuplicated(batch.labels)) {
        stop("names of batches should be unique")
    }
    ncells.per.batch <- vapply(batches, ncol, FUN.VALUE=0L)
    batch.names <- rep(batch.labels, ncells.per.batch)

    SingleCellExperiment(list(corrected=do.call(cbind, batches)), 
        colData=.create_unnamed_coldata(batch.names)) 
}

############################################
# A helpful dispatcher system to accommodate different matrix representations.

setGeneric(".unlog", function(x, log.base, pseudo.count) standardGeneric(".unlog"))

setMethod(".unlog", "ANY", function(x, log.base, pseudo.count) {
    log.base^x - pseudo.count
})

#' @importClassesFrom Matrix dsparseMatrix
setMethod(".unlog", "dsparseMatrix", function(x, log.base, pseudo.count) {
    if (pseudo.count!=1) {
        callNextMethod()
    } else {
        x@x <- log.base^x@x - pseudo.count
        x
    }
})

setGeneric(".relog", function(x, log.base, pseudo.count) standardGeneric(".relog"))

setMethod(".relog", "ANY", function(x, log.base, pseudo.count) {
    log(x + pseudo.count, log.base)
})

#' @importClassesFrom Matrix dsparseMatrix
setMethod(".relog", "dsparseMatrix", function(x, log.base, pseudo.count) {
    if (pseudo.count!=1) {
        callNextMethod()
    } else {
        x@x <- log(x@x + pseudo.count, log.base)
        x
    }
})
