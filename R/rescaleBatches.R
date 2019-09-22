#' Scale counts across batches
#'
#' Scale counts so that the average count within each batch is the same for each gene.
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
#' @param log.base A numeric scalar specifying the base of the log-transformation.
#' @param pseudo.count A numeric scalar specifying the pseudo-count used for the log-transformation.
#' @param subset.row A vector specifying which features to use for correction.
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes Deprecated, a logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
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
#' This function is equivalent to centering in log-expression space, the simplest application of linear regression methods for batch correction.
#' However, by scaling the raw counts, it avoids loss of sparsity that would otherwise result from centering.
#' It also mitigates issues with artificial differences in variance due to log-transformation.
#'
#' The output values are always re-log-transformed with the same \code{log.base} and \code{pseudo.count}.
#' These can be used directly in place of the input values for downstream operations.
#'
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs, though this is purely for convenience as each gene is processed independently of other genes.
#'
#' @section Using restriction:
#' It is possible to compute the correction using only a subset of cells in each batch, and then extrapolate that correction to all other cells.
#' This may be desirable in experimental designs where a control set of cells from the same source population were run on different batches.
#' Any difference in the controls must be artificial in origin and can be directly removed without making further biological assumptions.
#'
#' To do this, users should set \code{restrict} to specify the subset of cells in each batch to be used for correction.
#' This should be set to a list of length equal to the length of \code{...}, where each element is a subsetting vector to be applied to the columns of the corresponding batch.
#' A \code{NULL} element indicates that all the cells from a batch should be used.
#' In situations where one input object contains multiple batches, \code{restrict} is simply a list containing a single subsetting vector for that object.
#'
#' The function will compute the scaling differences using only the specified subset of cells.
#' However, the re-scaling will then be applied to all cells in each batch - hence the extrapolation.
#' This means that the output is always of the same dimensionality, regardless of whether \code{restrict} is specified.
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
rescaleBatches <- function(..., batch=NULL, restrict=NULL, log.base=2, pseudo.count=1, subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {
    originals <- batches <- list(...)
    checkBatchConsistency(batches)
    restrict <- checkRestrictions(batches, restrict)

    # Pulling out information from the SCE objects.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        sce.batches <- batches[is.sce]
        checkSpikeConsistency(sce.batches)
        subset.row <- .SCE_subset_genes(subset.row, sce.batches[[1]], get.spikes)
        batches[is.sce] <- lapply(sce.batches, assay, i=assay.type, withDimnames=FALSE)
    }

    # Subsetting by 'batch'. 
    do.split <- length(batches)==1L
    if (do.split) {
        divided <- divideIntoBatches(batches[[1]], batch=batch, byrow=FALSE, restrict=restrict[[1]])
        batches <- divided$batches
        restrict <- divided$restrict
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
    subset.row.m1 <- .row_subset_to_index(batches[[1]], subset.row) - 1L
    averages <- vector("list", nbatches)

    for (b in seq_along(batches)) {
        curbatch <- batches[[b]]

        currestrict <- restrict[[b]]
        if (!is.null(currestrict)) {
            curbatch <- curbatch[,currestrict,drop=FALSE]
        }

        averages[[b]] <- .Call(cxx_unlog_exprs_mean, curbatch, log.base, pseudo.count, subset.row.m1)
    }

    # Defining the reference.
    reference <- do.call(pmin, averages)
    for (b in seq_along(batches)) {
        rescale <- reference / averages[[b]] 
        rescale[!is.finite(rescale)] <- 0
        batches[[b]] <- .Call(cxx_unlog_exprs_scaled, batches[[b]], log.base, pseudo.count, subset.row.m1, rescale)
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
        colData=DataFrame(batch=batch.names)) 
}
