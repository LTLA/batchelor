#' Scale counts across batches
#'
#' Scale counts so that the average count within each batch is the same for each gene.
#'
#' @param ... Two or more log-expression matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
#' Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a count matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for matrix inputs.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param log.base A numeric scalar specifying the base of the log-transformation.
#' @param pseudo.count A numeric scalar specifying the pseudo-count used for the log-transformation.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} object containing the \code{corrected} assay.
#' This contains corrected log-expression values for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
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
rescaleBatches <- function(..., batch=NULL, log.base=2, pseudo.count=1, subset.row=NULL, assay.type="logcounts", get.spikes=FALSE) {
    batches <- list(...)

    # Pulling out information from the SCE objects.        
    if (.check_if_SCEs(batches)) {
        .check_batch_consistency(batches, byrow=TRUE)
        .check_spike_consistency(batches)
        subset.row <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes)
        batches <- lapply(batches, assay, i=assay.type, withDimnames=FALSE)
    } else {
        .check_batch_consistency(batches, byrow=TRUE)
    }

    # Subsetting by 'batch'. 
    do.split <- length(batches)==1L
    if (do.split) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }

        divided <- .divide_into_batches(batches[[1]], batch=batch, byrow=FALSE)
        batches <- divided$batches
    } 

    output <- do.call(.rescale_batches, c(batches, list(log.base=log.base, pseudo.count=pseudo.count, subset.row=subset.row)))

    # Reordering the output for correctness.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
    }

    return(output)
}

############################################
# Internal main function, to separate the input handling from the actual calculations.

#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom BiocGenerics cbind rowMeans
.rescale_batches <- function(..., log.base=2, pseudo.count=1, subset.row=NULL) {
    batches <- list(...)
    nbatches <- length(batches)
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }

    # Computing the unlogged means for each matrix.
    subset.row.m1 <- .row_subset_to_index(batches[[1]], subset.row) - 1L
    averages <- vector("list", nbatches)
    for (b in seq_along(batches)) {
        averages[[b]] <- .Call(cxx_unlog_exprs_mean, batches[[b]], log.base, pseudo.count, subset.row.m1)
    }

    # Defining the reference.
    reference <- do.call(pmin, averages)
    for (b in seq_along(batches)) {
        rescale <- reference / averages[[b]] 
        rescale[!is.finite(rescale)] <- 0
        batches[[b]] <- .Call(cxx_unlog_exprs_scaled, batches[[b]], log.base, pseudo.count, subset.row.m1, rescale)
    }

    ncells.per.batch <- vapply(batches, ncol, FUN.VALUE=0L)
    batch.names <- .create_batch_names(names(batches), ncells.per.batch)
    SingleCellExperiment(list(corrected=do.call(cbind, batches)), 
        colData=DataFrame(batch=batch.names$ids)) 
}
