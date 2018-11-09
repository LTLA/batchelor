#' Scale counts across batches
#'
#' Scale counts so that the average count within each batch is the same for each gene.
#'
#' @param ... Two or more count matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
#' Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a count matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for matrix inputs.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param use.size.factors A logical scalar specifying if size factors should be used when computing the batch-specific average with \code{\link{calculateAverage}}.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#'
#' @return
#' A \linkS4class{SummarizedExperiment} object containing the \code{corrected} assay.
#' This contains scaled count values for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
#' 
#' @details
#' This function will downscale counts in each batch so that the average count is equal across batches.
#' The assumption here is that each batch contains the same population composition.
#' Thus, any scaling difference between batches is technical and must be removed.
#'
#' This function is equivalent to centering in log-expression space, the simplest application of linear regression methods for batch correction.
#' However, by scaling the raw counts, it avoids loss of sparsity that would otherwise result from centering.
#' It also mitigates issues with artificial differences in variance due to log-transformation.
#'
#' For exploratory analyses, the corrected values can be treated as being equivalent to the raw counts.
#' Users can use them to compute size factors and log-transformed expression values prior to downstream procedures like dimensionality reduction and clustering.
#' However, the corrected values do not preserve the mean-variance relationship expected by count-based models, and should not be used in \pkg{edgeR} and \pkg{DESeq2}.
#'
#' @author Aaron Lun
#'
#' @examples
#' means <- 2^rgamma(1000, 2, 1)
#' B1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
#' B2 <- matrix(rpois(10000, lambda=means*runif(1000, 0, 2)), ncol=50) # Batch 2
#' out <- rescaleBatches(B1, B2) 
#' 
#' @export
#' @importFrom scater calculateAverage
#' @importFrom SummarizedExperiment assay
rescaleBatches <- function(..., batch=NULL, use.size.factors=TRUE, subset.row=NULL, assay.type="counts", get.spikes=FALSE) {
    batches <- list(...)
    
    # Pulling out information from the SCE objects.        
    if (.check_if_SCEs(batches)) {
        .check_batch_consistency(batches, byrow=TRUE)
        .check_spike_consistency(batches)
        subset.row <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes)
        averages <- lapply(batches, calculateAverage, exprs_values=assay.type, subset_row=subset.row, use_size_factors=use.size.factors)
        batches <- lapply(batches, assay, i=assay.type, withDimnames=FALSE)
    } else {
        .check_batch_consistency(batches, byrow=TRUE)
        averages <- lapply(batches, calculateAverage, subset_row=subset.row)
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

    output <- do.call(.rescale_batches, list(batches=batches, subset.row=subset.row, averages=averages))

    # Reordering the output for correctness.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
    }

    return(output)
}

############################################
# Internal main function, to separate the input handling from the actual calculations.

#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors Rle
#' @importFrom BiocGenerics cbind
.rescale_batches <- function(batches, subset.row, averages) {
    reference <- do.call(pmin, averages)
    for (b in seq_along(batches)) {
        rescale <- reference / averages[[b]] 
        rescale[!is.finite(rescale)] <- 0

        if (!is.null(subset.row)) {
            batches[[b]] <- batches[[b]][subset.row,,drop=FALSE]
        }
        batches[[b]] <- batches[[b]] * rescale # preserves sparsity.
    }

    ncells.per.batch <- vapply(batches, ncol, FUN.VALUE=0L)
    batch.names <- .create_batch_names(names(batches), ncells.per.batch)
    SummarizedExperiment(list(corrected=do.call(cbind, batches)), 
        colData=DataFrame(batch=batch.names$ids)) 
}
