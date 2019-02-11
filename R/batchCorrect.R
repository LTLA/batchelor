#' Batch correction methods
#'
#' A common interface for single-cell batch correction methods.
#'
#' @param ... Named data-dependent parameters to pass to the dispatched batch correction methods.
#' This should contain one or more matrix-like objects containing single-cell gene expression matrices.
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied.
#' @param batch A factor specifying the batch of origin for each cell if only one batch is supplied.
#' This will be ignored if two or more batches are supplied.
#' @param subset.row A vector specifying the subset of genes to use for correction.
#' Defaults to \code{NULL}, in which case all genes are used.
#' @param correct.all A logical scalar indicating whether to return corrected expression values for all genes, even if \code{subset.row} is set.
#' Used to ensure that the output is of the same dimensionality as the input.
#' @param assay.type A string or integer scalar specifying the assay to use for correction.
#' Only used for SingleCellExperiment inputs.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#' @param PARAM A \linkS4class{BatchelorParam} object specifying the batch correction method to dispatch to.
#' \linkS4class{ClassicMnnParam} will dispatch to \code{\link{mnnCorrect}};
#' \linkS4class{FastMnnParam} will dispatch to \code{\link{fastMNN}};
#' and \linkS4class{RescaleParam} will dispatch to \code{\link{rescaleBatches}}.
#' 
#' @return
#' A SingleCellExperiment where the first assay contains corrected gene expression values for all genes.
#' Corrected values should be returned for all genes if \code{subset.row=NULL} or if \code{correct.all=TRUE};
#' otherwise they should only be returned for the genes in the subset.
#'
#' Cells should be reported in the same order that they are supplied.
#' In cases with multiple batches, the cell identities are simply concatenated from successive objects in their specified order,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#'
#' The \code{colData} slot should contain \code{batch}, a vector specifying the batch of origin for each cell.
#'
#' @author Aaron Lun
#'
#' @details
#' Users can pass parameters to each method directly via \code{...} or via the constructors for \code{PARAM}.
#' While there is no restriction on which parameters go where, we recommend only passing data-agnostic and method-specific parameters to \code{PARAM}.
#' Data-dependent parameters - and indeed, the data themselves - should be passed in via \code{...}.
#' This means that different data sets can be used without modifying \code{PARAM}, and allows users to switch to a different algorithm by only changing \code{PARAM}.
#' 
#' Note that \code{get.spikes=FALSE} effectively modifies \code{subset.row} to exclude spike-in transcripts when SingleCellExperiment inputs are supplied.
#' This means that the reported SingleCellExperiment will not, by default, contain corrected expression values for spike-in transcripts unless \code{get.spikes=TRUE}.
#'
#' @seealso
#' \linkS4class{BatchelorParam} classes to determine dispatch.
#'
#' @examples
#' B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
#'
#' # Switching easily between batch correction methods.
#' m.out <- batchCorrect(B1, B2, PARAM=ClassicMnnParam())
#' f.out <- batchCorrect(B1, B2, PARAM=FastMnnParam(d=20))
#' r.out <- batchCorrect(B1, B2, PARAM=RescaleParam(pseudo.count=0))
#'
#' @rdname batchCorrect
#' @export
#' @import methods
setMethod("batchCorrect", "ClassicMnnParam", function(..., batch=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", get.spikes=FALSE, PARAM) {
    if (is.null(assay.type)) assay.type <- "logcounts"
    combined <- c(list(...), list(batch=batch, subset.row=subset.row, correct.all=correct.all, assay.type=assay.type, get.spikes=get.spikes), as.list(PARAM))
    do.call(mnnCorrect, combined)
})

#' @rdname batchCorrect
#' @export
setMethod("batchCorrect", "FastMnnParam", function(..., batch=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", get.spikes=FALSE, PARAM) {
    if (is.null(assay.type)) assay.type <- "logcounts"
    combined <- c(list(...), list(batch=batch, subset.row=subset.row, correct.all=correct.all, assay.type=assay.type, get.spikes=get.spikes), as.list(PARAM))
    combined$pc.input <- FALSE # Ensure we get an SCE out.
    combined$use.dimred <- NULL
    do.call(fastMNN, combined)
})

#' @rdname batchCorrect
#' @export
setMethod("batchCorrect", "RescaleParam", function(..., batch=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", get.spikes=FALSE, PARAM) {
    if (is.null(assay.type)) assay.type <- "logcounts"
    combined <- c(list(...), list(batch=batch, subset.row=subset.row, correct.all=correct.all, assay.type=assay.type, get.spikes=get.spikes), as.list(PARAM))
    do.call(rescaleBatches, combined)
})
