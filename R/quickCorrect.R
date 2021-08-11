#' Quickly perform batch correction
#'
#' Quickly perform batch correction by intersecting the gene features, normalizing and log-transforming,
#' modelling per-gene variances and identifying highly variable genes, and then applying the correction algorithm of choice.
#'
#' @inheritParams batchCorrect
#' @param multi.norm.args Named list of further arguments to pass to \code{\link{multiBatchNorm}}.
#' @param model.var.args Named list of further arguments to pass to \code{\link[scran]{modelGeneVar}}.
#' @param precomputed List of \link{DataFrame}s containing precomputed variance modelling results.
#' This should be a list of the same length as the number of entries in \code{...},
#' and each should have the same row names as the corresponding entry of \code{...}.
#' @param hvg.args Named list of further arguments to pass to \code{\link[scran]{getTopHVGs}}.
#' By default, we take the top 5000 genes with the highest variances.
#' 
#' @details
#' This function wraps the sequence of typical steps required to obtain corrected expression values.'
#' \enumerate{
#' \item Intersecting each batch to the universe of common features with \code{\link{intersectRows}}.
#' \item Applying normalization and log-transformation to the batches with \code{\link{multiBatchNorm}}.
#' \item Modelling the per-gene variance with \code{\link[scran]{modelGeneVar}}.
#' If \code{precomputed} is supplied, the precomputed results are used instead.
#' \item Identifying highly variable genes with \code{\link[scran]{getTopHVGs}}.
#' These genes will be used in the correction, though corrected values for all genes can be returned by setting \code{correct.all=TRUE}.
#' \item Applying the batch correction algorithm of choice with \code{\link{batchCorrect}}, as specified by \code{PARAM}.
#' }
#' 
#' The default of \code{correct.all=TRUE} differs from that of other functions.
#' This is because the subsetting to HVGs is done internally here, and we avoid surprises by returning results for all genes in the input object(s).
#' In contrast, the other functions require explicit subsetting via \code{subset.row=} and it is expected that users will set \code{correct.all=} if all genes are desired.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{dec}, a \linkS4class{DataFrame} containing the combined variance modelling results across all batches.
#' \item \code{hvgs}, a character vector of the genes selected to use in the correction.
#' \item \code{corrected}, a \linkS4class{SingleCellExperiment} containing the corrected values across all cells in all batches.
#' }
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{intersectRows}}, for the intersection to obtain a universe of genes.
#' 
#' \code{\link{multiBatchNorm}}, to perform the normalization.
#'
#' \code{\link[scran]{modelGeneVar}} and \code{\link[scran]{getTopHVGs}}, to identify the top HVGs for use in correction.
#'
#' \code{\link{batchCorrect}}, to dispatch to the desired correction algorithm.
#' 
#' @examples
#' d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
#' sce1 <- SingleCellExperiment(list(counts=d1))
#' sizeFactors(sce1) <- runif(ncol(d1))
#' rownames(sce1) <- paste0("GENE", 1:500)
#' 
#' d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
#' sce2 <- SingleCellExperiment(list(counts=d2))
#' sizeFactors(sce2) <- runif(ncol(d2))
#' rownames(sce2) <- paste0("GENE", 201:700)
#' 
#' # Fire and forget:
#' set.seed(1000)
#' output <- quickCorrect(sce1, sce2) 
#' output$corrected
#'
#' @export
quickCorrect <- function(..., 
    batch=NULL, 
    restrict=NULL, 
    correct.all=TRUE, 
    assay.type="counts", 
    PARAM=FastMnnParam(), 
    multi.norm.args=list(),
    precomputed=NULL,
    model.var.args=list(),
    hvg.args=list(n=5000))
{
    all.batches <- intersectRows(...)
    single.batch <- length(all.batches)==1L

    # Perform multibatchNorm.
    all.batches <- do.call(multiBatchNorm, c(
        all.batches, 
        list(batch=batch, assay.type=assay.type, preserve.single=TRUE),
        multi.norm.args    
    ))

    # Perform HVG calling within each batch.
    if (is.null(precomputed)) {
        if (single.batch) { 
            dec <- do.call(scran::modelGeneVar, c(list(all.batches, block=batch), model.var.args))
        } else {
            sub.dec <- vector("list", length(all.batches))
            for (i in seq_along(all.batches)) {
                sub.dec[[i]] <- do.call(scran::modelGeneVar, c(list(all.batches[[i]]), model.var.args))
            }
            dec <- do.call(scran::combineVar, sub.dec)
        }
    } else {
        if (single.batch) {
            if (length(precomputed)!=1L){ 
                stop("'length(precomputed)' is not the same as the number of objects in '...'")
            }
            dec <- precomputed[[1]]
        } else {
            if (length(precomputed)!=length(all.batches)) {
                stop("'length(precomputed)' is not the same as the number of objects in '...'")
            }
            universe <- rownames(all.batches[[1]])
            sub.dec <- lapply(precomputed, function(x) x[universe,,drop=FALSE])
            dec <- do.call(scran::combineVar, sub.dec)
        }
    }

    hvgs <- do.call(scran::getTopHVGs, c(list(dec), hvg.args))

    # Finally, get around to the batch correction.
    corrected <- batchCorrect(all.batches, batch=batch, restrict=restrict, correct.all=correct.all, subset.row=hvgs, PARAM=PARAM)

    list(dec=dec, hvgs=hvgs, corrected=corrected)
}
