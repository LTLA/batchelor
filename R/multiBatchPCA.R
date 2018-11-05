#' Multi-batch PCA
#'
#' Perform a principal components analysis across multiple gene expression matrices to project all cells to a common low-dimensional space.
#'
#' @param ... Two or more matrices containing expression values (usually log-normalized).
#' Each matrix is assumed to represent one batch.
#' Alternatively, two or more SingleCellExperiment objects containing these matrices.
#' @param d An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.
#' @param subset.row See \code{?"\link{scran-gene-selection}"}.
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes See \code{?"\link{scran-gene-selection}"}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, see \code{\link{runSVD}} for details.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the SVD should be parallelized.
#'
#' @details
#' This function is roughly equivalent to \code{cbind}ing all matrices in \code{...} and performing PCA on the merged matrix.
#' The main difference is that each sample contributes equally to the identification of the rotation vectors.
#' Specifically, the mean vector used for centering is defined as the grand mean of the mean vectors within each batch.
#' Each batch's contribution to the gene-gene covariance matrix is also divided by the number of cells.
#' 
#' In effect, we weight the cells in each batch to mimic the situation where all batches have the same number of cells.
#' This ensures that the low-dimensional space can distinguish subpopulations in smaller batches.
#' Otherwise, batches with a large number of cells would dominate the PCA; the mean vector and covariance matrix would be almost fully defined by those batches.
#' 
#' If \code{...} contains SingleCellExperiment objects, any spike-in transcripts should be the same across all batches.
#' These will be removed prior to PCA unless \code{get.spikes=TRUE}, see \code{?"\link{scran-gene-selection}"} for details.
#'
#' @return
#' A \linkS4class{List} of numeric matrices where each matrix corresponds to a batch and contains the first \code{d} PCs (columns) for all cells in the batch (rows).
#' 
#' The metadata also contains a matrix of rotation vectors, which can be used to construct a low-rank approximation of the input matrices.
#'
#' Aaron Lun
#'
#' @seealso \code{\link{runSVD}}
#'
#' @examples 
#' d1 <- matrix(rnorm(5000), ncol=100)
#' d1[1:10,1:10] <- d1[1:10,1:10] + 2 # unique population in d1
#' d2 <- matrix(rnorm(2000), ncol=40)
#' d2[11:20,1:10] <- d2[11:20,1:10] + 2 # unique population in d2
#' 
#' out <- multiBatchPCA(d1, d2)
#' 
#' xlim <- range(c(out[[1]][,1], out[[2]][,1]))
#' ylim <- range(c(out[[1]][,2], out[[2]][,2]))
#' plot(out[[1]][,1], out[[1]][,2], col="red", xlim=xlim, ylim=ylim)
#' points(out[[2]][,1], out[[2]][,2], col="blue") 
#' 
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay
multiBatchPCA <- function(..., d=50, subset.row=NULL, assay.type="logcounts", get.spikes=FALSE, BSPARAM=NULL, BPPARAM=SerialParam()) 
# Performs a multi-sample PCA (i.e., batches).
# Each batch is weighted inversely by the number of cells when computing the gene-gene covariance matrix.
# This avoids domination by samples with a large number of cells.
#
# written by Aaron Lun
# created 4 July 2018
{
    mat.list <- list(...)

    if (.check_if_SCEs(mat.list)) {
		.check_spike_consistency(mat.list)
        subset.row <- .SCE_subset_genes(subset.row, mat.list[[1]], get.spikes)
        mat.list <- lapply(mat.list, assay, i=assay.type, withDimnames=FALSE)
    }

    .check_batch_consistency(mat.list)
    if (length(mat.list)==0L) {
        stop("at least one batch must be specified") 
    }

    .multi_pca(mat.list, subset.row=subset.row, d=d, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
}

#' @importFrom BiocSingular runSVD
#' @importFrom DelayedArray DelayedArray t
#' @importFrom BiocGenerics rowMeans cbind
#' @importFrom BiocParallel SerialParam
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<-
.multi_pca <- function(mat.list, subset.row=NULL, d=50, BSPARAM=NULL, BPPARAM=SerialParam()) 
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    all.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- DelayedArray(mat.list[[idx]])
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }

        centers <- rowMeans(current)
        all.centers <- all.centers + centers
        mat.list[[idx]] <- current
    }
    all.centers <- all.centers/length(mat.list) # grand average of centers (not using batch-specific centers, which makes compositional assumptions).

    centered <- scaled <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        current <- current - all.centers # centering each batch by the grand average.
        centered[[idx]] <- current
        current <- current/sqrt(ncol(current)) # downweighting samples with many cells.
        scaled[[idx]] <- current
    }

    # Performing an SVD on the untransposed expression matrix.
    final <- do.call(cbind, scaled)
    svd.out <- runSVD(final, k=0, nv=0, nu=d, BSPARAM=BSPARAM)

    # Projecting the _unscaled_ matrices back into this space.
    final <- centered
    for (idx in seq_along(centered)) {
        final[[idx]] <- as.matrix(t(centered[[idx]]) %*% svd.out$u) # replace with DA crossprod.
    }

    final <- as(final, "List")
    metadata(final) <- list(rotation=svd.out$u)
    return(final)
}
