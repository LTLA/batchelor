#' Multi-batch PCA
#'
#' Perform a principal components analysis across multiple gene expression matrices to project all cells to a common low-dimensional space.
#'
#' @param ... Two or more matrices containing expression values (usually log-normalized).
#' Each matrix is assumed to represent one batch.
#'
#' Alternatively, two or more \linkS4class{SingleCellExperiment} objects containing these matrices.
#'
#' Alternatively, one matrix or SingleCellExperiment can be supplied containing cells from all batches. 
#' This requires \code{batch} to also be specified.
#' @param batch A factor specifying the batch identity of each cell in the input data.
#' Ignored if \code{...} contains more than one argument.
#' @param d An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#' @param rotate.all A logical scalar indicating whether the reported rotation vectors should include genes that are excluded by a non-\code{NULL} value of \code{subset.row}.
#' @param get.variance A logical scalar indicating whether to return the (weighted) variance explained by each PC.
#' @param preserve.single A logical scalar indicating whether to combine the results into a single matrix if only one object was supplied in \code{...}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, see \code{\link{runSVD}} for details.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the SVD should be parallelized.
#'
#' @details
#' This function is roughly equivalent to \code{cbind}ing all matrices in \code{...} and performing PCA on the merged matrix.
#' The main difference is that each sample is forced to contribute equally to the identification of the rotation vectors.
#' Specifically, the mean vector used for centering is defined as the grand mean of the mean vectors within each batch.
#' Each batch's contribution to the gene-gene covariance matrix is also divided by the number of cells in that batch.
#' 
#' Our approach is to effectively weight the cells in each batch to mimic the situation where all batches have the same number of cells.
#' This ensures that the low-dimensional space can distinguish subpopulations in smaller batches.
#' Otherwise, batches with a large number of cells would dominate the PCA, i.e., the definition of the mean vector and covariance matrix.
#' This may reduce resolution of unique subpopulations in smaller batches that differ in a different dimension to the subspace of the larger batches.
#' 
#' If \code{...} contains SingleCellExperiment objects, any spike-in transcripts should be the same across all batches.
#' These will be removed prior to PCA unless \code{get.spikes=TRUE}.
#' If \code{subset.row} is specified and \code{get.spikes=FALSE}, only the non-spike-in specified features will be used. 
#'
#' Setting \code{rotate.all=TRUE} will report rotation vectors that span all genes, even when only a subset of genes are used for the PCA.
#' This is done by projecting all non-used genes into the low-dimensional \dQuote{cell space} defined by the first \code{d} components.
#'
#' If \code{BSPARAM} is defined with \code{deferred=TRUE}, the per-gene centering and per-cell scaling will be manually deferred during matrix multiplication.
#' This can greatly improve speeds when the input matrices are sparse, as deferred operations avoids loss of sparsity (at the cost of numerical precision).
#'
#' @return
#' A \linkS4class{List} of numeric matrices is returned where each matrix corresponds to a batch and contains the first \code{d} PCs (columns) for all cells in the batch (rows).
#'
#' If \code{preserve.single=TRUE} and \code{...} contains a single object, the List will only contain a single matrix.
#' This contains the first \code{d} PCs (columns) for all cells in the same order as supplied in the single input object.
#' 
#' The metadata contains \code{rotation}, a matrix of rotation vectors, which can be used to construct a low-rank approximation of the input matrices.
#' This has number of rows equal to the number of genes after any subsetting, except if \code{rotate.all=TRUE}, where the number of rows is equal to the genes before subsetting.
#'
#' If \code{get.variance=TRUE}, the metadata will also contain \code{var.explained}, the weighted variance explained by each PC;
#' and \code{var.total}, the total variance after weighting.
#'
#' @author
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
#' @importFrom S4Vectors metadata List metadata<-
#' @importFrom BiocGenerics colnames<- rownames<- colnames rownames
#' @importFrom BiocSingular ExactParam
multiBatchPCA <- function(..., batch=NULL, d=50, subset.row=NULL, rotate.all=FALSE, get.variance=FALSE, preserve.single=FALSE,
    assay.type="logcounts", get.spikes=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Performs a multi-sample PCA (i.e., batches).
# Each batch is weighted inversely by the number of cells when computing the gene-gene covariance matrix.
# This avoids domination by samples with a large number of cells.
#
# written by Aaron Lun
# created 4 July 2018
{
    originals <- mat.list <- list(...)
    if (length(mat.list)==0L) {
        stop("at least one batch must be specified") 
    }
    .check_batch_consistency(mat.list, byrow=TRUE)

    if (.check_if_SCEs(mat.list)) {
		.check_spike_consistency(mat.list)
        subset.row <- .SCE_subset_genes(subset.row, mat.list[[1]], get.spikes)
        mat.list <- lapply(mat.list, assay, i=assay.type, withDimnames=FALSE)
    }

    # Subsetting by 'batch'.
    do.split <- length(mat.list)==1L
    if (do.split) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        divided <- .divide_into_batches(mat.list[[1]], batch=batch, byrow=FALSE)
        mat.list <- divided$batches
    }

    collected <- .multi_pca(mat.list, subset.row=subset.row, d=d, rotate.all=rotate.all, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 

    # Adding row names.
    if (do.split) {
        b <- as.factor(batch)
        for (i in names(collected)) {
            rownames(collected[[i]]) <- colnames(originals[[1]])[b==i] 
        }
    } else {
        for (i in seq_along(mat.list)) {
            rownames(collected[[i]]) <- colnames(originals[[i]])
        }
    }

    # Reordering the output for correctness.
    if (do.split && preserve.single) {
        combined <- do.call(rbind, as.list(collected))
        d.reo <- divided$reorder
        combined <- combined[d.reo,,drop=FALSE]

        tmp <- List(combined)
        metadata(tmp) <- metadata(collected)
        collected <- tmp
    }

    # Adding dimension names to the rotation matrix.        
    rotation <- metadata(collected)$rotation
    if (rotate.all || is.null(subset.row)) {
        rnames <- rownames(originals[[1]])
    } else {
        if (is.character(subset.row)) { 
            rnames <- subset.row
        } else {
            rnames <- rownames(originals[[1]])[subset.row]
        }
    }
    rownames(rotation) <- rnames
    metadata(collected)$rotation <- rotation

    return(collected)
}

#' @importFrom BiocSingular runSVD ExactParam bsdeferred
#' @importFrom DelayedArray DelayedArray t
#' @importFrom BiocGenerics cbind nrow
#' @importFrom BiocParallel SerialParam
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<- normalizeSingleBracketSubscript
.multi_pca <- function(mat.list, subset.row=NULL, d=50, rotate.all=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    if (bsdeferred(BSPARAM)) {
        processed <- .process_deferred_matrices_for_pca(mat.list, subset.row)
    } else {
        processed <- .process_delayed_matrices_for_pca(mat.list, subset.row)
    }
    centered <- processed$centered
    scaled <- processed$scaled

    # Performing an SVD on the untransposed expression matrix,
    # and projecting the _unscaled_ matrices back into this space.
    svd.out <- runSVD(scaled, 
        k=if (rotate.all || get.variance) d else 0,
        nu=d, 
        nv=if (rotate.all) d else 0, 
        BSPARAM=BSPARAM)

    output <- centered
    for (idx in seq_along(centered)) {
        output[[idx]] <- as.matrix(t(centered[[idx]]) %*% svd.out$u) # replace with DA crossprod.
    }
    output <- as(output, "List")

    # Recording the rotation vectors. Optionally inferring them for genes not in subset.row.
    if (rotate.all && !is.null(subset.row)) {
        subset.row <- normalizeSingleBracketSubscript(subset.row, mat.list[[1]])

        if (bsdeferred(BSPARAM)) {
            leftovers <- .process_deferred_matrices_for_pca(mat.list, -subset.row)$scaled
        } else {
            leftovers <- .process_delayed_matrices_for_pca(mat.list, -subset.row)$scaled
        }
        leftover.u <- as.matrix(sweep(leftovers %*% svd.out$v, 2, svd.out$d, FUN="/"))

        rotation <- matrix(0, nrow(mat.list[[1]]), d)
        rotation[subset.row,] <- svd.out$u
        rotation[-subset.row,] <- leftover.u
    } else {
        rotation <- svd.out$u
    }
    extra.info <- list(rotation=rotation)

    # Computing the amount of variation explained.
    # Note the second DelayedArray, as DeferredMatrix doesn't support arbitrary math.
    if (get.variance) {
        nbatches <- length(mat.list)
        extra.info$var.explained <- svd.out$d^2 / nbatches
        extra.info$var.total <- sum(DelayedArray(scaled)^2) / nbatches
    }
        
    metadata(output) <- extra.info
    return(output)
}

#' @importFrom DelayedArray DelayedArray 
#' @importFrom BiocGenerics rowMeans ncol
.process_delayed_matrices_for_pca <- function(mat.list, subset.row) {
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

    # grand average of centers (not using batch-specific centers, which makes compositional assumptions).
    all.centers <- all.centers/length(mat.list) 

    centered <- scaled <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        current <- current - all.centers # centering each batch by the grand average.
        centered[[idx]] <- current
        current <- current/sqrt(ncol(current)) # downweighting samples with many cells.
        scaled[[idx]] <- current
    }

    list(centered=centered, scaled=do.call(cbind, scaled))
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom BiocGenerics t
.process_deferred_matrices_for_pca <- function(mat.list, subset.row) {
    all.centers <- 0
    all.scale <- vector("list", length(mat.list))

    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }

        centers <- Matrix::rowMeans(current)
        all.centers <- all.centers + centers
        mat.list[[idx]] <- current

        all.scale[[idx]] <- rep(sqrt(ncol(current)), ncol(current))
    }

    # grand average of centers (not using batch-specific centers, which makes compositional assumptions).
    all.centers <- all.centers/length(mat.list) 

    centered <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- DelayedArray(mat.list[[idx]])
        current <- current - all.centers # centering each batch by the grand average.
        centered[[idx]] <- current
    }

    # Deferring the centering and scaling operations.
    m1 <- DeferredMatrix(Matrix::t(do.call(cbind, mat.list)), center=all.centers)
    m2 <- DeferredMatrix(t(m1), scale=unlist(all.scale))
   
    list(centered=centered, scaled=m2)
}
