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
#' @importFrom S4Vectors metadata metadata<-
#' @importClassesFrom S4Vectors List
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

    # Different function calls for different input modes.
    common.args <- list(subset.row=subset.row, d=d, rotate.all=rotate.all, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
    if (length(mat.list)==1L) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }

        collected <- do.call(.multi_pca_single, c(list(mat=mat.list[[1]], batch=batch), common.args)) 
        rownames(collected[[1]]) <- colnames(originals[[1]])

        if (!preserve.single) {
            output <- divideIntoBatches(collected[[1]], batch, byrow=TRUE)
            output <- as(output$batches, "List")
            metadata(output) <- metadata(collected)
            collected <- output
        }

    } else {
        collected <- do.call(.multi_pca_list, c(list(mat.list=mat.list), common.args))
        for (i in seq_along(mat.list)) {
            rownames(collected[[i]]) <- colnames(originals[[i]])
        }
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

###########################################

#' @importFrom BiocSingular ExactParam bsdeferred
#' @importFrom BiocParallel SerialParam
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<- 
#' @importFrom Matrix crossprod
.multi_pca_list <- function(mat.list, subset.row=NULL, d=50, rotate.all=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    processed <- .process_listed_matrices_for_pca(mat.list, subset.row=subset.row, deferred=bsdeferred(BSPARAM))
    centered <- processed$centered
    scaled <- processed$scaled

    # Performing an SVD on the untransposed expression matrix,
    # and projecting the _unscaled_ matrices back into this space.
    svd.out <- .run_scaled_SVD(scaled, d=d, rotate.all=rotate.all, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    output <- centered
    for (idx in seq_along(centered)) {
        output[[idx]] <- as.matrix(crossprod(centered[[idx]], svd.out$u))
    }
    output <- as(output, "List")

    # Recording the rotation vectors. 
    extra.info <- list(
        rotation=.get_rotation_vectors(
            mat=mat.list[[1]], 
            subset.row=subset.row, 
            svd.out=svd.out, 
            rotate.all=rotate.all,
            process.FUN=.process_listed_matrices_for_pca, 
            process.args=list(mat.list=mat.list, deferred=bsdeferred(BSPARAM))
        )
    )

    # Computing the amount of variation explained.
    if (get.variance) {
        extra.info <- c(extra.info, .compute_var_explained(D=svd.out$d, nbatches=length(mat.list), scaled=scaled))
    }
        
    metadata(output) <- extra.info
    return(output)
}

#' @importFrom BiocGenerics ncol
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom BiocGenerics t
.process_listed_matrices_for_pca <- function(mat.list, subset.row, deferred=FALSE) {
    all.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }

        centers <- Matrix::rowMeans(current)
        all.centers <- all.centers + centers
        mat.list[[idx]] <- current
    }

    # grand average of centers (not using batch-specific centers, which makes compositional assumptions).
    all.centers <- all.centers/length(mat.list) 
        
    centered <- mat.list
    for (idx in seq_along(mat.list)) {
        current <- DelayedArray(mat.list[[idx]])
        current <- current - all.centers # centering each batch by the grand average.
        centered[[idx]] <- current
    } 

    if (deferred) {
        transposed <- lapply(mat.list, t)
        all.scale <- lapply(mat.list, function(x) rep(sqrt(ncol(x)), ncol(x)))

        # Deferring the centering and scaling operations.
        tmp <- DeferredMatrix(do.call(rbind, transposed), center=all.centers)
        scaled <- DeferredMatrix(t(tmp), scale=unlist(all.scale))

    } else {
        scaled <- centered 
        for (idx in seq_along(scaled)) {
            current <- centered[[idx]] 
            current <- current/sqrt(ncol(current)) # downweighting samples with many cells.
            scaled[[idx]] <- current
        }
        scaled <- do.call(cbind, scaled)
    }

    list(centered=centered, scaled=scaled)
}

#' @importFrom BiocSingular runSVD
.run_scaled_SVD <- function(scaled, d, rotate.all=FALSE, get.variance=FALSE, ...) {
    runSVD(scaled, 
        k=if (rotate.all || get.variance) d else 0,
        nu=d, 
        nv=if (rotate.all) d else 0, 
        ...
    )
}

.get_rotation_vectors <- function(mat, subset.row, svd.out, process.FUN, process.args, rotate.all=FALSE)
# inferring the for genes not in subset.row, if rotate.all=TRUE.
# This involves projecting the unused genes into the space defined by the PCs.
{
    if (rotate.all && !is.null(subset.row)) {
        subset.row <- .row_subset_to_index(mat, subset.row)

        leftovers <- do.call(process.FUN, c(process.args, list(subset.row=-subset.row)))$scaled
        leftover.u <- as.matrix(sweep(leftovers %*% svd.out$v, 2, svd.out$d, FUN="/"))
    
        rotation <- matrix(0, nrow(mat), ncol(svd.out$u))
        rotation[subset.row,] <- svd.out$u
        rotation[-subset.row,] <- leftover.u
    } else {
        rotation <- svd.out$u
    }
    return(rotation)
}

#' @importFrom DelayedArray DelayedArray
.compute_var_explained <- function(D, nbatches, scaled) {
    var.explained <- D^2 / nbatches
    var.total <- sum(DelayedArray(scaled)^2) / nbatches # DelayedArray needed as DeferredMatrix doesn't support arbitrary math.
    list(var.explained=var.explained, var.total=var.total)
}

###########################################

#' @importFrom S4Vectors List
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam bsdeferred
#' @importFrom Matrix crossprod
.multi_pca_single <- function(mat, batch, subset.row=NULL, d=50, rotate.all=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Similar to .multi_pca_list, but avoids the unnecessary
# overhead of splitting 'mat' into batch-specific matrices
# when you end up having to put them back together again anyway.
{
    processed <- .process_single_matrix_for_pca(mat, batch=batch, subset.row=subset.row, deferred=bsdeferred(BSPARAM))
    centered <- processed$centered
    scaled <- processed$scaled

    # Performing an SVD on the untransposed expression matrix,
    # and projecting the _unscaled_ matrices back into this space.
    svd.out <- .run_scaled_SVD(scaled, d=d, rotate.all=rotate.all, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    output <- List(as.matrix(crossprod(centered, svd.out$u)))

    # Recording the rotation vectors. 
    extra.info <- list(
        rotation=.get_rotation_vectors(
            mat=mat,
            subset.row=subset.row, 
            svd.out=svd.out, 
            rotate.all=rotate.all,
            process.FUN=.process_single_matrix_for_pca, 
            process.args=list(mat=mat, batch=batch, deferred=bsdeferred(BSPARAM))
        )
    )

    # Computing the amount of variation explained.
    if (get.variance) {
        extra.info <- c(extra.info, .compute_var_explained(D=svd.out$d, nbatches=length(unique(batch)), scaled=scaled))
    }
        
    metadata(output) <- extra.info
    return(output)
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom BiocGenerics t
.process_single_matrix_for_pca <- function(x, batch, subset.row, deferred=FALSE) {
    batch <- factor(batch)
    tab <- table(batch)
    w <- sqrt(as.numeric(tab[batch]))

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    # Computing the grand average of the centers.
    all.centers <- 0
    for (b in levels(batch)) {
        current <- x[,batch==b,drop=FALSE]
        all.centers <- all.centers + Matrix::rowMeans(current)
    }
    all.centers <- all.centers/length(tab)

    # Applying the centering and scaling.
    if (deferred) {
        centered <- t(DeferredMatrix(t(x), center=all.centers))
        scaled <- DeferredMatrix(centered, scale=w)
    } else {
        centered <- DelayedArray(x)
        centered <- centered - all.centers
        scaled <- sweep(centered, 2, w, "/", check.margin=FALSE)
    }

    list(centered=centered, scaled=scaled)
}
