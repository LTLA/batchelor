#' Multi-batch PCA
#'
#' Perform a principal components analysis across multiple gene expression matrices to project all cells to a common low-dimensional space.
#'
#' @param ... Two or more matrices containing expression values (usually log-normalized).
#' Each matrix is assumed to represent one batch and should contain the same number of rows, corresponding to the same genes (in the same order).
#'
#' Alternatively, two or more \linkS4class{SingleCellExperiment} objects containing these matrices.
#' Note the same restrictions described above for gene expression matrix inputs.
#'
#' Alternatively, one matrix or SingleCellExperiment can be supplied containing cells from all batches. 
#' This requires \code{batch} to also be specified.
#' @param batch A factor specifying the batch identity of each cell in the input data.
#' Ignored if \code{...} contains more than one argument.
#' @param weights Numeric vector of length equal to the number of entries in \code{...}, specifying the scaling of the weight of each batch.
#' This defaults to 1 for all batches.
#' @param d An integer scalar specifying the number of dimensions to keep from the initial multi-sample PCA.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes Deprecated, a logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#' @param get.all.genes A logical scalar indicating whether the reported rotation vectors should include genes 
#' that are excluded by a non-\code{NULL} value of \code{subset.row}.
#' @param rotate.all A deprecated synonym for \code{get.all.genes}.
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
#' When \code{weights} is set, this will scale the weight of each batch by the specified value. 
#' For example, each batch may represent one replicate, with multiple replicates per study.
#' In such cases, it may be more appropriate to ensure that each \emph{study} has equal weight.
#' This is done by assigning a value of \code{weights} to each replicate that is inversely proportional to the number of replicates in the same study - see Examples.
#' 
#' If \code{...} contains SingleCellExperiment objects, any spike-in transcripts should be the same across all batches.
#' These will be removed prior to PCA unless \code{get.spikes=TRUE}.
#' If \code{subset.row} is specified and \code{get.spikes=FALSE}, only the non-spike-in specified features will be used. 
#'
#' Setting \code{get.all.genes=TRUE} will report rotation vectors that span all genes, even when only a subset of genes are used for the PCA.
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
#' # Examining results.
#' xlim <- range(c(out[[1]][,1], out[[2]][,1]))
#' ylim <- range(c(out[[1]][,2], out[[2]][,2]))
#' plot(out[[1]][,1], out[[1]][,2], col="red", xlim=xlim, ylim=ylim)
#' points(out[[2]][,1], out[[2]][,2], col="blue") 
#'
#' # Using the weighting scheme, assuming that 'd2' and 'd3'
#' # are replicates and should contribute the same combined
#' # weight as 'd1'.
#' d3 <- d2 + 5
#' out <- multiBatchPCA(d1, d2, d3, weights=c(1, 0.5, 0.5))
#' 
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata<-
#' @importClassesFrom S4Vectors List
#' @importFrom BiocGenerics colnames<- rownames<- colnames rownames
#' @importFrom BiocSingular ExactParam
multiBatchPCA <- function(..., batch=NULL, d=50, subset.row=NULL, weights=NULL,
    get.all.genes=FALSE, rotate.all=FALSE, get.variance=FALSE, preserve.single=FALSE, 
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
    checkBatchConsistency(mat.list)

    is.sce <- checkIfSCE(mat.list)
    if (any(is.sce)) {
        sce.in <- mat.list[is.sce]
        checkSpikeConsistency(sce.in)
        subset.row <- .SCE_subset_genes(subset.row, sce.in[[1]], get.spikes)
        mat.list[is.sce] <- lapply(sce.in, assay, i=assay.type)
    }

    if (rotate.all) {
        .Deprecated(msg="'rotate.all=TRUE' is deprecated.\nUse 'get.all.genes=TRUE' instead.")
        get.all.genes <- TRUE
    }

    # Different function calls for different input modes.
    common.args <- list(subset.row=subset.row, d=d, get.all.genes=get.all.genes, 
        weights=weights, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM) 
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

    # Adding dimension names to the metadata fields.
    if (get.all.genes || is.null(subset.row)) {
        rnames <- rownames(originals[[1]])
    } else {
        if (is.character(subset.row)) { 
            rnames <- subset.row
        } else {
            rnames <- rownames(originals[[1]])[subset.row]
        }
    }
    rownames(metadata(collected)$rotation) <- rnames
    names(metadata(collected)$centers) <- rnames

    collected
}

###########################################

#' @importFrom BiocSingular ExactParam bsdeferred
#' @importFrom BiocParallel SerialParam bpparam register bpisup bpstart bpstop
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<- 
#' @importFrom Matrix crossprod
.multi_pca_list <- function(mat.list, subset.row=NULL, d=50, weights=NULL,
    get.all.genes=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    FUN <- function(keep) {
        .process_listed_matrices_for_pca(mat.list, weights=weights, subset.row=keep, deferred=bsdeferred(BSPARAM))
    }
    processed <- FUN(subset.row)
    centered <- processed$centered
    scaled <- processed$scaled
    centers <- processed$centers

    # Setting up the parallelization environment (for crossprod,
    # but runSVD might as well take advantage of it).
    old <- bpparam()
    register(BPPARAM)
    on.exit(register(old))

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    # Performing an SVD on the untransposed _scaled_ expression matrix,
    # then projecting the _unscaled_ matrices back into this space.
    svd.out <- .run_scaled_SVD(scaled, d=d, get.all.genes=get.all.genes, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM)

    output <- centered
    for (idx in seq_along(centered)) {
        output[[idx]] <- as.matrix(crossprod(centered[[idx]], svd.out$u))
    }
    output <- as(output, "List")

    # Recording useful PCA metadata.
    metadata(output) <- .get_pca_metadata(mat=mat.list[[1]], subset.row=subset.row, 
        get.all.genes=get.all.genes, FUN=FUN, svd.out=svd.out, centers=centers, 
        get.variance=get.variance, nbatches=length(mat.list), scaled=scaled)

    output
}

#' @importFrom BiocGenerics ncol rbind
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom Matrix t rowMeans
.process_listed_matrices_for_pca <- function(mat.list, weights, subset.row, deferred=FALSE) {
    if (is.null(weights)) {
        weights <- rep(1, length(mat.list))
    } else {
        if (!is.null(names(weights)) && !identical(names(mat.list), names(weights))) {
            stop("'names(weights)' should be the same as names of '...'")
        } else if (length(mat.list)!=length(weights)) {
            stop("'length(weights)' should be the same as number of entries in '...'")
        }
    }

    # Computing the grand average of centers, and using that to center each batch.
    # (Not using batch-specific centers, which makes compositional assumptions.)
    grand.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }

        centers <- rowMeans(current)
        grand.centers <- grand.centers + centers * weights[idx]
        mat.list[[idx]] <- current
    }
    grand.centers <- grand.centers/sum(weights)

    # Creating a dummy list but with the same names.
    # Ensures that we catch errors from not centering.
    centered <- lapply(mat.list, function(x) NULL)

    if (deferred) {
        transposed <- all.scale <- vector("list", length(mat.list))
        for (idx in seq_along(mat.list)) {
            current <- t(mat.list[[idx]]) 
            transposed[[idx]] <- current
            current <- DeferredMatrix(current, center=grand.centers)
            centered[[idx]] <- t(current)
            w <- sqrt(nrow(current) / weights[idx])
            all.scale[[idx]] <- rep(w, nrow(current))
        } 

        # (Deferred) scaling to implicitly downweight samples with many cells.
        tmp <- DeferredMatrix(do.call(rbind, transposed), center=grand.centers)
        scaled <- DeferredMatrix(t(tmp), scale=unlist(all.scale))

    } else {
        for (idx in seq_along(mat.list)) {
            current <- DelayedArray(mat.list[[idx]])
            current <- current - grand.centers
            centered[[idx]] <- current
        }

        # (Delayed) scaling to implicitly downweight samples with many cells.
        scaled <- vector("list", length(mat.list))
        for (idx in seq_along(centered)) {
            current <- centered[[idx]] 
            w <- sqrt(ncol(current) / weights[idx])
            current <- current/w
            scaled[[idx]] <- current
        }
        scaled <- do.call(cbind, scaled)
    }

    list(centered=centered, scaled=scaled, centers=grand.centers)
}

###########################################

#' @importFrom BiocSingular runSVD
.run_scaled_SVD <- function(scaled, d, get.all.genes=FALSE, get.variance=FALSE, ...) {
    runSVD(scaled, 
        k=if (get.all.genes || get.variance) d else 0,
        nu=d, 
        nv=if (get.all.genes) d else 0, 
        ...
    )
}

#' @importFrom DelayedArray DelayedArray sweep
.get_pca_metadata <- function(mat, subset.row, 
    get.all.genes, FUN, svd.out, centers,
    get.variance, nbatches, scaled)
# Inferring rotation vectors for genes not in subset.row, if get.all.genes=TRUE.
# This involves projecting the unused genes into the space defined by the PCs.
# We also fill in the missing center values.
{
    if (get.all.genes && !is.null(subset.row)) {
        subset.row <- .row_subset_to_index(mat, subset.row)

        leftovers <- FUN(-subset.row)
        left.scaled <- leftovers$scaled
        leftover.u <- as.matrix(sweep(left.scaled %*% svd.out$v, 2, svd.out$d, FUN="/"))
    
        rotation <- matrix(0, nrow(mat), ncol(svd.out$u))
        rotation[subset.row,] <- svd.out$u
        rotation[-subset.row,] <- leftover.u

        all.centers <- numeric(nrow(mat))
        all.centers[subset.row] <- centers
        all.centers[-subset.row] <- leftovers$centers
    } else {
        rotation <- svd.out$u
        all.centers <- centers
    }

    output <- list(rotation=rotation, centers=all.centers)

    # Adding the variance explained.
    if (get.variance) {
        output$var.explained <- svd.out$d^2 / nbatches
        output$var.total <- sum(scaled^2) / nbatches 
    }

    output
}

###########################################

#' @importFrom S4Vectors List
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam bsdeferred
#' @importFrom Matrix crossprod
.multi_pca_single <- function(mat, batch, subset.row=NULL, weights=NULL, d=50, 
    get.all.genes=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
# Similar to .multi_pca_list, but avoids the unnecessary
# overhead of splitting 'mat' into batch-specific matrices
# when you end up having to put them back together again anyway.
{
    FUN <- function(keep) {
        .process_single_matrix_for_pca(mat, batch=batch, weights=weights, subset.row=keep, deferred=bsdeferred(BSPARAM))
    }
    processed <- FUN(subset.row)
    centered <- processed$centered
    scaled <- processed$scaled
    centers <- processed$centers

    # Performing an SVD on the untransposed expression matrix,
    # and projecting the _unscaled_ matrices back into this space.
    svd.out <- .run_scaled_SVD(scaled, d=d, get.all.genes=get.all.genes, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    output <- List(as.matrix(crossprod(centered, svd.out$u)))

    # Recording the rotation vectors. 
    metadata(output) <- .get_pca_metadata(mat=mat, subset.row=subset.row, 
        get.all.genes=get.all.genes, FUN=FUN, svd.out=svd.out, centers=centers, 
        get.variance=get.variance, nbatches=length(unique(batch)), scaled=scaled)

    output
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular DeferredMatrix
#' @importFrom Matrix t rowMeans
.process_single_matrix_for_pca <- function(x, batch, weights, subset.row, deferred=FALSE) {
    batch <- factor(batch)
    tab <- table(batch)
    nms <- names(tab)
    tab <- as.numeric(tab) 
    names(tab) <- nms

    if (!is.null(weights)) {
        weights <- weights[levels(batch)]
        if (any(is.na(weights))) {
            stop("'weights' should be named with all levels of 'batch'")
        }
    } else {
        weights <- tab
        weights[] <- 1
    }

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    # Computing the grand average of the centers.
    centers <- 0
    for (b in levels(batch)) {
        current <- x[,batch==b,drop=FALSE]
        centers <- centers + rowMeans(current) * weights[b]
    }
    centers <- centers/sum(weights)

    # Applying the centering and scaling.
    w <- sqrt(tab[batch] / weights[batch])
    if (deferred) {
        centered <- t(DeferredMatrix(t(x), center=centers))
        scaled <- DeferredMatrix(centered, scale=w)
    } else {
        centered <- DelayedArray(x)
        centered <- centered - centers
        scaled <- sweep(centered, 2, w, "/", check.margin=FALSE)
    }

    list(centered=centered, scaled=scaled, centers=centers)
}
