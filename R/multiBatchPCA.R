#' Multi-batch PCA
#'
#' Perform a principal components analysis across multiple gene expression matrices to project all cells to a common low-dimensional space.
#'
#' @param ... Two or more matrices containing expression values (usually log-normalized).
#' Each matrix is assumed to represent one batch and should contain the same number of rows, 
#' corresponding to the same genes in the same order.
#'
#' Alternatively, two or more \linkS4class{SingleCellExperiment} objects containing these matrices.
#' Note the same restrictions described above for gene expression matrix inputs.
#'
#' Alternatively, one matrix or SingleCellExperiment can be supplied containing cells from all batches. 
#' This requires \code{batch} to also be specified.
#'
#' Alternatively, one or more lists of matrices or SingleCellExperiments can be provided;
#' this is flattened as if the objects inside were passed directly to \code{...}.
#' @param batch A factor specifying the batch identity of each cell in the input data.
#' Ignored if \code{...} contains more than one argument.
#' @param weights Numeric vector, logical scalar or list specifying the weighting scheme to use, see below for details.
#' @param d An integer scalar specifying the number of dimensions to keep from the PCA.
#' Alternatively \code{NA}, in which case the PCA step is omitted entirely - see details below.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.all.genes A logical scalar indicating whether the reported rotation vectors should include genes 
#' that are excluded by a non-\code{NULL} value of \code{subset.row}.
#' @param get.variance A logical scalar indicating whether to return the (weighted) variance explained by each PC.
#' @param preserve.single A logical scalar indicating whether to combine the results into a single matrix if only one object was supplied in \code{...}.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA, see \code{\link{runSVD}} for details.
#' @param deferred A logical scalar used to overwrite the \code{deferred} status of \code{BSPARAM} for greater speed.
#' Set to \code{NULL} to use the supplied status in \code{BSPARAM} directly.
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
#' It is usually recommended to set \code{subset.row} to a subset of interesting features (e.g., highly variable genes).
#' This reduces computational time and eliminates uninteresting noise that could interfere with identification of the most relevant axes of variation.
#' Setting \code{get.all.genes=TRUE} will report rotation vectors that span all genes, even when only a subset of genes are used for the PCA.
#' This is done by projecting all non-used genes into the low-dimensional \dQuote{cell space} defined by the first \code{d} components.
#'
#' @section Tuning the weighting:
#' By default, \code{weights=NULL} or \code{TRUE} will use the default weights,
#' i.e., the reciprocal of the number of cells in each batch.
#' This equalizes the contribution of batches with different numbers of cells as described above.
#'
#' If \code{weights=FALSE}, no weighting is performed.
#' This means that larger batches will drive the PCA,
#' which may be desirable when dealing with technical replicates where there is no concern about unique subpopulations in smaller batches.
#'
#' If \code{weights} is a numeric vector, it is expected to be of the same length 
#' (and, if named, have the same names) as the entries in \code{...}.
#' Each entry of the vector is used to scale the default weight of the corresponding batch.
#' This allows users to fine-tune the contribution of each batch in situations with multiple levels of heterogeneity.
#' For example, consider merging data from multiple donors where each donor contains a variable number of batches.
#' In such cases, it may be more appropriate to ensure that each donor has equal weight, rather than each batch.
#' This is done by assigning a value of \code{weights} to each replicate that is inversely proportional to the number of batches for the same donor - see Examples.
#'
#' Alternatively, \code{weights} can be a list representing a tree-like structure, 
#' identical to the tree controlling the merge order in \code{\link{fastMNN}}.
#' Here, weights are scaled so that each partition in the tree yields subtrees with identical scaled weights (summed across each subtree's children).
#' This allows us to easily adjust the weighting scheme for hierarchical batch structures like the replicate-study scenario described above.
#'
#' @section Performing the PCA:
#' Most of the parameters related to the PCA are controlled by \code{BSPARAM}, which determines the choice and parameterization of SVD algorithms.
#' The default is to use \code{\link{IrlbaParam}} though \code{\link{RandomParam}} is often faster (at the cost of some accuracy).
#' Most choices of \code{BSPARAM} will interactly sanely with \code{BPPARAM} to achieve parallelization during the PCA itself.
#'
#' With the default \code{deferred=TRUE}, the per-gene centering and per-cell scaling will be deferred during matrix multiplication in approximate SVD algorithms.
#' This is much faster when the input matrices are sparse, as deferred operations avoids loss of sparsity at the cost of numerical precision.
#' If \code{deferred=NULL}, the use of deferred scaling is determined by the setting within \code{BSPARAM} itself - see \code{?\link{bsdeferred}} for details.
#'
#' If \code{d=NA}, no PCA is performed.
#' Instead, the centered matrices are transposed and returned directly, while the rotation matrix is set to an identity matrix.
#' This allows developers to easily switch between PCA-based approximations versus the underlying dataset in their functions by simply changing \code{d}.
#' (In some settings, one can even interpret \code{d=NA} as the maximum \code{d}, as Euclidean distance calculations between cells will be identical to a full-rank PC matrix.)
#' Note that, in this mode, no projection will be done with \code{get.all.genes=TRUE}; rather, genes not in \code{subset.row} will simply have rotation values of zero.
#' If \code{get.variance=TRUE}, the weighted variances of the individual genes are returned.
#'
#' @return
#' A \linkS4class{List} of numeric matrices is returned where each matrix corresponds to a batch and contains the first \code{d} PCs (columns) for all cells in the batch (rows).
#'
#' If \code{preserve.single=TRUE} and \code{...} contains a single object, the List will only contain a single matrix.
#' This contains the first \code{d} PCs (columns) for all cells in the same order as supplied in the single input object.
#' 
#' The metadata contains \code{rotation}, a matrix of rotation vectors, which can be used to construct a low-rank approximation of the input matrices.
#' This has number of rows equal to the number of genes after any subsetting, except if \code{get.all.genes=TRUE}, where the number of rows is equal to the genes before subsetting.
#'
#' If \code{get.variance=TRUE}, the metadata will also contain \code{var.explained}, a numeric vector containing the weighted variance explained by each of \code{d} PCs;
#' and \code{var.total}, the total variance across all components after weighting.
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
#' # PCA defaults to IRLBA, so we need to set the seed.
#' set.seed(10)
#' out <- multiBatchPCA(d1, d2, d=10)
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
#' set.seed(10)
#' out <- multiBatchPCA(d1, d2, d3, d=10, weights=c(1, 0.5, 0.5))
#'
#' set.seed(10)
#' alt <- multiBatchPCA(d1, d2, d3, d=10, weights=list(1, list(2, 3))) 
#' stopifnot(all.equal(out, alt)) # As they are the same.
#' 
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata<-
#' @importClassesFrom S4Vectors List
#' @importFrom BiocGenerics colnames<- rownames<- colnames rownames
#' @importFrom BiocSingular ExactParam
#' @importFrom scuttle .bpNotSharedOrUp .unpackLists
multiBatchPCA <- function(..., batch=NULL, d=50, subset.row=NULL, weights=NULL,
    get.all.genes=FALSE, get.variance=FALSE, preserve.single=FALSE, assay.type="logcounts", 
    BSPARAM=IrlbaParam(), deferred=TRUE, BPPARAM=SerialParam()) 
{
    mat.list <- originals <- .unpackLists(...)
    if (length(mat.list)==0L) {
        stop("at least one batch must be specified") 
    }
    checkBatchConsistency(mat.list)

    is.sce <- checkIfSCE(mat.list)
    if (any(is.sce)) {
        mat.list[is.sce] <- lapply(mat.list[is.sce], assay, i=assay.type)
    }

    # Setting up the parallelization environment.
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    BSPARAM <- .set_deferred(BSPARAM, deferred)

    # Different function calls for different input modes.
    common.args <- list(subset.row=subset.row, d=d, get.all.genes=get.all.genes, 
        weights=weights, get.variance=get.variance, BSPARAM=BSPARAM, BPPARAM=BPPARAM)

    if (length(mat.list)==1L) {
        .check_valid_batch(mat.list[[1]], batch)

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
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom methods as
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors metadata<- 
#' @importFrom Matrix crossprod
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.multi_pca_list <- function(mat.list, subset.row=NULL, d=50, weights=NULL,
    get.all.genes=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam())
# Internal function that uses DelayedArray to do the centering and scaling,
# to avoid actually realizing the matrices in memory.
{
    # NOTE: need to set this up here, as fastMNN() calls this function, not multiBatchPCA().
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    FUN <- function(keep) {
        .process_listed_matrices_for_pca(mat.list, weights=weights, subset.row=keep, deferred=bsdeferred(BSPARAM))
    }
    processed <- FUN(subset.row)
    centered <- processed$centered
    scaled <- processed$scaled
    centers <- processed$centers

    if (!is.na(d)) { 
        # Performing an SVD on the untransposed _scaled_ expression matrix,
        # then projecting the _unscaled_ matrices back into this space.
        svd.out <- .run_scaled_SVD(scaled, d=d, 
            get.all.genes=get.all.genes, get.variance=get.variance, 
            BSPARAM=BSPARAM, BPPARAM=BPPARAM)

        output <- centered
        for (idx in seq_along(centered)) {
            output[[idx]] <- as.matrix(crossprod(centered[[idx]], svd.out$u))
        }
        output <- as(output, "List")

        metadata(output) <- .make_pca_metadata(mat=mat.list[[1]], subset.row=subset.row, 
            get.variance=get.variance, get.all.genes=get.all.genes, FUN=FUN, 
            svd.out=svd.out, centers=centers, nbatches=length(mat.list), scaled=scaled)
    } else {
        # Returning the scaled matrix directly if d=NA.
        output <- centered
        for (idx in seq_along(centered)) {
            output[[idx]] <- as.matrix(t(centered[[idx]]))
        }
        output <- as(output, "List")

        metadata(output) <- .make_fake_metadata(mat=mat.list[[1]], subset.row=subset.row,
            get.variance=get.variance, get.all.genes=get.all.genes, scaled=scaled)
    }

    output
}

#' @importFrom BiocGenerics ncol rbind
#' @importFrom DelayedArray DelayedArray
#' @importFrom ScaledMatrix ScaledMatrix
#' @importFrom Matrix t rowMeans
#' @importFrom beachmat realizeFileBackedMatrix
.process_listed_matrices_for_pca <- function(mat.list, weights, subset.row, deferred=FALSE) {
    weights <- .construct_weight_vector(vapply(mat.list, ncol, 0L), weights)

    # Computing the grand average of centers, and using that to center each batch.
    # (Not using batch-specific centers, which makes compositional assumptions.)
    grand.centers <- 0
    for (idx in seq_along(mat.list)) {
        current <- mat.list[[idx]]
        if (!is.null(subset.row)) {
            current <- current[subset.row,,drop=FALSE]
        }
        current <- realizeFileBackedMatrix(current)

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
            current <- ScaledMatrix(current, center=grand.centers)
            centered[[idx]] <- t(current)
            w <- sqrt(nrow(current) / weights[idx])
            all.scale[[idx]] <- rep(w, nrow(current))
        } 

        # (Deferred) scaling to implicitly downweight samples with many cells.
        tmp <- ScaledMatrix(do.call(rbind, transposed), center=grand.centers)
        scaled <- ScaledMatrix(t(tmp), scale=unlist(all.scale))

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

# Returns a vector of weights for each batch. If 'ncells' is named, the output
# is also expected to be named, as this is used in the single matrix
# implementation below.
.construct_weight_vector <- function(ncells, weights) {
    if (is.null(weights) || isTRUE(weights)) {
        ncells[] <- 1
        ncells
    } else if (is.list(weights)) {
        # Assume that we have a tree, and each split is of equal weight.
        w <- unlist(.get_list_weights(weights))
        ids <- unlist(weights)

        if (is.numeric(ids)) {
            if (!isTRUE(all.equal(sort(ids), seq_along(ncells)))) {
                stop("invalid integer indices in tree-like 'weights'")
            }
        } else if (is.character(ids)) {
            if (!identical(unname(sort(ids)), names(ncells))) {
                stop("names in tree-like 'weights' do not match names in '...'")
            }
        } else {
            stop("invalid children type in tree-like 'weights'")
        }

        # Piggy-backing off a named 'ncells' vector.
        ncells[ids] <- w
        ncells
    } else if (isFALSE(weights)) {
        ncells
    } else {
        if (!is.null(names(weights)) && !identical(names(ncells), names(weights))) {
            stop("'names(weights)' should be the same as names of '...'")
        } else if (length(ncells)!=length(weights)) {
            stop("'length(weights)' should be the same as number of entries in '...'")
        }
        weights
    }
}

# Recurse through the tree and assign equal weight at each split.
# This allows us to handle the weighting for hierarchical batches.
.get_list_weights <- function(tree, current=1) {
    reweight <- current/length(tree)
    if (!is.list(tree)) {
        rep(reweight, length(tree))
    } else {
        for (i in seq_along(tree)) {
            if (is.list(tree[[i]]) || length(tree[[i]]) > 1L) {
                tree[[i]] <- .get_list_weights(tree[[i]], current=reweight)
            } else {
                tree[[i]] <- reweight
            }
        }
        tree
    }
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
.make_pca_metadata <- function(mat, subset.row, get.variance, get.all.genes, FUN, svd.out, centers, nbatches, scaled)
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

    if (get.variance) {
        # Adding the variance explained. We actually don't need to divide by the
        # number of cells here, because the weighting already divided by the number
        # of cells, such that each batch contributes a weight of 1... so we just
        # divide by the number of batches.
        output$var.explained <- svd.out$d^2 / nbatches

        # 'scaled' should already be centered, so we dispense with that and just
        # compute the sum of squares.
        output$var.total <- sum(scaled^2) / nbatches
    }

    output
}

#' @importFrom Matrix sparseMatrix Diagonal
#' @importFrom DelayedMatrixStats rowVars
.make_fake_metadata <- function(mat, subset.row, get.variance, get.all.genes, scaled) {
    if (get.all.genes && !is.null(subset.row)) {
        subset.row <- .row_subset_to_index(mat, subset.row)
        output <- list(
             rotation=sparseMatrix(x=rep(1, length(subset.row)),
                 i=subset.row, j=seq_along(subset.row), 
                 dims=c(nrow(mat), length(subset.row))),
             centers=numeric(nrow(mat))
         )
    } else {
        output <- list(
             rotation=Diagonal(nrow(scaled)),
             centers=numeric(nrow(scaled))
         )
    }

    if (get.variance) {
        output$var.explained <- rowVars(scaled)
        output$var.total <- sum(output$var.explained)
    }

    output
}

###########################################

#' @importFrom S4Vectors List
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam
#' @importFrom Matrix crossprod t
#' @importFrom DelayedArray getAutoBPPARAM setAutoBPPARAM
.multi_pca_single <- function(mat, batch, subset.row=NULL, weights=NULL, d=50,
    get.all.genes=FALSE, get.variance=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam())
# Similar to .multi_pca_list, but avoids the unnecessary
# overhead of splitting 'mat' into batch-specific matrices
# when you end up having to put them back together again anyway.
{
    # NOTE: need to set this up here, as fastMNN() calls this function, not multiBatchPCA().
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

    FUN <- function(keep) {
        .process_single_matrix_for_pca(mat, batch=batch, weights=weights, subset.row=keep, deferred=bsdeferred(BSPARAM))
    }
    processed <- FUN(subset.row)
    centered <- processed$centered
    scaled <- processed$scaled
    centers <- processed$centers

    # Performing an SVD on the untransposed expression matrix,
    # and projecting the _unscaled_ matrices back into this space.
    if (!is.na(d)) {
        svd.out <- .run_scaled_SVD(scaled, d=d, get.all.genes=get.all.genes, get.variance=get.variance,
            BSPARAM=BSPARAM, BPPARAM=BPPARAM)

        output <- List(as.matrix(crossprod(centered, svd.out$u)))

        # Recording the rotation vectors.
        metadata(output) <- .make_pca_metadata(mat=mat, subset.row=subset.row,
            get.all.genes=get.all.genes, FUN=FUN, svd.out=svd.out, centers=centers, 
            get.variance=get.variance, nbatches=length(unique(batch)), scaled=scaled)
    } else {
        output <- List(as.matrix(t(centered)))
        metadata(output) <- .make_fake_metadata(mat=mat, subset.row=subset.row, 
            get.variance=get.variance, get.all.genes=get.all.genes, scaled=scaled)
    }

    output
}

#' @importFrom DelayedArray DelayedArray
#' @importFrom ScaledMatrix ScaledMatrix
#' @importFrom Matrix t rowMeans
.process_single_matrix_for_pca <- function(x, batch, weights, subset.row, deferred=FALSE) {
    batch <- factor(batch)
    tab <- table(batch)
    nms <- names(tab)
    tab <- as.numeric(tab) 
    names(tab) <- nms

    weights <- .construct_weight_vector(tab, weights)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }
    x <- realizeFileBackedMatrix(x)

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
        centered <- t(ScaledMatrix(t(x), center=centers))
        scaled <- ScaledMatrix(centered, scale=w)
    } else {
        centered <- DelayedArray(x)
        centered <- centered - centers
        scaled <- sweep(centered, 2, w, "/", check.margin=FALSE)
    }

    list(centered=centered, scaled=scaled, centers=centers)
}

###########################################

.set_deferred <- function(BSPARAM, deferred) {
    if (!is.null(deferred)) {
        # A bit naughty, but oh well.
        BSPARAM@deferred <- deferred 
    }
    BSPARAM
}
