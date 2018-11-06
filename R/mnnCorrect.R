#' Mutual nearest neighbors correction
#' 
#' Correct for batch effects in single-cell expression data using the mutual nearest neighbors method.
#' 
#' @param ... Two or more expression matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain cells from the same batch; multiple matrices represent separate batches of cells.
#' Each matrix should contain the same number of rows, corresponding to the same genes (in the same order).
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a log-expression matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for gene expression matrix inputs.
#' @param k An integer scalar specifying the number of nearest neighbors to consider when identifying mutual nearest neighbors.
#' @param sigma A numeric scalar specifying the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell.
#' @param cos.norm.in A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.
#' @param cos.norm.out A logical scalar indicating whether cosine normalization should be performed prior to computing corrected expression values.
#' @param svd.dim An integer scalar specifying the number of dimensions to use for summarizing biological substructure within each batch.
#' @param var.adj A logical scalar indicating whether variance adjustment should be performed on the correction vectors.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param correct.all A logical scalar specifying whether correction should be applied to all genes, even if only a subset is used for the MNN calculations.
#' @param order An integer vector specifying the order in which batches are to be corrected.
#' @param assay.type A string or integer scalar specifying the assay containing the expression values, if SingleCellExperiment objects are present in \code{...}.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the SVD algorithm to use.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbor search algorithm to use.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying the parallelization scheme to use.
#' 
#' @return
#' A \linkS4class{SummarizedExperiment} object containing the \code{corrected} assay.
#' This contains corrected expression values for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
#'
#' Cells (i.e., columns) are always ordered in the same manner as supplied in \code{...}, regardless of whether \code{order} is specified.
#' In cases with multiple objects in \code{...}, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#'
#' The metadata of the SummarizedExperiment contains:
#' \itemize{
#' \item{\code{pairs}: a list of DataFrames specifying which pairs of cells in \code{corrected} were identified as MNNs at each step.} 
#' \item{\code{order}: a vector of batch names or indices, specifying the order in which batches were merged.}
#' }
#' 
#' @details
#' This function is designed for batch correction of single-cell RNA-seq data where the batches are partially confounded with biological conditions of interest.
#' It does so by identifying pairs of mutual nearest neighbors (MNN) in the high-dimensional expression space.
#' Each MNN pair represents cells in different batches that are of the same cell type/state, assuming that batch effects are mostly orthogonal to the biological manifold.
#' Correction vectors are calculated from the pairs of MNNs and corrected expression values are returned for use in clustering and dimensionality reduction.
#' 
#' The threshold to define nearest neighbors is defined by \code{k}, which is passed to \code{\link{findMutualNN}} to identify MNN pairs.
#' The size of \code{k} can be interpreted as the minimum size of a subpopulation in each batch.
#' Values that are too small will not yield enough MNN pairs, while values that are too large will ignore substructure within each batch.
#' The algorithm is generally robust to various choices of \code{k}.
#' 
#' For each MNN pair, a pairwise correction vector is computed based on the difference in the expression profiles.
#' The correction vector for each cell is computed by applying a Gaussian smoothing kernel with bandwidth \code{sigma} is the pairwise vectors.
#' This stabilizes the vectors across many MNN pairs and extends the correction to those cells that do not have MNNs.
#' The choice of \code{sigma} determines the extent of smoothing - a value of 0.1 is used by default, corresponding to 10\% of the radius of the space after cosine normalization.
#' 
#' % We would consider 20 cells involved in MNN pairs to be the minimum number required for stable batch correction.
#'
#' @section Choosing the gene set:
#' Distances between cells are calculated with all genes if \code{subset.row=NULL}.
#' However, users can set \code{subset.row} to perform the distance calculation on a subset of genes, e.g., highly variable genes or marker genes.
#' This may provide more meaningful identification of MNN pairs by reducing the noise from irrelevant genes.
#' 
#' If \code{subset.row} is specified and \code{correct.all=TRUE}, corrected values are returned for \emph{all} genes.
#' This is possible as \code{subset.row} is only used to identify the MNN pairs and other cell-based distance calculations.
#' Correction vectors between MNN pairs can then be computed in the original space involving all genes in the supplied matrices.
#' 
#' @section Expected type of input data:
#' The input expression values should generally be log-transformed, e.g., log-counts, see \code{\link[scater]{normalize}} for details.
#' They should also be normalized within each data set to remove cell-specific biases in capture efficiency and sequencing depth.
#' By default, a further cosine normalization step is performed on the supplied expression data to eliminate gross scaling differences between data sets.
#' \itemize{
#' \item When \code{cos.norm.in=TRUE}, cosine normalization is performed on the matrix of expression values used to compute distances between cells.
#' This can be turned off when there are no scaling differences between data sets. 
#' \item When \code{cos.norm.out=TRUE}, cosine normalization is performed on the matrix of values used to calculate correction vectors (and on which those vectors are applied).
#' This can be turned off to obtain corrected values on the log-scale, similar to the input data.
#' }
#' The cosine normalization is achieved using the \code{\link{cosineNorm}} function.
#' 
#' Users should note that the order in which batches are corrected will affect the final results.
#' The first batch in \code{order} is used as the reference batch against which the second batch is corrected.
#' Corrected values of the second batch are added to the reference batch, against which the third batch is corrected, and so on.
#' This strategy maximizes the chance of detecting sufficient MNN pairs for stable calculation of correction vectors in subsequent batches.
#' 
#' @section Further options:
#' The function depends on a shared biological manifold, i.e., one or more cell types/states being present in multiple batches.
#' If this is not true, MNNs may be incorrectly identified, resulting in over-correction and removal of interesting biology.
#' Some protection can be provided by removing components of the correction vectors that are parallel to the biological subspaces in each batch.
#' The biological subspace in each batch is identified with a SVD on the expression matrix to obtain \code{svd.dim} dimensions.
#' (By default, this option is turned off by setting \code{svd.dim=0}.)
#' 
#' If \code{var.adj=TRUE}, the function will adjust the correction vector to equalize the variances of the two data sets along the batch effect vector.
#' In particular, it avoids \dQuote{kissing} effects whereby MNN pairs are identified between the surfaces of point clouds from different batches.
#' Naive correction would then bring only the surfaces into contact, rather than fully merging the clouds together.
#' The adjustment ensures that the cells from the two batches are properly intermingled after correction.
#' This is done by identifying each cell's position on the correction vector, identifying corresponding quantiles between batches, 
#' and scaling the correction vector to ensure that the quantiles are matched after correction.
#' 
#' @author
#' Laleh Haghverdi,
#' with modifications by Aaron Lun
#' 
#' @seealso
#' \code{\link{fastMNN}} for a faster equivalent.
#' 
#' @examples
#' B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
#' out <- mnnCorrect(B1, B2) # corrected values
#' 
#' @references
#' Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018).
#' Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.
#' \emph{Nat. Biotechnol.} 36(5):421
#' 
#' @export
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors metadata metadata<-
mnnCorrect <- function(..., batch=NULL, k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, svd.dim=0L, var.adj=TRUE, 
    subset.row=NULL, correct.all=FALSE, order=NULL, 
    assay.type="logcounts", get.spikes=FALSE, BSPARAM=NULL, BNPARAM=NULL, BPPARAM=SerialParam())
{
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

    output <- do.call(.mnn_correct, c(batches, 
        list(k=k, sigma=sigma, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, svd.dim=svd.dim, 
            var.adj=var.adj, subset.row=subset.row, correct.all=correct.all, order=order, 
            BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)))

    # Reordering the output for correctness.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
        metadata(output)$pairs <- .reindex_pairings(metadata(output)$pairs, d.reo)
    }

    return(output)
}

####################################
# Internal main function, to separate data handling from the actual calculations.

#' @importFrom S4Vectors DataFrame 
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocGenerics t rbind
#' @importFrom DelayedArray DelayedArray
#' @importFrom SummarizedExperiment SummarizedExperiment
.mnn_correct <- function(..., k=20, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, svd.dim=0L, var.adj=TRUE, 
    subset.row=NULL, correct.all=FALSE, order=NULL, BSPARAM=NULL, BNPARAM=NULL, BPPARAM=SerialParam())
{
    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
    
    prep.out <- .prepare_input_data(batches, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, subset.row=subset.row, correct.all=correct.all)
    in.batches <- prep.out$In
    out.batches <- prep.out$Out
    subset.row <- prep.out$Subset
    same.set <- prep.out$Same

    # Setting up the order.
    use.order <- !is.null(order)
    if (!use.order) {
        order <- seq_len(nbatches)
    } else {
        order <- as.integer(order)
        if (!identical(seq_len(nbatches), sort(order))) { 
            stop(sprintf("'order' should contain values in 1:%i", nbatches))
        }
    }
   
    # Setting up the variables.
    ref <- order[1]
    ref.batch.in <- t(in.batches[[ref]])
    if (!same.set) { 
        ref.batch.out <- t(out.batches[[ref]])
    }
    mnn.pairings <- vector("list", nbatches-1L)

    # Looping through the batches.
    for (b in 2:nbatches) { 
        target <- order[b]
        other.batch.in.untrans <- in.batches[[target]]
        other.batch.in <- t(other.batch.in.untrans)
        if (!same.set) { 
            other.batch.out.untrans <- out.batches[[target]]
            other.batch.out <- t(other.batch.out.untrans)
        }
        
        # Finding pairs of mutual nearest neighbours.
        sets <- findMutualNN(ref.batch.in, other.batch.in, k1=k, k2=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        s1 <- sets$first
        s2 <- sets$second      
        mnn.pairings[[b-1L]] <- DataFrame(first=s1, second=s2 + nrow(ref.batch.in))

        # Computing the correction vector.
        correction.in <- .compute_correction_vectors(ref.batch.in, other.batch.in, s1, s2, other.batch.in.untrans, sigma)
        if (!same.set) {
            correction.out <- .compute_correction_vectors(ref.batch.out, other.batch.out, s1, s2, other.batch.in.untrans, sigma)
            # NOTE: use of 'other.batch.in.untrans' here is intentional, 
            # as the distances should be the same as the MNN distances.
        }

        # Removing any component of the correction vector that's parallel to the biological basis vectors in either batch.
        # Only using cells involved in MNN pairs, to avoid undercorrection in directions that weren't problematic anyway
        # (given that using SVDs was intended to mitigate the effect of identifying the wrong MNN pairs).
        if (svd.dim>0) {
            u1 <- unique(s1)
            u2 <- unique(s2)

            # Computing the biological subspace in both batches, and subtract it from the batch correction vector.
            in.span1 <- .get_bio_span(t(ref.batch.in[u1,,drop=FALSE]), ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
            in.span2 <- .get_bio_span(other.batch.in.untrans[,u2,drop=FALSE], ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
            correction.in <- .subtract_bio(correction.in, in.span1, in.span2)

            # Repeating for the output values.
            if (!same.set) { 
                out.span1 <- .get_bio_span(t(ref.batch.out[u1,,drop=FALSE]), subset.row=subset.row, ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
                out.span2 <- .get_bio_span(other.batch.out.untrans[,u2,drop=FALSE], subset.row=subset.row, ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
                correction.out <- .subtract_bio(correction.out, out.span1, out.span2, subset.row=subset.row)
            }
        } 
       
        # Adjusting the shift variance; done after any SVD so that variance along the correction vector is purely technical.
        if (var.adj) { 
            correction.in <- .adjust_shift_variance(ref.batch.in, other.batch.in, correction.in, sigma=sigma)
            if (!same.set) {
                correction.out <- .adjust_shift_variance(ref.batch.out, other.batch.out, correction.out, sigma=sigma, subset.row=subset.row) 
            }
        }

        # Applying the correction and expanding the reference batch. 
        other.batch.in <- other.batch.in + correction.in
        ref.batch.in <- rbind(ref.batch.in, other.batch.in)
        if (!same.set) {
            other.batch.out <- other.batch.out + correction.out
            ref.batch.out <- rbind(ref.batch.out, other.batch.out)
        }
    }

    # Formatting the output.
    if (same.set) {
        ref.batch.out <- ref.batch.in
    }

    ncells.per.batch <- vapply(batches, FUN=ncol, FUN.VALUE=0L)
    batch.names <- .create_batch_names(names(batches), ncells.per.batch)

    # Adjusting the output back to the input order in '...'.
    if (use.order) {
        ordering <- .restore_original_order(order, ncells.per.batch)
        ref.batch.out <- ref.batch.out[ordering,,drop=FALSE]
        mnn.pairings <- .reindex_pairings(mnn.pairings, ordering)
    }
   
	SummarizedExperiment(list(corrected=t(ref.batch.out)), colData=DataFrame(batch=batch.names$ids),
        metadata=list(pairs=mnn.pairings, order=batch.names$labels[order]))
}

####################################
# Input/output functions.

#' @importFrom S4Vectors normalizeSingleBracketSubscript
.prepare_input_data <- function(batches, cos.norm.in, cos.norm.out, subset.row, correct.all) {
    nbatches <- length(batches)

    # Checking for identical number of rows (and rownames).
    first <- batches[[1]]
    ref.nrow <- nrow(first)
    ref.rownames <- rownames(first)
    for (b in 2:nbatches) {
        current <- batches[[b]]
        if (!identical(nrow(current), ref.nrow)) {
            stop("number of rows is not the same across batches")
        } else if (!identical(rownames(current), ref.rownames)) {
            stop("row names are not the same across batches")
        }
    }

    # Subsetting to the desired subset of genes.
    in.batches <- out.batches <- batches
    same.set <- TRUE
    if (!is.null(subset.row)) { 
        subset.row <- normalizeSingleBracketSubscript(subset.row, batches[[1]])
        if (identical(subset.row, seq_len(ref.nrow))) { 
            subset.row <- NULL
        } else {
            in.batches <- lapply(in.batches, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
            if (correct.all) {
                same.set <- FALSE
            } else {
                out.batches <- in.batches
            }
        }
    }

    # Applying cosine normalization for MNN identification. 
    # We use the L2 norm for the subsetted input to adjust the output, 
    # to ensure that results are consistent regardless of the manner of subsetting.
    if (cos.norm.in) { 
        norm.scaling <- vector("list", nbatches)
        for (b in seq_len(nbatches)) { 
            current.in <- in.batches[[b]]
            cos.out <- cosineNorm(current.in, mode="all")
            in.batches[[b]] <- cos.out$matrix
            norm.scaling[[b]] <- cos.out$l2norm
        }
    }
    if (cos.norm.out) { 
        if (!cos.norm.in) { 
            norm.scaling <- lapply(in.batches, cosineNorm, mode="l2norm")
        }
        for (b in seq_len(nbatches)) { 
            out.batches[[b]] <- sweep(out.batches[[b]], 2, norm.scaling[[b]], "/")
        }
    }
    if (cos.norm.out!=cos.norm.in) { 
        same.set <- FALSE
    }

    return(list(In=in.batches, Out=out.batches, Subset=subset.row, Same=same.set))
}

####################################
# Internal functions for correction.

.compute_correction_vectors <- function(data1, data2, mnn1, mnn2, original2, sigma) 
# Computes the batch correction vector for each cell in 'data2'.
# 'original2' should also be supplied to compute distances 
# (this may not be the same as 't(data2)' due to normalization, subsetting).
{      
    vect <- data1[mnn1,,drop=FALSE] - data2[mnn2,,drop=FALSE]    
    cell.vect <- .Call(cxx_smooth_gaussian_kernel, vect, mnn2-1L, original2, sigma)
    return(t(cell.vect)) 
}

.adjust_shift_variance <- function(data1, data2, correction, sigma, subset.row=NULL) 
# Performs variance adjustment to avoid kissing effects.    
{
    cell.vect <- correction 
    if (!is.null(subset.row)) { 
        # Only using subsetted genes to compute locations, consistent with our policy in SVD. 
        cell.vect <- cell.vect[,subset.row,drop=FALSE]
        data1 <- data1[,subset.row,drop=FALSE]
        data2 <- data2[,subset.row,drop=FALSE]
    }
    scaling <- .Call(cxx_adjust_shift_variance, t(data1), t(data2), cell.vect, sigma)
    scaling <- pmax(scaling, 1)
    return(scaling * correction)
}

#' @importFrom BiocSingular runSVD
#' @importFrom BiocGenerics rowMeans
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam
.get_bio_span <- function(exprs, ndim, subset.row=NULL, BSPARAM=NULL, BPPARAM=SerialParam())
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
{
    centered <- DelayedArray(exprs) - rowMeans(exprs)

    if (!is.null(subset.row)) {
        leftovers <- centered[-subset.row,,drop=FALSE]
        centered <- centered[subset.row,,drop=FALSE]
        nv <- ndim
    } else {
        nv <- 0
    }

    # Returning the basis vectors directly if there was no subsetting. 
    ndim <- min(ndim, dim(centered))
    nv <- min(nv, ndim)
    S <- runSVD(centered, k=ndim, nu=ndim, nv=nv, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    if (is.null(subset.row)) { 
        return(S$u)
    } 

    # Computing the basis values for the unused genes.
    output <- matrix(0, nrow(exprs), ndim)
    output[subset.row,] <- S$u
    leftovers <- as.matrix(leftovers %*% S$v)
    leftovers <- sweep(leftovers, 2, S$d, "/")
    output[-subset.row,] <- leftovers
    return(output)
}

.subtract_bio <- function(correction, span1, span2, subset.row=NULL) 
# Computes the component parallel biological basis vectors in the correction vectors,
# and subtracts them. Note that the order of span1 and span2 does not matter.
{ 
    for (span in list(span1, span2)) { 
        if (is.null(subset.row)) { 
            bio.mag <- correction %*% span
        } else { 
            # Only using the subset to compute the magnitude that is parallel.
            bio.mag <- correction[,subset.row,drop=FALSE] %*% span[subset.row,,drop=FALSE]
        }
        bio.comp <- bio.mag %*% t(span)
        correction <- correction - bio.comp
    }
    return(correction)
}  
