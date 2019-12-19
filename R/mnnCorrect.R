#' Mutual nearest neighbors correction
#' 
#' Correct for batch effects in single-cell expression data using the mutual nearest neighbors method.
#' 
#' @param sigma A numeric scalar specifying the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell.
#' @param cos.norm.in A logical scalar indicating whether cosine normalization should be performed on the input data prior to calculating distances between cells.
#' @param cos.norm.out A logical scalar indicating whether cosine normalization should be performed prior to computing corrected expression values.
#' @param svd.dim An integer scalar specifying the number of dimensions to use for summarizing biological substructure within each batch.
#' @param var.adj A logical scalar indicating whether variance adjustment should be performed on the correction vectors.
#' @param correct.all A logical scalar specifying whether correction should be applied to all genes, even if only a subset is used for the MNN calculations.
#' @inheritParams fastMNN
#' 
#' @return
#' A \linkS4class{SingleCellExperiment} object containing the \code{corrected} assay.
#' This contains corrected expression values for each gene (row) in each cell (column) in each batch.
#' A \code{batch} field is present in the column data, specifying the batch of origin for each cell.
#'
#' Cells in the output object are always ordered in the same manner as supplied in \code{...}.
#' For a single input object, cells will be reported in the same order as they are arranged in that object.
#' In cases with multiple input objects, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#'
#' The metadata of the SingleCellExperiment contains \code{merge.info}, a DataFrame where each row corresponds to a merge step.
#' See \dQuote{Merge diagnostics} for more information.
#'
#' @details
#' This function is designed for batch correction of single-cell RNA-seq data where the batches are partially confounded with biological conditions of interest.
#' It does so by identifying pairs of mutual nearest neighbors (MNN) in the high-dimensional log-expression space.
#' Each MNN pair represents cells in different batches that are of the same cell type/state, assuming that batch effects are mostly orthogonal to the biological manifold.
#' Correction vectors are calculated from the pairs of MNNs and corrected (log-)expression values are returned for use in clustering and dimensionality reduction.
#' 
#' For each MNN pair, a pairwise correction vector is computed based on the difference in the log-expression profiles.
#' The correction vector for each cell is computed by applying a Gaussian smoothing kernel with bandwidth \code{sigma} is the pairwise vectors.
#' This stabilizes the vectors across many MNN pairs and extends the correction to those cells that do not have MNNs.
#' The choice of \code{sigma} determines the extent of smoothing - a value of 0.1 is used by default, corresponding to 10\% of the radius of the space after cosine normalization.
#' 
#' % We would consider 20 cells involved in MNN pairs to be the minimum number required for stable batch correction.
#'
#' @section Choosing the gene set:
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs to highly variable genes or marker genes.
#' This may provide more meaningful identification of MNN pairs by reducing the noise from irrelevant genes.
#' Note that users should not be too restrictive with subsetting, as high dimensionality is required to satisfy the orthogonality assumption in MNN detection.
#'
#' If \code{subset.row} is specified and \code{correct.all=TRUE}, corrected values are returned for \emph{all} genes.
#' This is possible as \code{subset.row} is only used to identify the MNN pairs and other cell-based distance calculations.
#' Correction vectors between MNN pairs can then be computed in for all genes in the supplied matrices.
#' Note that setting \code{correct.all=TRUE} will not alter the corrected expression values for the subsetted genes.
#' 
#' @section Expected type of input data:
#' The input expression values should generally be log-transformed, e.g., log-counts, see \code{\link[scater]{logNormCounts}} for details.
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
#' @inheritSection fastMNN Specifying the number of neighbors
#' 
#' @inheritSection fastMNN Controlling the merge order
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
#' See \code{?"\link{batchelor-restrict}"} for a description of the \code{restrict} argument.
#' Specifically, \code{mnnCorrect} will only use the restricted subset of cells in each batch to identify MNN pairs (and to perform variance adjustment, if \code{var.adj=TRUE}), and then apply the correction to all cells in each batch. 
#'
#' @section Merge diagnostics:
#' Each merge step combines two mutually exclusive sets of cells, a \dQuote{left} set and \dQuote{right} set.
#' The metadata thus contains the following fields:
#' \itemize{
#' \item \code{left}, a \linkS4class{List} of integer or character vectors.
#' Each vector specifies the batches in the left set at a given merge step. 
#' \item \code{right}, a similar List of integer or character vectors.
#' Each vector specifies the batches in the right set at a given merge step. 
#' \item \code{pairs}, a List of DataFrames specifying which pairs of cells were identified as MNNs at each step.
#' In each DataFrame, each row corresponds to a single MNN pair and specifies the
#' paired cells that were in the left and right sets, respectively.
#' Note that the indices refer to those paired cells in the \emph{output} ordering of cells,
#' i.e., users can identify the paired cells at each step by column-indexing the output of the \code{mnnCorrect} function.
#' }
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
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocNeighbors KmknnParam
mnnCorrect <- function(..., batch=NULL, restrict=NULL, k=20, prop.k=NULL, sigma=0.1, 
    cos.norm.in=TRUE, cos.norm.out=TRUE, svd.dim=0L, var.adj=TRUE, 
    subset.row=NULL, correct.all=FALSE, merge.order=NULL, auto.merge=FALSE, 
    assay.type="logcounts", BSPARAM=ExactParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam())
{
    original <- batches <- .unpack_batches(...)
    checkBatchConsistency(batches)
    restrict <- checkRestrictions(batches, restrict)
    
    # Pulling out information from the SCE objects.        
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type, withDimnames=FALSE)
    }

    # Subsetting by 'batch'.
    do.split <- length(batches)==1L
    if (do.split) {
        divided <- divideIntoBatches(batches[[1]], batch=batch, restrict=restrict[[1]])
        batches <- divided$batches
        restrict <- divided$restrict
    }

    # Setting up the parallelization environment.
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    output <- do.call(.mnn_correct, c(batches, 
        list(k=k, prop.k=prop.k, sigma=sigma, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, svd.dim=svd.dim, 
            var.adj=var.adj, subset.row=subset.row, correct.all=correct.all, restrict=restrict,
            merge.order=merge.order, auto.merge=auto.merge, 
            BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)))

    # Reordering the output for correctness.
    if (do.split) {
        d.reo <- divided$reorder
        output <- output[,d.reo,drop=FALSE]
        metadata(output)$merge.info$pairs <- .reindex_pairings(metadata(output)$merge.info$pairs, d.reo)
    }

    .rename_output(output, original, subset.row=subset.row)
}

####################################
# Internal main functions, to separate input/output munging from the actual calculations.
# This is split into two functions - .mnn_correct() sets up the core function, which
# defines how to handle predefined vs automatically determined merge trees.

#' @importFrom Matrix t
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocNeighbors KmknnParam
.mnn_correct <- function(..., k=20, prop.k=NULL, sigma=0.1, cos.norm.in=TRUE, cos.norm.out=TRUE, svd.dim=0L, var.adj=TRUE, 
    subset.row=NULL, correct.all=FALSE, restrict=NULL, 
    merge.order=NULL, auto.merge=FALSE, 
    BSPARAM=ExactParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam())
{
    batches <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }
    
    # Setting up the variables.
    prep.out <- .prepare_input_data(batches, cos.norm.in=cos.norm.in, cos.norm.out=cos.norm.out, 
        subset.row=subset.row, correct.all=correct.all)
    in.batches <- prep.out$In
    out.batches <- prep.out$Out
    subset.row <- prep.out$Subset
    same.set <- prep.out$Same

    in.batches <- lapply(in.batches, t)
    if (!same.set) {
        out.batches <- lapply(out.batches, t)
    }

    # Choosing the merge order.
    if (!auto.merge) {
        merge.tree <- .create_tree_predefined(in.batches, restrict, merge.order)
        NEXT <- .get_next_merge
        UPDATE <- .update_tree

        # Putting the out.batches information in the 'extras'.
        merge.tree <- .add_out_batches_to_tree(merge.tree, if (same.set) NULL else out.batches)
    } else {
        mnn.args <- list(k=k, prop.k=prop.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM, orthogonalize=FALSE)
        merge.tree <- do.call(.initialize_auto_search, c(list(in.batches, restrict), mnn.args))

        NEXT <- .pick_best_merge
        UPDATE <- function(remainders, chosen, ...) {
            .update_remainders(remainders, chosen, ..., mnn.args=mnn.args)
        }

        # Putting the out.batches information in the 'extras'.
        for (i in seq_along(merge.tree)) {
            merge.tree[[i]]@extras <- if (same.set) list(NULL) else out.batches[.get_node_index(merge.tree[[i]])]
        }
    }
    output <- .mnn_correct_core(merge.tree, NEXT=NEXT, UPDATE=UPDATE, k=k, prop.k=prop.k, sigma=sigma,
        same.set=same.set, svd.dim=svd.dim, var.adj=var.adj, subset.row=subset.row, 
        BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    nms <- names(batches)
    if (!is.null(nms)) {
        if (anyDuplicated(nms)) {
            stop("names of batches should be unique")
        }
        output$batch <- nms[output$batch]
        output <- .fix_names_in_merge_info(output, nms)
    }
    output
}

.add_out_batches_to_tree <- function(merge.tree, out.batches) {
    if (!is.list(merge.tree)) {
        merge.tree@extras <- list(out.batches[[.get_node_index(merge.tree)]])
        return(merge.tree)
    }
    merge.tree[[1]] <- .add_out_batches_to_tree(merge.tree[[1]], out.batches)
    merge.tree[[2]] <- .add_out_batches_to_tree(merge.tree[[2]], out.batches)
    merge.tree
}

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular ExactParam
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom S4Vectors DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment 
#' @importFrom Matrix t
.mnn_correct_core <- function(merge.tree, NEXT, UPDATE,
    same.set=FALSE, k=20, prop.k=NULL, sigma=0.1, svd.dim=0L, var.adj=TRUE, subset.row=NULL,
    BSPARAM=ExactParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    nbatches <- length(unlist(merge.tree))
    nmerges <- nbatches - 1L
    mnn.pairings <- left.set <- right.set <- vector("list", nmerges)

    for (mdx in seq_len(nmerges)) {
        # Traversing the merge tree to find the next two batches to merge.
        next.step <- NEXT(merge.tree) 
        left <- next.step$left
        right <- next.step$right

        left.data <- .get_node_data(left)
        left.restrict <- .get_node_restrict(left)
        left.index <- .get_node_index(left)
        left.origin <- .get_node_origin(left)
        left.extras <- .get_node_extras(left)[[1]]

        right.data <- .get_node_data(right)
        right.restrict <- .get_node_restrict(right)
        right.index <- .get_node_index(right)
        right.origin <- .get_node_origin(right)
        right.extras <- .get_node_extras(right)[[1]]

        # Computing the correction vector for each cell.
        mnn.sets <- .restricted_mnn(left.data, left.restrict, right.data, right.restrict,
            k=k, prop.k=prop.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        s1 <- mnn.sets$first
        s2 <- mnn.sets$second

        mnn.pairings[[mdx]] <- DataFrame(left=s1, right=s2)
        left.set[[mdx]] <- left.index
        right.set[[mdx]] <- right.index

        trans.right <- t(right.data)
        correction.in <- .compute_correction_vectors(left.data, right.data, s1, s2, trans.right, sigma)
        if (!same.set) {
            # NOTE: use of 'trans.right' here is intentional, as the distances
            # used for the kernel weighting should be comparable to 'sigma';
            # this protects against the case where cos.in=TRUE and cos.out=FALSE.
            correction.out <- .compute_correction_vectors(left.extras, right.extras, s1, s2, trans.right, sigma)
        }

        # Removing any component of the correction vector that's parallel to
        # the biological basis vectors in either batch.  Only using cells
        # involved in MNN pairs, to avoid undercorrection in directions that
        # weren't problematic anyway (given that using SVDs was intended to
        # mitigate the effect of identifying the wrong MNN pairs).
        if (svd.dim>0) {
            u1 <- unique(s1)
            u2 <- unique(s2)

            # Computing the biological subspace in both batches, and subtract it from the batch correction vector.
            in.span1 <- .get_bio_span(t(left.data[u1,,drop=FALSE]), ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
            in.span2 <- .get_bio_span(t(right.data[u2,,drop=FALSE]), ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
            correction.in <- .subtract_bio(correction.in, in.span1, in.span2)

            if (!same.set) { 
                out.span1 <- .get_bio_span(t(left.extras[u1,,drop=FALSE]), subset.row=subset.row, 
                    ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
                out.span2 <- .get_bio_span(t(right.extras[u2,,drop=FALSE]), subset.row=subset.row, 
                    ndim=svd.dim, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
                correction.out <- .subtract_bio(correction.out, out.span1, out.span2, subset.row=subset.row)
            }
        } 
       
        # Adjusting the shift variance; done after any SVD so that variance
        # along the correction vector is purely technical.
        if (var.adj) {
            args <- list(sigma=sigma, restrict1=left.restrict, restrict2=right.restrict)

            correction.in <- do.call(.adjust_shift_variance, 
                c(list(t(left.data), t(right.data), correction.in), args)
            ) 
            if (!same.set) {
                correction.out <- do.call(.adjust_shift_variance, 
                    c(list(t(left.extras), t(right.extras), correction.out, subset.row=subset.row), args)
                )
            }
        }

        # Applying the correction and expanding the reference batch. 
        right.data <- right.data + correction.in
        if (!same.set) {
            right.extras <- right.extras + correction.out
        }

        merge.tree <- UPDATE(merge.tree, next.step$chosen, 
            data=rbind(left.data, right.data),  
            index=c(left.index, right.index),
            restrict=.combine_restrict(left.data, left.restrict, right.data, right.restrict),
            origin=c(left.origin, right.origin),
            extras=list(rbind(left.extras, right.extras)))
    }

    full.order <- .get_node_index(merge.tree)
    full.origin <- .get_node_origin(merge.tree)
    if (same.set) {
        full.data <- .get_node_data(merge.tree)
    } else {
        full.data <- .get_node_extras(merge.tree)[[1]]
    }

    # Re-indexing all of the pairing indices to account for the final position of each cell.
    for (mdx in seq_along(mnn.pairings)) {
        bonus1 <- match(left.set[[mdx]][1], full.origin) - 1L
        mnn.pairings[[mdx]]$left <- mnn.pairings[[mdx]]$left + bonus1
        bonus2 <- match(right.set[[mdx]][1], full.origin) - 1L
        mnn.pairings[[mdx]]$right <- mnn.pairings[[mdx]]$right + bonus2
    }

    # Adjusting the output back to the input order in 'batches'.
    if (is.unsorted(full.order)) {
        ncells.per.batch <- tabulate(full.origin)
        ordering <- .restore_original_order(full.order, ncells.per.batch)
        full.data <- full.data[ordering,,drop=FALSE]
        full.origin <- full.origin[ordering]
        mnn.pairings <- .reindex_pairings(mnn.pairings, ordering)
    }

	SingleCellExperiment(list(corrected=t(full.data)), 
        colData=DataFrame(batch=full.origin),
        metadata=list(
            merge.info=DataFrame(
                left=I(as(left.set, "List")),
                right=I(as(right.set, "List")),
                pairs=I(as(mnn.pairings, "List"))
            ) 
        )
    )
}

####################################
# Input/output functions.

.prepare_input_data <- function(batches, cos.norm.in, cos.norm.out, subset.row, correct.all) {
    nbatches <- length(batches)
    in.batches <- out.batches <- batches
    same.set <- TRUE

    # Subsetting to the desired subset of genes.
    if (!is.null(subset.row)) { 
        subset.row <- .row_subset_to_index(batches[[1]], subset.row)
        if (identical(subset.row, seq_len(nrow(batches[[1]])))) { 
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
        out.batches <- mapply(FUN=.apply_cosine_norm, out.batches, norm.scaling, SIMPLIFY=FALSE)
    }

    if (cos.norm.out!=cos.norm.in) { 
        same.set <- FALSE
    }

    list(In=in.batches, Out=out.batches, Subset=subset.row, Same=same.set)
}

####################################
# Internal functions for correction.

#' @importFrom Matrix t
.compute_correction_vectors <- function(data1, data2, mnn1, mnn2, tdata2, sigma) 
# Computes the batch correction vector for each cell in 'data2'.
# 'tdata2' should also be supplied to compute distances 
# (this may not be the same as 't(data2)' due to normalization, subsetting).
{      
    vect <- data1[mnn1,,drop=FALSE] - data2[mnn2,,drop=FALSE]    
    cell.vect <- .Call(cxx_smooth_gaussian_kernel, vect, mnn2-1L, tdata2, sigma)
    t(cell.vect)
}

.adjust_shift_variance <- function(data1, data2, correction, sigma, subset.row=NULL, restrict1=NULL, restrict2=NULL) 
# Performs variance adjustment to avoid kissing effects.    
{
    cell.vect <- correction 

    if (!is.null(subset.row)) { 
        # Only using subsetted genes to compute locations, consistent with our policy in SVD. 
        cell.vect <- cell.vect[,subset.row,drop=FALSE]
        data1 <- data1[subset.row,,drop=FALSE]
        data2 <- data2[subset.row,,drop=FALSE]
    }

    restrict1 <- .col_subset_to_index(data1, restrict1) - 1L
    restrict2 <- .col_subset_to_index(data2, restrict2) - 1L
    scaling <- .Call(cxx_adjust_shift_variance, data1, data2, cell.vect, sigma, restrict1, restrict2)

    scaling <- pmax(scaling, 1)
    scaling * correction
}

#' @importFrom BiocSingular runSVD ExactParam
#' @importFrom Matrix rowMeans
#' @importFrom DelayedArray DelayedArray getAutoBPPARAM setAutoBPPARAM
#' @importFrom BiocParallel SerialParam
.get_bio_span <- function(exprs, ndim, subset.row=NULL, BSPARAM=ExactParam(), BPPARAM=SerialParam())
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
{
    # Protect against DA parallelization.
    old <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(old))

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
