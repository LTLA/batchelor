#' Fast mutual nearest neighbors correction
#'
#' Correct for batch effects in single-cell expression data using a fast version of the mutual nearest neighbors (MNN) method.
#'
#' @param ... One or more log-expression matrix where genes correspond to rows and cells correspond to columns, if \code{pc.input=FALSE}.
#' Each matrix should contain the same number of rows, corresponding to the same genes in the same order.
#' 
#' Alternatively, one or more matrices of low-dimensional representations can be supplied if \code{pc.input=TRUE}, where rows are cells and columns are dimensions.
#' Each object should contain the same number of columns, corresponding to the same dimensions.
#'
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a log-expression matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for matrix inputs.
#'
#' Alternatively, the SingleCellExperiment objects can contain reduced dimension coordinates in the \code{reducedDims} slot if \code{use.dimred} is specified.
#' Note the same restrictions described above for low-dimensional matrix inputs.
#' 
#' In all cases, each object contains cells from a single batch; multiple objects represent separate batches of cells.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param restrict A list of length equal to the number of objects in \code{...}.
#' Each entry of the list corresponds to one batch and specifies the cells to use when computing the correction.
#' @param k An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.
#' @param cos.norm A logical scalar indicating whether cosine normalization should be performed on the input data prior to PCA.
#' @param ndist A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors.
#' Each threshold is defined as a multiple of the number of median distances.
#' @param d Number of dimensions to use for dimensionality reduction in \code{\link{multiBatchPCA}}.
#' @param auto.order Logical scalar indicating whether re-ordering of batches should be performed to maximize the number of MNN pairs at each step.
#' 
#' Alternatively, an integer vector containing a permutation of \code{1:N} where \code{N} is the number of batches.
#' @param compute.variances Logical scalar indicating whether the percentage of variance lost due to non-orthogonality should be computed.
#' @param subset.row A vector specifying which features to use for correction. 
#' Only relevant for gene expression inputs (i.e., \code{pc.input=FALSE} and \code{use.dimred=NULL}).
#' @param correct.all Logical scalar indicating whether a rotation matrix should be computed for genes not in \code{subset.row}.
#' Only used for gene expression inputs, i.e., when \code{pc.input=FALSE}.
#' @param pc.input Logical scalar indicating whether the values in \code{...} are already low-dimensional, e.g., the output of \code{\link{multiBatchPCA}}.
#' Only used when \code{...} does \emph{not} contain SingleCellExperiment objects - in those cases, set \code{use.dimred} instead.
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values.
#' Only used for SingleCellExperiment inputs with \code{use.dimred=NULL}.
#' @param get.spikes A logical scalar indicating whether to retain rows corresponding to spike-in transcripts.
#' Only used for SingleCellExperiment inputs with \code{use.dimred=NULL}.
#' @param use.dimred A string or integer scalar specifying which reduced dimension result to use, if any.
#' Only used for SingleCellExperiment inputs.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA and nearest-neighbor searches should be parallelized.
#' 
#' @return
#' The output of this function depends on whether a PCA is performed on the input \code{...}.
#' This will be the case if \code{pc.input=FALSE} for matrix inputs or if \code{use.dimred=NULL} for SingleCellExperiment inputs.
#' 
#' If a PCA is performed, a \linkS4class{SingleCellExperiment} is returned where each row is a gene and each column is a cell. 
#' This contains:
#' \itemize{
#' \item A \code{corrected} matrix in the \code{reducedDims} slot, containing corrected low-dimensional coordinates for each cell.
#' This has number of columns equal to \code{d} and number of rows equal to the total number of cells in \code{...}.
#' \item A \code{batch} column in the \code{colData} slot, containing the batch of origin for each row (i.e., cell) in \code{corrected}.
#' \item A \code{rotation} column the \code{rowData} slot, containing the rotation matrix used for the PCA.
#' This has \code{d} columns and number of rows equal to the number of genes to report (see the \dQuote{Choice of genes} section).
#' \item A \code{reconstructed} matrix in the \code{assays} slot, containing the low-rank reconstruction of the original expression matrix.
#' This can be interpreted as per-gene corrected log-expression values (after cosine normalization, if \code{cos.norm=TRUE}) but should not be used for quantitative analyses.
#' }
#' 
#' Otherwise, a \linkS4class{DataFrame} is returned where each row corresponds to a cell.
#' It contains:
#' \itemize{
#' \item \code{corrected}, the matrix of corrected low-dimensional coordinates for each cell.
#' \item \code{batch}, the Rle specifying the batch of origin for each row. 
#' }
#'
#' Cells in the output object are always ordered in the same manner as supplied in \code{...}.
#' For a single input object, cells will be reported in the same order as they are arranged in that object.
#' In cases with multiple input objects, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#' This is true regardless of the value of \code{auto.order}, which only affects the internal merge order of the batches.
#' 
#' The metadata of the output object contains:
#' \itemize{
#' \item \code{pairs}, a list of DataFrames specifying which pairs of cells in \code{corrected} were identified as MNNs at each step. 
#' \item \code{order}, a vector of batch names or indices, specifying the order in which batches were merged.
#' \item \code{lost.var}, a numeric vector containing the proportion of lost variance from each batch supplied in \code{...}.
#' Only returned when \code{compute.variances=TRUE}.
#' }
#' 
#' @details
#' This function provides a variant of the \code{\link{mnnCorrect}} function, modified for speed and more robust performance.
#' In particular:
#' \itemize{
#' \item It performs a multi-sample PCA via \code{\link{multiBatchPCA}} and subsequently performs all calculations in the PC space.
#' This reduces computational work and provides some denoising for improved neighbour detection. 
#' As a result, though, the corrected output cannot be interpreted on a gene level and is useful only for cell-level comparisons, e.g., clustering and visualization.
#' \item The correction vector for each cell is directly computed from its \code{k} nearest neighbours in the same batch.
#' Specifically, only the \code{k} nearest neighbouring cells that \emph{also} participate in MNN pairs are used.
#' Each MNN-participating neighbour is weighted by distance from the current cell, using a tricube scheme with bandwidth equal to the median distance multiplied by \code{ndist}.
#' This ensures that the correction vector only uses information from the closest cells, improving the fidelity of local correction.
#' \item Issues with \dQuote{kissing} are avoided with a two-step procedure that removes variation along the batch effect vector.
#' First, the average correction vector across all MNN pairs is computed.
#' Cell coordinates are adjusted such that all cells in a single batch have the same position along this vector.
#' The correction vectors are then recalculated with the adjusted coordinates (but the same MNN pairs).
#' }
#' 
#' The default setting of \code{cos.norm=TRUE} provides some protection against differences in scaling between log-expression matrices from batches that are normalized separately
#' (see \code{\link{cosineNorm}} for details).
#' However, if possible, we recommend using the output of \code{\link{multiBatchNorm}} as input to \code{fastMNN}.
#' This will equalize coverage on the count level before the log-transformation, which is a more accurate rescaling than cosine normalization on the log-values.
#' 
#' If \code{compute.variances=TRUE}, the function will compute the percentage of variance that is lost from each batch during orthogonalization.
#' This represents the variance in each batch that is parallel to the average correction vectors (and hence removed during orthogonalization) at each merge step.
#' Large proportions suggest that there is biological structure that is parallel to the batch effect, 
#' corresponding to violations of the assumption that the batch effect is orthogonal to the biological subspace.
#'
#' The \code{batch} argument allows users to easily perform batch correction when all cells have already been combined into a single object.
#' This avoids the need to manually split the matrix or SingleCellExperiment object into separate objects for input into \code{fastMNN}.
#' In this situation, the order of input batches is defined by the order of levels in \code{batch}.
#'
#' @section Controlling the merge order:
#' By default, batches are merged in the user-supplied order.
#' However, if \code{auto.order=TRUE}, batches are ordered to maximize the number of MNN pairs at each step.
#' The aim is to improve the stability of the correction by first merging more similar batches with more MNN pairs.
#' This can be somewhat time-consuming as MNN pairs need to be iteratively recomputed for all possible batch pairings.
#' It is often more convenient for the user to specify an appropriate ordering based on prior knowledge about the batches.
#' 
#' If \code{auto.order} is an integer vector, it is treated as an ordering permutation with which to merge batches.
#' For example, if \code{auto.order=c(4,1,3,2)}, batches 4 and 1 in \code{...} are merged first, followed by batch 3 and then batch 2.
#' This is often more convenient than changing the order manually in \code{...}, which would alter the order of batches in the output \code{corrected} matrix.
#' Indeed, no matter what the setting of \code{auto.order} is, the order of cells in the output corrected matrix is \emph{always} the same.
#' 
#' Further control of the merge order can be achieved by performing the multi-sample PCA outside of this function with \code{\link{multiBatchPCA}}.
#' Batches can then be progressively merged by repeated calls to \code{fastMNN} with \code{pc.input=TRUE} or by setting \code{use.dimred} (see \dQuote{Using low-dimensional inputs}).
#' This is useful in situations where the batches need to be merged in a hierarhical manner, e.g., combining replicate samples before merging them across different conditions.
#' For example, we could merge batch 1 with 4 to obtain a corrected 1+4; and then batch 2 with 3 to obtain a corrected 2+3;
#' before merging the corrected 1+4 and 2+3 to obtain the final set of corrected values.
#'
#' @section Choice of genes:
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs to highly variable genes or marker genes.
#' This improves the quality of the PCA and identification of MNN pairs by reducing the noise from irrelevant genes.
#' Note that users should not be too restrictive with subsetting, as high dimensionality is required to satisfy the orthogonality assumption in MNN detection.
#'
#' For \linkS4class{SingleCellExperiment} inputs, spike-in transcripts are automatically removed unless \code{get.spikes=TRUE}.
#' If \code{subset.row} is specified and \code{get.spikes=FALSE}, only the non-spike-in specified features will be used. 
#' All SingleCellExperiment objects should have the same set of spike-in transcripts.
#'
#' By default, only the selected genes are used to compute rotation vectors and a low-rank representation of the input matrix.
#' However, rotation vectors can be obtained that span all genes in the supplied input data with \code{correct.all=TRUE}.
#' This will not affect the corrected low-dimension coordinates or the output for the selected genes. 
#'
#' Note that these settings for the choice of genes are completely ignored when using low-dimensional inputs (see below).
#'
#' @section Using low-dimensional inputs:
#' Low-dimensional inputs can be supplied directly to \code{fastMNN} if the PCA (or some other projection to low-dimensional space) is performed outside the function.
#' This intructs the function to skip the \code{\link{multiBatchPCA}} step. 
#' To enable this, set \code{pc.input=TRUE} for matrix-like inputs in \code{...}, or specify \code{use.dimred} with SingleCellExperiment inputs.
#' We only recommend this mode for advanced users, and note that:
#' \itemize{
#'     \item \code{\link{multiBatchPCA}} will not perform cosine-normalization, 
#' so it is the responsibility of the user to cosine-normalize each batch beforehand with \code{\link{cosineNorm}} to recapitulate results with \code{cos.norm=TRUE}.
#'     \item \code{\link{multiBatchPCA}} must be run on all samples at once, to ensure that all cells are projected to the same low-dimensional space.
#'     \item Setting \code{pc.input=TRUE} is criticial to avoid unnecessary (and incorrect) cosine-normalization and PCA within each step of the merge.
#' }
#' Users are referred to the Examples for a demonstration of this functionality.
#'
#' @section Using restriction:
#' It is possible to compute the correction using only a subset of cells in each batch, and then extrapolate that correction to all other cells.
#' This may be desirable in experimental designs where a control set of cells from the same source population were run on different batches.
#' Any difference in the controls must be artificial in origin and can be directly removed without making further biological assumptions.
#'
#' To do this, users should set \code{restrict} to specify the subset of cells in each batch to be used for correction.
#' This should be set to a list of length equal to the length of \code{...}, where each element is a subsetting vector to be applied to the columns of the corresponding batch.
#' A \code{NULL} element indicates that all the cells from a batch should be used.
#' In situations where one input object contains multiple batches, \code{restrict} is simply a list containing a single subsetting vector for that object.
#' 
#' \code{fastMNN} will only use the restricted subset of cells in each batch to identify MNN pairs and the center of the orthogonalization.
#' However, it will apply the correction to all cells in each batch - hence the extrapolation.
#' This means that the output is always of the same dimensionality, regardless of whether \code{restrict} is specified.
#'
#' Note that \emph{all} cells are used to perform the PCA, regardless of whether \code{restrict} is set.
#' Constructing the projection vectors with only control cells will not guarantee resolution of unique non-control populations in each batch.
#' The function will only completely ignore cells that are not in \code{restrict} if \code{pc.input=TRUE} or, for SingleCellExperiment inputs, \code{use.dimred} is set.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{cosineNorm}} and \code{\link{multiBatchPCA}} to obtain the values to be corrected.
#'
#' \code{\link{mnnCorrect}} for the \dQuote{classic} version of the MNN correction algorithm.
#'
#' @examples
#' B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
#' out <- fastMNN(B1, B2)
#' str(reducedDim(out)) # corrected values
#' 
#' # An equivalent approach with PC input.
#' cB1 <- cosineNorm(B1)
#' cB2 <- cosineNorm(B2)
#' pcs <- multiBatchPCA(cB1, cB2)
#' out2 <- fastMNN(pcs[[1]], pcs[[2]], pc.input=TRUE)
#'
#' all.equal(out, out2) # should be TRUE
#' 
#' # Extracting corrected expression values for gene 10.
#' summary(assay(out)[10,])
#'
#' @references
#' Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018).
#' Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.
#' \emph{Nat. Biotechnol.} 36(5):421
#'
#' Lun ATL (2018).
#' Further MNN algorithm development.
#' \url{https://github.com/MarioniLab/FurtherMNN2018}
#'
#' @rdname fastMNN
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom BiocSingular ExactParam
fastMNN <- function(..., batch=NULL, k=20, restrict=NULL, cos.norm=TRUE, ndist=3, d=50, auto.order=FALSE, compute.variances=FALSE, 
        subset.row=NULL, correct.all=FALSE, pc.input=FALSE, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL, 
        BSPARAM=ExactParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    originals <- batches <- list(...)
    
    # Pulling out information from the SCE objects.        
    if (checkIfSCE(batches)) {
        checkBatchConsistency(batches)
        checkSpikeConsistency(batches)
        restrict <- checkRestrictions(batches, restrict)

        pc.input <- !is.null(use.dimred)
        if (pc.input) {
            batches <- lapply(batches, reducedDim, type=use.dimred, withDimnames=FALSE)
            checkBatchConsistency(batches, cells.in.columns=FALSE) # re-checking dimensions of reddim matrices.
        } else {
            subset.row <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes)
            batches <- lapply(batches, assay, i=assay.type, withDimnames=FALSE)
        }
    } else {
        checkBatchConsistency(batches, cells.in.columns=!pc.input)
        restrict <- checkRestrictions(batches, restrict, cells.in.columns=!pc.input)
    }

    # Subsetting by 'batch'.
    common.args <-list(k=k, cos.norm=cos.norm, ndist=ndist, d=d, subset.row=subset.row, 
        correct.all=correct.all, auto.order=auto.order, pc.input=pc.input, 
        compute.variances=compute.variances, BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    if (length(batches)==1L) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        output <- do.call(.fast_mnn_single, c(list(x=batches[[1]], batch=batch, restrict=restrict[[1]]), common.args))
    } else {
        output <- do.call(.fast_mnn_list, c(list(batch.list=batches, restrict=restrict), common.args))
    }

    if (pc.input) {
        output$corrected <- .rename_output(output$corrected, originals, cells.in.columns=FALSE)
        rownames(output) <- rownames(output$corrected)
    } else {
        output <- .rename_output(output, originals, subset.row=subset.row)
    }
    output
}

############################################
# Internal main function, to separate the input handling from the actual calculations.

#' @importFrom BiocSingular ExactParam 
#' @importFrom BiocParallel SerialParam
.fast_mnn_list <- function(batch.list, ..., subset.row=NULL, cos.norm=TRUE, pc.input=FALSE, d=50, correct.all=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
    nbatches <- length(batch.list) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }

    # Creating the PCA input, if it is not already low-dimensional.
    if (!pc.input) {
        if (!is.null(subset.row)) {
            batch.list <- lapply(batch.list, "[", i=subset.row, , drop=FALSE) # Need the extra comma!
        }
        if (cos.norm) { 
            batch.list <- lapply(batch.list, FUN=cosineNorm, mode="matrix")
        }
        pc.mat <- .multi_pca_list(batch.list, d=d, rotate.all=correct.all, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    } else {
        pc.mat <- batch.list
    }

    out <- .fast_mnn(batches=pc.mat, ..., BPPARAM=BPPARAM)
    if (!pc.input) {
        out <- .convert_to_SCE(out, pc.mat)
    }
    out
}

#' @importFrom BiocSingular ExactParam
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors List metadata metadata<-
.fast_mnn_single <- function(x, batch, restrict=NULL, ..., subset.row=NULL, cos.norm=TRUE, pc.input=FALSE, d=50, correct.all=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
    batch <- factor(batch)
    if (nlevels(batch) < 2L) { 
        stop("at least two batches must be specified") 
    }

    # Creating the PCA input, if it is not already low-dimensional.
    if (!pc.input) {
        if (!is.null(subset.row)) {
            x <- x[subset.row,,drop=FALSE]
        }
        if (cos.norm) { 
            x <- cosineNorm(x, mode="matrix")
        }
        mat <- .multi_pca_single(x, batch=batch, d=d, rotate.all=correct.all, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    } else {
        mat <- List(x)
    }

    divided <- divideIntoBatches(mat[[1]], batch=batch, restrict=restrict, byrow=TRUE)
    output <- .fast_mnn(batches=divided$batches, restrict=divided$restricted, ..., BPPARAM=BPPARAM)

    # Reordering by the input order.        
    d.reo <- divided$reorder
    output <- output[d.reo,,drop=FALSE]

    rev.order <- integer(length(d.reo))
    rev.order[d.reo] <- seq_along(d.reo)
    pairings <- metadata(output)$pairs
    for (x in seq_along(pairings)) {
        pairings[[x]][,1] <- rev.order[pairings[[x]][,1]]
        pairings[[x]][,2] <- rev.order[pairings[[x]][,2]]
    }
    metadata(output)$pairs <- pairings

    if (!pc.input) {
        output <- .convert_to_SCE(output, mat)
    }
    output
}

#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame metadata<- 
#' @importFrom BiocNeighbors KmknnParam
.fast_mnn <- function(batches, k=20, restrict=NULL, ndist=3, auto.order=FALSE, compute.variances=FALSE, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    nbatches <- length(batches)
    var.kept <- rep(1, nbatches)

    # Defining ordering vectors.
    if (!is.logical(auto.order)) {
        re.order <- as.integer(auto.order)
        if (!identical(sort(re.order), seq_len(nbatches))) {
            stop("integer 'auto.order' must contain a permutation of 1:nbatches") 
        }
        mnn.store <- MNN_supplied_order(batches, restrict, re.order)
    } else if (auto.order) {
        mnn.store <- MNN_auto_order(batches, restrict)
    } else {
        mnn.store <- MNN_supplied_order(batches, restrict)
    }

    mnn.pairings <- vector("list", nbatches-1L)
    for (bdx in 2:nbatches) {
        mnn.store <- .advance(mnn.store, k=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        refdata <- .get_reference(mnn.store)
        curdata <- .get_current(mnn.store)

        # Estimate the overall batch vector (implicitly 'restrict'd by definition of MNNs).
        mnn.sets <- .get_mnn_result(mnn.store) 
        ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        # Remove variation along the batch vector, which responds to 'restrict'.
        # Also recording the lost variation if desired, which does not respond to 'restrict'.
        if (compute.variances) {
            var.before <- .compute_intra_var(refdata, curdata, batches, .get_reference_indices(mnn.store))
        }

        refdata <- .center_along_batch_vector(refdata, overall.batch, restrict=.get_reference_restrict(mnn.store))
        curdata <- .center_along_batch_vector(curdata, overall.batch, restrict=.get_current_restrict(mnn.store))

        if (compute.variances) {
            var.after <- .compute_intra_var(refdata, curdata, batches, .get_reference_indices(mnn.store))
            var.kept[seq_len(bdx)] <- var.kept[seq_len(bdx)] * var.after/var.before
        }

        # Recompute correction vectors and apply them to all cells (hence, no restriction).
        re.ave.out <- .average_correction(refdata, mnn.sets$first, curdata, mnn.sets$second)
        curdata <- .tricube_weighted_correction(curdata, re.ave.out$averaged, re.ave.out$second, k=k, ndist=ndist, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        mnn.pairings[[bdx-1L]] <- DataFrame(first=mnn.sets$first, second=mnn.sets$second + nrow(refdata))
        mnn.store <- .compile(mnn.store, new.reference=refdata, curdata)
    }

    refdata <- .get_reference(mnn.store)
    merge.order <- .get_reference_indices(mnn.store)

    # Reporting the batch identities.
    ncells.per.batch <- vapply(batches, FUN=nrow, FUN.VALUE=0L)
    batch.names <- .create_batch_names(names(batches), ncells.per.batch)

    # Adjusting the output back to the input order in 'batches'.
    if (is.unsorted(merge.order)) {
        ordering <- .restore_original_order(merge.order, ncells.per.batch)
        refdata <- refdata[ordering,,drop=FALSE]
        mnn.pairings <- .reindex_pairings(mnn.pairings, ordering)
        var.kept[merge.order] <- var.kept
    }
    
    # Formatting the output.
    output <- DataFrame(corrected=I(refdata), batch=batch.names$ids)
    m.out <- list(pairs=mnn.pairings, order=batch.names$labels[merge.order])
    if (compute.variances) {
        m.out$lost.var <- 1 - var.kept
    }
    metadata(output) <- m.out
    return(output)
}

############################################
# Correction-related functions.

.average_correction <- function(refdata, mnn1, curdata, mnn2)
# Computes correction vectors for each MNN pair, and then
# averages them for each MNN-involved cell in the second batch.
{
    corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
    corvec <- rowsum(corvec, mnn2)
    npairs <- table(mnn2)

    second.names <- as.character(names(npairs))
    stopifnot(identical(second.names, rownames(corvec)))

    corvec <- unname(corvec)/as.vector(npairs)
    list(averaged=corvec, second=as.integer(second.names))
}

.center_along_batch_vector <- function(mat, batch.vec, restrict=NULL) 
# Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
# This removes any variation along the overall batch vector within each matrix.
{
    batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
    batch.loc <- as.vector(mat %*% batch.vec)

    if (is.null(restrict)) {
        central.loc <- mean(batch.loc) 
    } else { 
        central.loc <- mean(batch.loc[restrict])
    }

    mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
    return(mat)
}

#' @importFrom BiocNeighbors queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, BNPARAM=KmknnParam(), BPPARAM=SerialParam())
# Computing tricube-weighted correction vectors for individual cells,
# using the nearest neighbouring cells _involved in MNN pairs_.
{
    cur.uniq <- curdata[in.mnn,,drop=FALSE]
    safe.k <- min(k, nrow(cur.uniq))
    closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
    weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
    curdata + weighted.correction
}

############################################
# Variance calculation functions.

#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colVars
.compute_intra_var <- function(reference, current, all.pcs, order) {
    all.var <- numeric(length(order))

    last <- 0L
    dm <- DelayedArray(reference)
    for (i in seq_len(length(order)-1L)) {
        cur.ncells <- nrow(all.pcs[[order[i]]])
        chosen <- last + seq_len(cur.ncells)
        all.var[i] <- sum(colVars(dm, rows=chosen))
        last <- last + cur.ncells
    }

    all.var[length(order)] <- sum(colVars(DelayedArray(current)))
    return(all.var)
}

############################################
# Output formatting functions.

#' @importFrom S4Vectors metadata DataFrame
#' @importFrom BiocSingular LowRankMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment
.convert_to_SCE <- function(corrected.df, pc.mat) {
    rot <- metadata(pc.mat)$rotation
    mat <- LowRankMatrix(rot, corrected.df$corrected)
    SingleCellExperiment(list(reconstructed=mat),
        colData=DataFrame(batch=corrected.df$batch),
        rowData=DataFrame(rotation=I(rot)),
        metadata=metadata(corrected.df),
        reducedDims=list(corrected=corrected.df$corrected))
}
