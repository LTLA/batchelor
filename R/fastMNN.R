#' Fast mutual nearest neighbors correction
#'
#' Correct for batch effects in single-cell expression data using a fast version of the mutual nearest neighbors (MNN) method.
#'
#' @param ... One or more log-expression matrices where genes correspond to rows and cells correspond to columns.
#' Each matrix should contain the same number of rows, corresponding to the same genes in the same order.
#' 
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied containing a log-expression matrix in the \code{assay.type} assay.
#' Note the same restrictions described above for gene expression matrix inputs.
#'
#' If multiple objects are supplied, each object is assumed to contain all and only cells from a single batch.
#' Objects of different types can be mixed together.
#' If a single object is supplied, \code{batch} should also be specified.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param restrict A list of length equal to the number of objects in \code{...}.
#' Each entry of the list corresponds to one batch and specifies the cells to use when computing the correction.
#' @param k An integer scalar specifying the number of nearest neighbors to consider when identifying MNNs.
#' @param prop.k A numeric scalar in (0, 1) specifying the proportion of cells in each dataset to use for mutual nearest neighbor searching.
#' If set, \code{k} for the search in each batch is redefined as \code{max(k, prop.k*N)} where \code{N} is the number of cells in that batch.
#' @param cos.norm A logical scalar indicating whether cosine normalization should be performed on the input data prior to PCA.
#' @param ndist A numeric scalar specifying the threshold beyond which neighbours are to be ignored when computing correction vectors.
#' Each threshold is defined as a multiple of the number of median distances.
#' @param d Number of dimensions to use for dimensionality reduction in \code{\link{multiBatchPCA}}.
#' @param weights Numeric scalar of weights to use in \code{\link{multiBatchPCA}}.
#' @param merge.order An integer vector containing the linear merge order of batches in \code{...}.
#' Alternatively, a list of lists representing a tree structure specifying a hierarchical merge order.
#' @param auto.merge Logical scalar indicating whether to automatically identify the \dQuote{best} merge order.
#' @param auto.order Deprecated, use \code{merge.order} or \code{auto.merge} instead.
#' @param min.batch.skip Numeric scalar specifying the minimum relative magnitude of the batch effect, 
#' below which no correction will be performed at a given merge step.
#' @param subset.row A vector specifying which features to use for correction. 
#' @param correct.all Logical scalar indicating whether a rotation matrix should be computed for genes not in \code{subset.row}.
#' @param pc.input Deprecated, use \code{\link{reducedMNN}} instead.
#' @param assay.type A string or integer scalar specifying the assay containing the log-expression values.
#' Only used for SingleCellExperiment inputs. 
#' @param get.spikes Deprecated.
#' @param use.dimred Deprecated, use \code{\link{reducedMNN}} instead.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.
#' This uses a fast approximate algorithm from \pkg{irlba} by default, see \code{\link{multiBatchPCA}} for details.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the nearest neighbor algorithm.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the PCA and nearest-neighbor searches should be parallelized.
#' 
#' @return
#' A \linkS4class{SingleCellExperiment} is returned where each row is a gene and each column is a cell. 
#' This contains:
#' \itemize{
#' \item A \code{corrected} matrix in the \code{reducedDims} slot, containing corrected low-dimensional coordinates for each cell.
#' This has number of columns equal to \code{d} and number of rows equal to the total number of cells in \code{...}.
#' \item A \code{batch} column in the \code{colData} slot, containing the batch of origin for each row (i.e., cell) in \code{corrected}.
#' \item A \code{rotation} column the \code{rowData} slot, containing the rotation matrix used for the PCA.
#' This has \code{d} columns and number of rows equal to the number of genes to report (see the \dQuote{Choice of genes} section).
#' \item A \code{reconstructed} matrix in the \code{assays} slot, containing the low-rank reconstruction of the expression matrix.
#' This can be interpreted as per-gene corrected log-expression values (after cosine normalization, if \code{cos.norm=TRUE}) but should not be used for quantitative analyses.
#' This has number of rows equal to the number of input genes if \code{subset.row=NULL} or \code{correct.all=TRUE},
#' otherwise each row corresponds to a gene in \code{subset.row}.
#' }
#' 
#' Cells in the output object are always ordered in the same manner as supplied in \code{...}.
#' For a single input object, cells will be reported in the same order as they are arranged in that object.
#' In cases with multiple input objects, the cell identities are simply concatenated from successive objects,
#' i.e., all cells from the first object (in their provided order), then all cells from the second object, and so on.
#' This is true regardless of the value of \code{merge.order} and \code{auto.merge},
#' which only affects the internal merge order of the batches.
#' 
#' The metadata of the output object contains \code{merge.info}, 
#' a \linkS4class{DataFrame} of diagnostic information about each merge step.
#' See the \dQuote{Merge diagnostics} section for more details.
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
#' The \code{batch} argument allows users to easily perform batch correction when all cells have already been combined into a single object.
#' This avoids the need to manually split the matrix or SingleCellExperiment object into separate objects for input into \code{fastMNN}.
#' In this situation, the order of input batches is defined by the order of levels in \code{batch}.
#'
#' @section Controlling the merge order:
#' By default, batches are merged in the user-supplied order, 
#' i.e., the first batch is merged with the second batch, 
#' the third batch is merged with the combined first-second batch,
#' the fourth batch is merged with the combined first-second-third batch and so on.
#' We refer to this approach as a progressive merge.
#'
#' If \code{merge.order} is an integer vector, it is treated as an ordering permutation with which to perform a progressive merge.
#' For example, if \code{merge.order=c(4,1,3,2)}, batches 4 and 1 in \code{...} are merged first;
#' batch 3 is merged with the combined 4+1 batch; 
#' and then batch 2 is merged with the combined 4+1+3 batch.
#' This is often more convenient than changing the order manually in \code{...}, 
#' which would alter the order of batches in the output \code{corrected} matrix.
#'
#' If \code{merge.order} is a character vector, it is treated as an ordering permutation for named batches.
#'
#' If \code{merge.order} is a nested list, it is treated as a tree that specifies a hierarchical merge.
#' Each element of the list should either be a string or integer scalar, corresponding to a leaf node that specifies a batch;
#' or another list, corresponding to an internal node that should contain at least two children;
#' or an integer or character vector of length 2 or more, again corresponding to an internal node.
#' \itemize{
#' \item For example, \code{list(list(1,2), list(3,4))} indicates that batch 1 should be merged with batch 2;
#' batch 3 should be merged with batch 4; and that, finally, the combined batches 1+2 and 3+4 should be merged.
#' \item More than two children per node are supported and will result in a progressive merge within that node.
#' For example, \code{list(list(1,2,3), list(4,5,6))} will merge batch 1 with 2, then 1+2 with 3;
#' batch 4 with 5, and then 4+5 with 6; and finally, 1+2+3 with 4+5+6.
#' \item The same approach can be used for integer or character vectors, e.g., \code{list(1:3, 4:6)} has the same effect as above.
#' }
#' Note that, while batches can be specified by name (character) or index (integer), users cannot use both in the same tree.
#'
#' The merge order may occasionally be important as it determines the number of MNN pairs available at each merge step.
#' MNN pairs results in greater stability of the batch vectors and increased likelihood of identifying shared subpopulations,
#' which are important to the precision and accuracy of the MNN-based correction, respectively.
#' \itemize{
#' \item  In a progressive merge, the reference increases in size at each step,
#' ensuring that more cells are available to identify MNN pairs in later merges.
#' We suggest setting the largest, most heterogeneous batch as the first reference,
#' which favors detection of sufficient MNN pairs between the first and other batches.
#' Conversely, if two small batches without shared populations are supplied first, 
#' the wrong MNN pairs will be detected and the result of the merge will be incorrect.
#' \item A merge tree is useful for merging together batches that are known to be more closely related (e.g., replicates)
#' before attempting difficult merges involving more dissimilar batches.
#' The idea is to increase the number of cells and thus MNN pairs prior to merging batches with few shared subpopulations.
#' By comparison, performing the more difficult merges first is more likely to introduce errors whereby distinct subpopulations are incorrectly placed together, which is propagated to later steps as the initial merge is used as a reference for subsequent merges.
#' \item If \code{auto.merge=TRUE}, merge steps are chosen to maximize the number of MNN pairs at each step.
#' The aim is to improve the stability of the correction by first merging more similar batches with more MNN pairs.
#' This can be somewhat time-consuming as MNN pairs need to be iteratively recomputed for all possible batch pairings.
#' }
#'
#' The order of cells in the output is \emph{never} affected by the setting of \code{merge.order} or \code{auto.order}.
#' It depends only on the order of objects in \code{...} and the order of cells within each object.
#'
#' @section Choice of genes:
#' All genes are used with the default setting of \code{subset.row=NULL}.
#' Users can set \code{subset.row} to subset the inputs to highly variable genes or marker genes.
#' This improves the quality of the PCA and identification of MNN pairs by reducing the noise from irrelevant genes.
#' Note that users should not be too restrictive with subsetting, 
#' as high dimensionality is required to satisfy the orthogonality assumption in MNN detection.
#'
#' By default, only the selected genes are used to compute rotation vectors and a low-rank corrected expression matrix.
#' However, setting \code{correct.all=TRUE} will return rotation vectors that span all genes in the supplied input data.
#' This is useful for ensuring that corrected values are returned for all input genes, e.g., in \code{\link{correctExperiments}}.
#' Note that this setting will not affect the corrected low-dimension coordinates or the rotation values for the selected genes. 
#'
#' @section Using restriction:
#' See \code{?"\link{batchelor-restrict}"} for a description of the \code{restrict} argument.
#' Specifically, \code{fastMNN} will only use the restricted subset of cells in each batch to identify MNN pairs and the center of the orthogonalization.
#' It will then extrapolate the correction to all cells in each batch.
#'
#' Note that \emph{all} cells are used to perform the PCA, regardless of whether \code{restrict} is set.
#' This is generally desirable in applications where \code{restrict} is useful.
#' For example, constructing the projection vectors with only control cells will not guarantee resolution of unique non-control populations in each batch.
#' 
#' However, this also means that the corrected values for the restricted cells will differ from the output when the inputs are directly subsetted to only contain the restricted cells.
#' If this is not desirable, users can perform the PCA manually and apply \code{\link{reducedMNN}} instead.
#'
#' @section Merge diagnostics:
#' We can consider \code{fastMNN}'s operation in terms of pairwise merge steps.
#' Each merge step involves two mutually exclusive sets of cells, a \dQuote{left} set and \dQuote{right} set.
#' Each set may consist of cells from different batches if those batches were merged in a previous step.
#' The merge will then create a new set of cells that combines the left and right sets.
#' Iteratively repeating this process with the newly formed sets will eventually merge all batches together.
#'
#' The output metadata contains \code{merge.info}, a DataFrame where each row corresponds to a merge step.
#' It contains the following fields:
#' \itemize{
#' \item \code{left}, a \linkS4class{List} of integer or character vectors.
#' Each vector specifies the batches in the left set at a given merge step. 
#' \item \code{right}, a similar List of integer or character vectors.
#' Each vector specifies the batches in the right set at a given merge step. 
#' \item \code{pairs}, a List of DataFrames specifying which pairs of cells were identified as MNNs at each step.
#' In each DataFrame, each row corresponds to a single MNN pair and specifies the
#' paired cells that were in the left and right sets, respectively.
#' Note that the indices refer to those paired cells in the \emph{output} ordering of cells,
#' i.e., users can identify the paired cells at each step by column-indexing the output of the \code{fastMNN} function.
#' \item \code{batch.size}, a numeric vector specifying the relative magnitude of the batch effect at each merge,
#' see \dQuote{Orthogonalization details}.
#' \item \code{skipped}, a logical vector indicating whether the correction was skipped 
#' if the magnitude of the batch effect was below \code{min.batch.skip}.
#' \item \code{lost.var}, a numeric matrix specifying the percentage of variance lost due to orthogonalization at each merge step.
#' This is reported separately for each batch (columns, ordered according to the input order, \emph{not} the merge order).
#' }
#'
#' @section Specifying the number of neighbors:
#' The threshold to define nearest neighbors is defined by \code{k}, which is passed to \code{\link{findMutualNN}} to identify MNN pairs.
#' The size of \code{k} can be roughly interpreted as the anticipated minimum size of a shared subpopulation in each batch.
#' If a batch has fewer than \code{k} cells of a shared subpopulation, there is an increased risk that its counterparts in other batches will form incorrect MNN pairs.
#' 
#' From the perspective of the algorithm, larger values allow for more MNN pairs to be obtained, which improves the stability of the correction vectors.
#' Larger values also increase robustness against non-orthogonality, by ignoring a certain level of biological variation when identifying pairs.
#' This can be used to avoid the kissing problem where MNN pairs are only detected on the \dQuote{surface} of the distribution.
#' However, values of \code{k} should not be too large, as this would result in MNN pairs being inappropriately identified between biologically distinct populations.
#'
#' In practice, increasing \code{k} will generally result in more aggressive merging as the algorithm is more generous in matching subpopulations across batches.
#' We suggest starting with the default \code{k} and increasing it if one is confident that the same cell types are not adequately merged across batches.
#' This is better than starting with a large \code{k} as incorrect merging is much harder to diagnose than insufficient merging.
#'
#' An additional consideration is that the effect of any given \code{k} will vary with the number of cells in each batch.
#' With more cells, a larger \code{k} may be preferable to achieve better merging in the presence of non-orthogonality.
#' We can achieve this by setting \code{prop.k}, which allows the choice of \code{k} to adapt to the size of each batch at each merge step.
#' This also handles asymmetry in batch sizes via the \code{k1} and \code{k2} arguments in \code{\link{findMutualNN}}.
#' 
#' @section Orthogonalization details:
#' \code{fastMNN} will compute the percentage of variance that is lost from each batch during orthogonalization at each merge step.
#' This represents the variance in each batch that is parallel to the average correction vectors (and hence removed during orthogonalization) at each merge step.
#' Large proportions suggest that there is biological structure that is parallel to the batch effect, 
#' corresponding to violations of the assumption that the batch effect is orthogonal to the biological subspace.
#' 
#' If \code{fastMNN} is called with DataFrame inputs, each DataFrame is assumed to be the result of a previous \code{fastMNN} call
#' and have a set of vectors used for orthogonalization in the merge steps of that previous call.
#' In the current call, \code{fastMNN} will gather all such batch vectors across all DataFrame inputs.
#' Each batch is then re-orthogonalized with respect to each of these vectors.
#' This ensures that the same variation is removed from each batch prior to merging.
#' The variance lost due to this pre-correction orthogonalization is reported in the \code{pre.orthog} field in the output metadata.
#'
#' Orthogonalization may cause problems if there is actually no batch effect, resulting in large losses of variance.
#' To avoid this, \code{fastMNN} will not perform any correction if the relative magnitude of the batch effect is less than \code{min.batch.skip}.
#' The relative magnitude is defined as the L2 norm of the average correction vector divided by the root-mean-square of the L2 norms of the per-MNN pair correction vectors.
#' This will be large when the per-pair vectors are all pointing in the same direction, 
#' and small when the per-pair vectors point in random directions due to the absence of a consistent batch effect.
#' If a large loss of variance is observed along with a small batch effect in a given merge step, users can set \code{min.batch.skip} to simply skip correction in that step.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{cosineNorm}} and \code{\link{multiBatchPCA}}, to obtain the values to be corrected.
#'
#' \code{\link{reducedMNN}}, for a version of the function that operates in low-dimensional space.
#'
#' \code{\link{mnnCorrect}} for the \dQuote{classic} version of the MNN correction algorithm.
#'
#' @examples
#' B1 <- matrix(rnorm(10000, -1), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000, 1), ncol=50) # Batch 2
#' out <- fastMNN(B1, B2)
#'
#' # Corrected values for use in clustering, etc.
#' str(reducedDim(out)) 
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
#' \url{https://MarioniLab.github.io/FurtherMNN2018/theory/description.html}
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom BiocSingular IrlbaParam
#' @importClassesFrom S4Vectors List
#' @importFrom S4Vectors DataFrame metadata<-
#' @importFrom methods as
fastMNN <- function(..., batch=NULL, k=20, prop.k=NULL, restrict=NULL, cos.norm=TRUE, ndist=3, d=50, weights=NULL,
    merge.order=NULL, auto.merge=FALSE, auto.order=NULL, min.batch.skip=0,
    subset.row=NULL, correct.all=FALSE, pc.input=FALSE, assay.type="logcounts", get.spikes=FALSE, use.dimred=NULL, 
    BSPARAM=IrlbaParam(deferred=TRUE), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    batches <- list(...)
    is.sce <- checkIfSCE(batches)
    is.df <- vapply(batches, is, class2="DataFrame", FUN.VALUE=TRUE)
    needs.reorth <- any(is.df)

    # Checking whether the PC status is consistent.
    if (any(is.df)) {
        pc.input <- TRUE
        if (any(is.sce) && is.null(use.dimred)) {
            stop("cannot mix low- and high-dimensional inputs")
        }
    } else if (any(is.sce) && !all(is.sce)) {
        sce.pcs <- !is.null(use.dimred)
        if (sce.pcs != pc.input) {
            stop("cannot mix low- and high-dimensional inputs")
        }
    } else if (all(is.sce)) {
        pc.input <- !is.null(use.dimred)
    }
    if (pc.input) {
        .Deprecated(msg="'pc.input=TRUE' and 'use.dimred=TRUE' are deprecated.\nUse 'reducedMNN' instead.")
        return(reducedMNN(..., batch=batch, k=k, prop.k=prop.k, restrict=restrict, ndist=ndist,
            merge.order=merge.order, auto.merge=auto.merge, auto.order=auto.order,
            min.batch.skip=min.batch.skip, BNPARAM=KmknnParam(), BPPARAM=SerialParam()))
    }

    # Extracting information from SCEs.
    if (any(is.sce)) {
        sce.batches <- batches[is.sce]
        checkBatchConsistency(sce.batches)
        checkSpikeConsistency(sce.batches)
        subset.row <- .SCE_subset_genes(subset.row, sce.batches[[1]], get.spikes)

        # NOTE: do NOT set withDimnames=FALSE, this will break the consistency check.
        sce.batches <- lapply(sce.batches, assay, i=assay.type)
        batches[is.sce] <- sce.batches
    }

    # Batch consistency checks.
    checkBatchConsistency(batches, cells.in.columns=TRUE)
    restrict <- checkRestrictions(batches, restrict, cells.in.columns=TRUE)

    # Setting up the parallelization environment.
    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    # Performing the MNN search.
    common.args <-list(k=k, prop.k=prop.k, cos.norm=cos.norm, ndist=ndist, 
        d=d, weights=weights, subset.row=subset.row, correct.all=correct.all, 
        min.batch.skip=min.batch.skip, 
        merge.order=merge.order, auto.merge=auto.merge, auto.order=auto.order,
        BSPARAM=BSPARAM, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    if (length(batches)==1L) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }
        output <- do.call(.fast_mnn_single, c(list(x=batches[[1]], batch=batch, restrict=restrict[[1]]), common.args))
    } else {
        output <- do.call(.fast_mnn_list, c(list(batch.list=batches, restrict=restrict), common.args))
    }

    # Adding names.
    feat.names <- rownames(batches[[1]])
    if (!is.null(subset.row) && !correct.all) {
        feat.names <- feat.names[.row_subset_to_index(batches[[1]], subset.row)]
    }
    rownames(output) <- feat.names
 
    output
}

############################################
# Functions to prepare data (specifically, generate the PCs) from a list or a single batch.

#' @importFrom BiocSingular ExactParam 
#' @importFrom BiocParallel SerialParam
.fast_mnn_list <- function(batch.list, ..., subset.row=NULL, cos.norm=TRUE, d=50, weights=NULL,
    correct.all=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
    nbatches <- length(batch.list) 
    if (nbatches < 2L) { 
        stop("at least two batches must be specified") 
    }

    if (cos.norm) { 
        all.l2s <- lapply(batch.list, FUN=cosineNorm, mode="l2norm", subset.row=subset.row, BPPARAM=BPPARAM)
        batch.list <- mapply(FUN=.apply_cosine_norm, batch.list, all.l2s, SIMPLIFY=FALSE) 
    }

    pc.mat <- .multi_pca_list(batch.list, d=d, weights=weights, subset.row=subset.row,
        get.all.genes=correct.all, BSPARAM=BSPARAM, BPPARAM=BPPARAM)

    out <- .fast_mnn(batches=pc.mat, ..., BPPARAM=BPPARAM)
    .convert_to_SCE(out, pc.mat)
}

#' @importFrom BiocSingular ExactParam
#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors List metadata metadata<-
.fast_mnn_single <- function(x, batch, restrict=NULL, ..., subset.row=NULL, cos.norm=TRUE, 
    d=50, weights=NULL, correct.all=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
    batch <- factor(batch)
    if (nlevels(batch) < 2L) { 
        stop("at least two batches must be specified") 
    }

    if (cos.norm) { 
        l2 <- cosineNorm(x, mode="l2norm", subset.row=subset.row, BPPARAM=BPPARAM)
        x <- .apply_cosine_norm(x, l2)
    }

    mat <- .multi_pca_single(x, batch=batch, d=d, weights=weights, subset.row=subset.row,
        get.all.genes=correct.all, BSPARAM=BSPARAM, BPPARAM=BPPARAM)

    divided <- divideIntoBatches(mat[[1]], batch=batch, restrict=restrict, byrow=TRUE)
    output <- .fast_mnn(batches=divided$batches, restrict=divided$restricted, ..., BPPARAM=BPPARAM)

    # Reordering by the input order.        
    d.reo <- divided$reorder
    output <- output[d.reo,,drop=FALSE]
    metadata(output)$merge.info$pairs <- .reindex_pairings(metadata(output)$merge.info$pairs, d.reo)
    .convert_to_SCE(output, mat)
}

############################################
# Functions to perform the correction, given a set of PC coordinates.
# Split into a wrapper around a core function to distinguish between 
# calculations with a pre-defined vs automatically determined merge tree. 

#' @importFrom BiocParallel SerialParam
#' @importFrom BiocNeighbors KmknnParam
#' @importClassesFrom S4Vectors List
.fast_mnn <- function(batches, k=20, prop.k=NULL, restrict=NULL, ndist=3, 
    merge.order=NULL, auto.merge=FALSE, auto.order=NULL, 
    min.batch.skip=0, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    if (!is.null(auto.order)) {
        .Deprecated(old="auto.order", new="auto.merge")
        if (isTRUE(auto.order)) {
            auto.merge <- TRUE
        } else if (is.numeric(auto.order)) {
            merge.order <- auto.order
        }
    }

    if (!auto.merge) {
        merge.tree <- .create_tree_predefined(batches, restrict, merge.order)
        UPDATE <- .update_tree
        NEXT <- .get_next_merge
    } else {
        mnn.args <- list(k=k, prop.k=prop.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
        merge.tree <- do.call(.initialize_auto_search, c(list(batches, restrict), mnn.args))
        UPDATE <- function(remainders, chosen, ...) {
            .update_remainders(remainders, chosen, ..., mnn.args=mnn.args)
        }
        NEXT <- .pick_best_merge
    }

    output <- .fast_mnn_core(merge.tree, k=k, prop.k=prop.k, restrict=restrict, ndist=ndist, 
        min.batch.skip=min.batch.skip, BNPARAM=BNPARAM, BPPARAM=BPPARAM, 
        NEXT=NEXT, UPDATE=UPDATE)

    nms <- names(batches)
    if (!is.null(nms)) {
        if (anyDuplicated(nms)) {
            stop("names of batches should be unique")
        }
        output$batch <- nms[output$batch]
        output <- .fix_names_in_merge_info(output, nms)
        colnames(metadata(output)$merge.info$lost.var) <- nms
    }
    output
}

#' @importFrom BiocParallel SerialParam
#' @importFrom S4Vectors DataFrame metadata<- 
#' @importFrom BiocNeighbors KmknnParam
#' @importClassesFrom S4Vectors List
#' @importFrom methods as
.fast_mnn_core <- function(merge.tree, k=20, prop.k=NULL, restrict=NULL, ndist=3, 
    min.batch.skip=0, BNPARAM=KmknnParam(), BPPARAM=SerialParam(),
    NEXT, UPDATE)
{
    # Filling in output containers.
    nbatches <- length(unlist(merge.tree))
    nmerges <- nbatches - 1L
    mnn.pairings <- left.set <- right.set <- vector("list", nmerges)
    batch.size <- rep(NA_real_, nmerges)
    skipped <- logical(nmerges)
    var.kept <- matrix(1, nmerges, nbatches) 

    for (mdx in seq_len(nmerges)) {
        # Traversing the merge tree to find the next two batches to merge.
        next.step <- NEXT(merge.tree) 
        left <- next.step$left
        right <- next.step$right

        left.data <- .get_node_data(left)
        left.restrict <- .get_node_restrict(left)
        left.index <- .get_node_index(left)
        left.origin <- .get_node_origin(left)
        left.extras <- .get_node_extras(left)

        right.data <- .get_node_data(right)
        right.restrict <- .get_node_restrict(right)
        right.index <- .get_node_index(right)
        right.origin <- .get_node_origin(right)
        right.extras <- .get_node_extras(right)

        # Computing the variance *before* attempting the merge.
        left.old.var <- .compute_perbatch_var(left.data, left.index, left.origin)
        right.old.var <- .compute_perbatch_var(right.data, right.index, right.origin)
        left.set[[mdx]] <- left.index
        right.set[[mdx]] <- right.index

        # Orthogonalizing and obtaining all MNNs.
        right.data <- .orthogonalize_other(right.data, right.restrict, left.extras)
        left.data <- .orthogonalize_other(left.data, left.restrict, right.extras)

        mnn.sets <- .restricted_mnn(left.data, left.restrict, right.data, right.restrict, 
            k=k, prop.k=prop.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        # Estimate the overall batch vector (implicitly 'restrict'd by definition of MNNs).
        ave.out <- .average_correction(left.data, mnn.sets$first, right.data, mnn.sets$second)
        overall.batch <- colMeans(ave.out$averaged)

        do.correct <- TRUE
        if (!is.na(min.batch.skip)) {
            batch.magnitude <- .get_batch_magnitude(ave.out$averaged, overall.batch)
            batch.size[mdx] <- batch.magnitude

            if (batch.magnitude < min.batch.skip) {
                do.correct <- FALSE
                skipped[mdx] <- TRUE
            }
        }

        if (do.correct) {
            # Remove variation along the batch vector, which responds to 'restrict'.
            left.data <- .center_along_batch_vector(left.data, overall.batch, restrict=left.restrict)
            right.data <- .center_along_batch_vector(right.data, overall.batch, restrict=right.restrict)

            # Also recording the lost variation if desired, which does not respond to 'restrict'.
            left.new.var <- .compute_perbatch_var(left.data, left.index, left.origin)
            right.new.var <- .compute_perbatch_var(right.data, right.index, right.origin)
            to.add <- list(overall.batch)

            # Recompute correction vectors and apply them to all cells (hence, no restriction).
            re.ave.out <- .average_correction(left.data, mnn.sets$first, right.data, mnn.sets$second)
            right.data <- .tricube_weighted_correction(right.data, re.ave.out$averaged, re.ave.out$second,
                k=.choose_k(k, prop.k, nrow(right.data)), ndist=ndist, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

        } else {
            to.add <- list()
            left.new.var <- .compute_perbatch_var(left.data, left.index, left.origin)
            right.new.var <- .compute_perbatch_var(right.data, right.index, right.origin)
        }

        # Recompute correction vectors and apply them to all cells (hence, no restriction).
        var.kept[mdx,left.index] <- left.new.var/left.old.var
        var.kept[mdx,right.index] <- right.new.var/right.old.var
        mnn.pairings[[mdx]] <- DataFrame(left=mnn.sets$first, right=mnn.sets$second)

        merge.tree <- UPDATE(merge.tree, next.step$chosen, 
            data=rbind(left.data, right.data),  
            index=c(left.index, right.index),
            restrict=.combine_restrict(left.data, left.restrict, right.data, right.restrict),
            origin=c(left.origin, right.origin),
            extras=c(left.extras, right.extras, to.add))
    }

    full.data <- .get_node_data(merge.tree)
    full.order <- .get_node_index(merge.tree)
    full.origin <- .get_node_origin(merge.tree)

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
    
    # Formatting the output.
    output <- DataFrame(corrected=I(full.data), batch=full.origin)
    mdf <- DataFrame(
        left=I(as(left.set, "List")),
        right=I(as(right.set, "List")),
        pairs=I(as(mnn.pairings, "List")), 
        batch.size=batch.size,
        skipped=skipped,
        lost.var=I(1 - var.kept)
    )

    metadata(output)$merge.info <- mdf
    output
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

.get_batch_magnitude <- function(correction, ave=colMeans(correction)) 
# Standardizes the magnitude of the average batch vector by the
# magnitude of the per-pair batch vectors making it up. Consistent
# per-pair directions and magnitudes should give a standardized
# magnitude of 1, inconsistencies will result in a magnitude of zero.
{
    ave.l2sq <- sum(colMeans(correction^2))
    if (ave.l2sq == 0) {
        0
    } else {
        l2sq <- sum(ave^2)
        sqrt(l2sq / ave.l2sq)
    }
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
# Orthogonalization-related functions.

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

    mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
}

.orthogonalize_other <- function(data, restrict, vectors) {
    for (vec in vectors) {
        data <- .center_along_batch_vector(data, vec, restrict=restrict)
    }
    data
}

#' @importFrom DelayedMatrixStats colVars
#' @importFrom DelayedArray DelayedArray
.compute_perbatch_var <- function(data, index, origin) {
    ref.var <- numeric(length(index))
    for (i in seq_along(index)) {
        rows <- origin==index[i]
        ref.var[i] <- sum(colVars(DelayedArray(data), rows=rows))
    }
    ref.var
}

############################################
# Output formatting functions.

.combine_restrict <- function(left.data, left.restrict, right.data, right.restrict) {
    if (is.null(left.restrict) && is.null(right.restrict)) {
        NULL
    } else {
        if (is.null(left.restrict)) {
            left.restrict <- seq_len(nrow(left.data))
        }
        if (is.null(right.restrict)) {
            right.restrict <- seq_len(nrow(right.data))
        }
        c(left.restrict, right.restrict + nrow(left.data))
    }
}

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
