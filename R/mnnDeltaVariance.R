#' Computes the variance of the paired MNN deltas
#'
#' @inheritParams fastMNN
#' @param pairs A \linkS4class{DataFrame} or list of \linkS4class{DataFrame}s containing MNN pairing information.
#' Each row of each DataFrame specifies an MNN pair; each DataFrame should have two columns containing the column indices of paired cells.
#' This is typically produced by \code{\link{fastMNN}}, see that documentation for more information.
#' @param trend.args Named list of further arguments to pass to \code{fitTrendVar} from the \pkg{scran} package.
#' 
#' @details
#' The \dQuote{MNN delta} is defined as the difference in the \emph{uncorrected} expression values of MNN-paired cells.
#' We compute the variance of these values across all pairs for each gene; genes with highly variable deltas are strongly affected by the correction beyond a simple shift.
#' By looking through the top genes, we may conclude that the correction did something unusual if, e.g., we find important markers in the top set.
#' Conversely, if key genes do not have unusually high variances in their delta, we may gain some confidence in the reliability of the correction -
#' at least with respect to the behavior of those genes.
#'
#' We compute the variance of the deltas rather than using their magnitude, 
#' as a shift in expression space is not particularly concerning (and is, in fact, par for the course for batch correction).
#' Rather, we are interested in identifying non-linear changes that might be indicative of correction errors, e.g., inappropriate merging of two different subpopulations.
#' This would manifest as high variances for the relevant marker genes.
#' Of course, whether or not this is actually an error can be a subjective decision - loss of some biological heterogeneity is often an acceptable cost for flexible correction.
#'
#' @return
#' A numeric vector of variances for each row in the input objects (or as specified by \code{subset.row}).
#'
#' @author Aaron Lun
#'
#' @examples
#' means <- 2^rexp(200) * 10
#' B1 <- matrix(rpois(10000, means), ncol=50) # Batch 1 
#' B2 <- matrix(rpois(10000, means), ncol=50) # Batch 2
#'
#' # Spiking in some population structure. Each batch is split into two halves 
#' # corresponding to 2 subpopulations. These populations match up across batches
#' # but the group in batch 1 also upregulates the first gene.
#' B1[1:2,1:25] <- B1[1:2,1:25] * 50 
#' B2[2,1:25] <- B2[2,1:25] * 50 
#'
#' # Quick log-transformation and correction. We'll omit the multiBatchNorm as the
#' # sequencing depth is the same across batches anyway.
#' library(scuttle)
#' B1 <- normalizeCounts(B1)
#' B2 <- normalizeCounts(B2)
#' out <- fastMNN(B1, B2)
#'
#' # First gene should have high variance of deltas, as the upregulation in B1 is 
#' # removed to enable the alignment of subpopulations across batches.
#' mnnDeltaVariance(B1, B2, pairs=metadata(out)$merge.info$pairs[[1]])
#' 
#' @export
#' @importFrom DelayedArray DelayedArray blockApply rowAutoGrid
mnnDeltaVariance <- function(..., pairs, subset.row=NULL, BPPARAM=SerialParam(), trend.args=list()) {
    batches <- .unpackLists(...)
    checkBatchConsistency(batches, cells.in.columns=TRUE)

    # Extracting information from SCEs.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }

    # Creating a common matrix.
    x <- lapply(batches, DelayedArray)
    x <- do.call(cbind, x)

    if (!is.null(subset.row)) {
        x <- x[subset.row,,drop=FALSE]
    }

    # Subsetting to all cells involved in pairs, to avoid unnecessary extraction.
    if (is(pairs, "DataFrame")) {
        pairs <- list(pairs)
    }

    universe <- lapply(pairs, function(x) union(x[,1], x[,2]))
    universe <- sort(Reduce(union, universe))
    x <- x[,universe,drop=FALSE]

    for (p in seq_along(pairs)) {
        pairs[[p]][,1] <- match(pairs[[p]][,1], universe)
        pairs[[p]][,2] <- match(pairs[[p]][,2], universe)
    }

    # Computing the variance of the deltas at each merge step.
    stats <- blockApply(x, FUN=.compute_mnn_variance, pairs=pairs, grid=rowAutoGrid(x), BPPARAM=BPPARAM)

    output <- vector("list", length(pairs))
    for (i in seq_along(pairs)) {
        xmean <- unlist(lapply(stats, FUN=function(s) s$mean[[i]]), use.names=FALSE)
        xmean <- unname(xmean)
        xvar <- unlist(lapply(stats, FUN=function(s) s$var[[i]]), use.names=FALSE)
        xvar <- unname(xvar)

        # Just putting 'p.value' and 'FDR' as placeholders for combineBlocks.
        fit <- do.call(scran::fitTrendVar, c(list(xmean, xvar), trend.args))
        df <- DataFrame(mean = xmean, total = xvar, tech = fit$trend(xmean), p.value=1, FDR=1)
        df$bio <- df$total - df$tech
        output[[i]] <- df
    }

    npairs <- vapply(pairs, nrow, 0L)
    combined <- scran::combineBlocks(output, 
        ave.fields=c("mean", "total", "tech", "bio"), 
        pval.field="p.value", method="fisher", 
        geometric=FALSE, equiweight=FALSE, 
        weights=npairs, valid=npairs >= 2L)

    combined$p.value <- NULL
    combined$FDR <- NULL
    for (i in seq_along(combined$per.block)) {
        combined$per.block[[i]]$p.value <- NULL
        combined$per.block[[i]]$FDR <- NULL
    }

    combined
}

#' @importFrom methods is as 
#' @importFrom DelayedMatrixStats rowVars
#' @importClassesFrom DelayedArray SparseArraySeed
#' @importClassesFrom Matrix sparseMatrix
#' @importFrom Matrix rowMeans
.compute_mnn_variance <- function(block, pairs) {
    if (is(block, "SparseArraySeed")) {
        block <- as(block, "sparseMatrix")
    }

    all.vars <- all.means <- vector("list", length(pairs))
    for (i in seq_along(pairs)) {
        cpairs <- pairs[[i]]
        b1 <- block[,cpairs[,1],drop=FALSE]
        b2 <- block[,cpairs[,2],drop=FALSE]
        all.vars[[i]] <- rowVars(b1 - b2)
        all.means[[i]] <- (rowMeans(b1) + rowMeans(b2))/2
    }

    list(var=all.vars, mean=all.means)
}
