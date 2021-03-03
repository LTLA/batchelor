#' Computes the variance of the paired MNN deltas
#'
#' @inheritParams fastMNN
#' @param pairs A \linkS4class{DataFrame} or list of \linkS4class{DataFrame}s containing MNN pairing information.
#' Each row of each DataFrame specifies an MNN pair; each DataFrame should have two columns containing the column indices of paired cells.
#' This is typically produced by \code{\link{fastMNN}}, see that documentation for more information.
#' @param trend.args Named list of further arguments to pass to \code{fitTrendVar} from the \pkg{scran} package.
#' @param subset.row A vector specifying which genes to use for the calculation.
#' This is typically the same as the \code{subset.row} passed to \code{\link{fastMNN}}, if any.
#' @param cos.norm A logical scalar indicating whether cosine normalization should be performed on the expression values.
#' @param compute.all Logical scalar indicating whether statistics should be returned for all genes when \code{subset.row} is specified, see DEtails.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether the calculations should be parallelized.
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
#' To eliminate the mean-variance relationship, we fit a trend to the variances across all genes with respect to the mean log-expression in each MNN pair.
#' We then define an adjusted variance based on the residuals of the fitted curve; 
#' this is a more appropriate measure to use for ranking affected genes, as it ensures that genes are not highly ranked due to their abundance alone.
#' Fitting is done using the \code{fitTrendVar} function from the \pkg{scran} package - further arguments can be passed via \code{trend.args}.
#' 
#' If a list is passed to \code{pairs}, each entry of the list is assumed to correspond to a merge step.
#' The variance and trend fitting is done separately for each merge step, and the results are combined by computing the average variance across all steps.
#' This avoids considering the variance of the deltas across merge steps, which is not particularly concerning, e.g., if different batches require different translations.
#'
#' If \code{cos.norm=TRUE}, cosine normalization is performed on each cell in \code{...} with \code{\link{cosineNorm}}.
#' This mimics what is done inside \code{\link{fastMNN}} for greater consistency with what the batch correction operates on.
#' While the scale of the variances are no longer comparable to that of log-expression data, this is not a major concern as we are mostly interested in the rankings anyway.
#'
#' If \code{subset.row} is specified, variances are only computed for the requested subset of genes, most typically the set of highly variable genes used in \code{\link{fastMNN}}.
#' This also implies that the normalization and trend fitting is limited to the specified subset.
#' However, if \code{compute.all=TRUE}, the scaling factor and fitted trend are extrapolated to compute adjusted variances for all other genes.
#' This is useful for picking up genes outside of the subset used in the correction.
#'
#' @return
#' A \linkS4class{DataFrame} with one row per gene in \code{...} (or as specified by \code{subset.row}),
#' containing the following columns:
#' \itemize{
#' \item \code{mean}, the mean of the mean log-expressions across all MNN pairs.
#' This may contain repeated contributions from the same cell if it is involved in many MNN pairs.
#' \item \code{total}, the total variance of the deltas across all MNN pairs.
#' \item \code{trend}, the fitted values of the trend in \code{total} with respect to \code{mean}.
#' \item \code{adjusted}, the adjusted variance, i.e., the residuals of the fitted trend.
#' }
#' The \code{\link{metadata}} contains the trend fitting statistics returned by \code{fitTrendVar}.
#'
#' If \code{pairs} is a list of length greater than 1, the returned DataFrame will also contain \code{per.block}, a nested DataFrames of nested DataFrames.
#' Each one of these contains statistics for the individual merge steps and has the same structure as that described above.
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
#' @seealso
#' \code{\link{fastMNN}} and related functions, to obtain the MNN pairs in the first place.
#' 
#' @export
#' @importFrom DelayedArray DelayedArray blockApply rowAutoGrid
#' @importFrom scuttle .bpNotSharedOrUp .unpackLists .subset2index
#' @importFrom BiocParallel bpstart bpstop SerialParam
mnnDeltaVariance <- function(..., pairs, cos.norm=TRUE, subset.row=NULL, compute.all=FALSE, assay.type="logcounts", BPPARAM=SerialParam(), trend.args=list()) {
    batches <- .unpackLists(...)
    checkBatchConsistency(batches, cells.in.columns=TRUE)

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    # Extracting information from SCEs.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }

    # Creating a common matrix.
    x <- lapply(batches, DelayedArray)

    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x[[1]], byrow=TRUE)
        if (!compute.all) {
            x <- lapply(x, function(y) y[subset.row,,drop=FALSE])
            subset.row <- NULL
        }
    }

    if (cos.norm) {
        l2 <- lapply(x, cosineNorm, BPPARAM=BPPARAM, subset.row=subset.row, mode="l2norm")
        x <- mapply(FUN=.apply_cosine_norm, x=x, l2=l2)
    }

    x <- do.call(cbind, x)

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

        xargs <- list(xmean, xvar)
        if (!is.null(subset.row)) {
            xargs <- lapply(xargs, function(X) X[subset.row])
        }
        fit <- do.call(scran::fitTrendVar, c(xargs, trend.args))

        # Just putting 'p.value' and 'FDR' as placeholders for combineBlocks.
        df <- DataFrame(mean = xmean, total = xvar, trend = fit$trend(xmean), p.value=1, FDR=1)
        df$adjusted <- df$total - df$trend
        metadata(df) <- fit
        output[[i]] <- df
    }

    npairs <- vapply(pairs, nrow, 0L)
    combined <- scran::combineBlocks(output, 
        ave.fields=c("mean", "total", "trend", "adjusted"), 
        pval.field="p.value", method="fisher", 
        geometric=FALSE, equiweight=FALSE, 
        weights=npairs, valid=npairs >= 2L)

    combined$p.value <- NULL
    combined$FDR <- NULL
    colnames(combined)[colnames(combined) %in% "per.block"] <- "per.step"
    for (i in seq_along(combined$per.step)) {
        combined$per.step[[i]]$p.value <- NULL
        combined$per.step[[i]]$FDR <- NULL
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
