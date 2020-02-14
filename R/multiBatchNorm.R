#' Per-batch scaling normalization
#'
#' Perform scaling normalization within each batch to provide comparable results to the lowest-coverage batch.
#' 
#' @param ... One or more \linkS4class{SingleCellExperiment} objects containing counts and size factors.
#' Each object should contain the same number of rows, corresponding to the same genes in the same order.
#'
#' If multiple objects are supplied, each object is assumed to contain all and only cells from a single batch.
#' If a single object is supplied, \code{batch} should also be specified.
#'
#' Alternatively, one or more lists of SingleCellExperiments can be provided;
#' this is flattened as if the objects inside were passed directly to \code{...}.
#' @param batch A factor specifying the batch of origin for all cells when only a single object is supplied in \code{...}.
#' This is ignored if multiple objects are present.
#' @param assay.type A string specifying which assay values contains the counts.
#' @param norm.args A named list of further arguments to pass to \code{\link[scater]{logNormCounts}}.
#' @param min.mean A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
#' @param subset.row A vector specifying which features to use for normalization.
#' @param normalize.all A logical scalar indicating whether normalized values should be returned for all genes.
#' @param preserve.single A logical scalar indicating whether to combine the results into a single matrix if only one object was supplied in \code{...}.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether calculations should be parallelized. 
#' 
#' @details
#' When performing integrative analyses of multiple batches, it is often the case that different batches have large differences in sequencing depth.
#' This function removes systematic differences in coverage across batches to simplify downstream comparisons.
#' It does so by resaling the size factors using median-based normalization on the ratio of the average counts between batches.
#' This is roughly equivalent to the between-cluster normalization described by Lun et al. (2016).
#' 
#' This function will adjust the size factors so that counts in high-coverage batches are scaled \emph{downwards} to match the coverage of the most shallow batch.
#' The \code{\link{logNormCounts}} function will then add the same pseudo-count to all batches before log-transformation.
#' By scaling downwards, we favour stronger squeezing of log-fold changes from the pseudo-count, mitigating any technical differences in variance between batches.
#' Of course, genuine biological differences will also be shrunk, but this is less of an issue for upregulated genes with large counts.
#' 
#' Only genes with library size-adjusted average counts greater than \code{min.mean} will be used for computing the rescaling factors.
#' This improves precision and avoids problems with discreteness.
#' By default, we use \code{min.mean=1}, which is usually satisfactory but may need to be lowered for very sparse datasets. 
#'
#' Users can also set \code{subset.row} to restrict the set of genes used for computing the rescaling factors.
#' By default, normalized values will only be returned for genes specified in the subset.
#' Setting \code{normalize.all=TRUE} will return normalized values for all genes.
#'
#' @section Comparison to other normalization strategies:
#' For comparison, imagine if we ran \code{\link{logNormCounts}} separately in each batch prior to correction.
#' Size factors will be computed within each batch, and batch-specific application in \code{\link{logNormCounts}} will not account for scaling differences between batches.
#' In contrast, \code{multiBatchNorm} will rescale the size factors so that they are comparable across batches.
#' This removes at least one difference between batches to facilitate easier correction.
#'
#' \code{\link{cosineNorm}} performs a similar role of equalizing the scale of expression values across batches. 
#' However, the advantage of \code{\link{multiBatchNorm}} is that its output is more easily interpreted - 
#' the normalized values remain on the log-scale and differences can still be interpreted (roughly) as log-fold changes.
#' The output can then be fed into downstream analysis procedures (e.g., HVG detection) in the same manner as typical log-normalized values from \code{\link{logNormCounts}}.
#' 
#' @return
#' A list of SingleCellExperiment objects with normalized log-expression values in the \code{"logcounts"} assay (depending on values in \code{norm.args}).
#' Each object contains cells from a single batch.
#'
#' If \code{preserve.single=TRUE} and \code{...} contains only one SingleCellExperiment, that object is returned with an additional \code{"logcounts"} assay containing normalized log-expression values.
#' The order of cells is not changed.
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{mnnCorrect}} and \code{\link{fastMNN}}, for methods that can benefit from rescaling.
#'
#' \code{\link[scater]{logNormCounts}} for the calculation of log-transformed normalized expression values.
#' 
#' @examples
#' d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
#' sce1 <- SingleCellExperiment(list(counts=d1))
#' sizeFactors(sce1) <- runif(ncol(d1))
#' 
#' d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
#' sce2 <- SingleCellExperiment(list(counts=d2))
#' sizeFactors(sce2) <- runif(ncol(d2))
#' 
#' out <- multiBatchNorm(sce1, sce2)
#' summary(sizeFactors(out[[1]]))
#' summary(sizeFactors(out[[2]]))
#' 
#' @references
#' Lun ATL (2018).
#' Further MNN algorithm development.
#' \url{https://MarioniLab.github.io/FurtherMNN2018/theory/description.html}
#'
#' @export
#' @importFrom BiocGenerics sizeFactors sizeFactors<- cbind
#' @importFrom scater logNormCounts librarySizeFactors
multiBatchNorm <- function(..., batch=NULL, assay.type="counts", norm.args=list(), 
    min.mean=1, subset.row=NULL, normalize.all=FALSE, preserve.single=TRUE, BPPARAM=SerialParam())
{
    batches <- .unpack_batches(...)
    checkBatchConsistency(batches)
    if (length(batches)==0L) { 
        stop("at least one SingleCellExperiment must be supplied") 
    }

    norm.args <- c(norm.args, list(exprs_values=assay.type, center_size_factors=FALSE))
    needs.subset <- !normalize.all && !is.null(subset.row)

    # Setting up the parallelization environment.
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    if (length(batches)==1L) {
        sce <- batches[[1]]
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }

        # We have to split them up for calcAverage anyway,
        # so there's no point writing special code here (unlike multiBatchPCA).
        by.batch <- split(seq_along(batch), batch)
        batches <- by.batch
        for (i in seq_along(by.batch)) {
            batches[[i]] <- sce[,by.batch[[i]],drop=FALSE]
        }
    } else {
        preserve.single <- FALSE
    }

    batches <- .rescale_size_factors(batches, assay.type=assay.type, 
        subset.row=subset.row, min.mean=min.mean, BPPARAM=BPPARAM)

    if (preserve.single) {
        # Reuse sce here to avoid unnecessary cbind to reform a single object.
        if (needs.subset) {
            sce <- sce[subset.row,]
        }

        all.sf <- unlist(lapply(batches, sizeFactors), use.names=FALSE)
        reorder <- order(unlist(by.batch))
        sizeFactors(sce) <- all.sf[reorder]

        do.call(logNormCounts, c(list(x=sce), norm.args))
    } else {
        if (needs.subset) {
            batches <- lapply(batches, FUN="[", i=subset.row, ) # empty argument is important!
        }

        mapply(FUN=logNormCounts, batches, MoreArgs=norm.args, SIMPLIFY=FALSE)
    }
}

#' @importFrom scater librarySizeFactors
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
.rescale_size_factors <- function(batches, assay.type, subset.row, min.mean, BPPARAM) {
    # Centering the endogenous size factors.
    for (b in seq_along(batches)) {
        current <- batches[[b]]

        cursf <- sizeFactors(current)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current, exprs_values=assay.type, subset_row=subset.row, BPPARAM=BPPARAM)
        }

        centering <- mean(cursf)
        sizeFactors(current) <- cursf / centering
        batches[[b]] <- current
    }

    # Adjusting size factors.
    rescaling <- .compute_batch_rescaling(batches, subset.row=subset.row,
        assay.type=assay.type, min.mean=min.mean, BPPARAM=BPPARAM)

    for (idx in seq_along(batches)) {
        current <- batches[[idx]]
        cursf <- sizeFactors(current)

        # Adjusting the size factors for the new reference.
        cursf <- cursf / rescaling[idx] 
        sizeFactors(current) <- cursf
        batches[[idx]] <- current
    }

    batches
}

#' @importFrom scater calculateAverage 
#' @importFrom stats median
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
.compute_batch_rescaling <- function(batches, subset.row, assay.type, min.mean, BPPARAM) 
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    nbatches <- length(batches)
    collected.ave <- lapply(batches, FUN=calculateAverage, exprs_values=assay.type, 
        subset_row=subset.row, BPPARAM=BPPARAM)

    collected.ratios <- matrix(1, nbatches, nbatches)
    for (first in seq_len(nbatches-1L)) {
        first.ave <- collected.ave[[first]]
        first.sum <- sum(first.ave)

        for (second in first + seq_len(nbatches - first)) {
            second.ave <- collected.ave[[second]]
            second.sum <- sum(second.ave)

            # Mimic calcAverage(cbind(first.ave, second.ave)).
            grand.mean <- (first.ave/first.sum + second.ave/second.sum)/2 * (first.sum + second.sum)/2
            keep <- grand.mean >= min.mean

            # Computing it twice for exactly equal results when order of batches is rearranged.
            kept.f <- first.ave[keep]
            kept.s <- second.ave[keep]
            curratio1 <- median(kept.s/kept.f)
            curratio2 <- median(kept.f/kept.s)

            if (!is.finite(curratio1) || curratio1==0 || !is.finite(curratio2) || curratio2==0) {
                stop("median ratio of averages between batches is not finite")
            }

            collected.ratios[first,second] <- curratio1
            collected.ratios[second,first] <- curratio2
        }
    }

    # Rescaling to the lowest coverage batch.
    smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
    collected.ratios[,smallest]
}
