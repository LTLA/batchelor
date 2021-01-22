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
#' @param norm.args A named list of further arguments to pass to \code{\link[scuttle]{logNormCounts}}.
#' @param min.mean A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
#' @param subset.row A vector specifying which features to use for normalization.
#' @param normalize.all A logical scalar indicating whether normalized values should be returned for all genes.
#' @param preserve.single A logical scalar indicating whether to combine the results into a single matrix if only one object was supplied in \code{...}.
#' @param as.altexp String or integer scalar indicating the alternative Experiment to use in the function (see below for details).
#' All entries of \code{...} must contain the specified entry in their \code{\link{altExps}}.
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
#' @section Handling alternative Experiments:
#' If \code{as.altexp} is specified, this function is applied to the specified alternative Experiment from each entry of \code{...}.
#' This provides a convenient way to normalize multiple feature sets via repeated calls to \code{\link{multiBatchNorm}} with different \code{as.altexp}.
#' The result is equivalent to extracting one alternative Experiment from each input SingleCellExperiment and running \code{multiBatchNorm} on that set.
#' No information is shared between main and alternative Experiments, or between the different alternative Experiments.
#'
#' When \code{as.altexp} is specified, any non-\code{NULL} value of \code{subset.row} applies to the features of the alternative Experiment.
#' Similarly, filtering by \code{min.mean} applies to the specified alternative Experiment.
#' Note that the latter may not always be appropriate depending on the scale of counts in the alternative Experiments;
#' the default is based on typical RNA count data, and a lower threshold may be required for features that are less deeply sequenced.
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
#' \code{\link[scuttle]{logNormCounts}} for the calculation of log-transformed normalized expression values.
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
#' @importFrom scuttle logNormCounts librarySizeFactors .unpackLists
#' @importFrom SingleCellExperiment altExp altExp<-
multiBatchNorm <- function(..., batch=NULL, norm.args=list(), 
    min.mean=1, subset.row=NULL, normalize.all=FALSE, preserve.single=TRUE, 
    assay.type="counts", as.altexp=NULL, BPPARAM=SerialParam())
{
    batches <- .unpackLists(...)
    if (!is.null(as.altexp)) {
        originals <- batches
        batches <- lapply(batches, altExp, e=as.altexp)
    }
    checkBatchConsistency(batches)

    # Handling the batch= and preserve.single= options.
    if (length(batches)==1L) {
        if (is.null(batch)) { 
            stop("'batch' must be specified if '...' has only one object")
        }

        if (!preserve.single) {
            by.batch <- split(seq_along(batch), batch)
            FRAGMENT <- function(target) {
                batches <- by.batch
                for (i in seq_along(by.batch)) {
                    batches[[i]] <- target[,by.batch[[i]],drop=FALSE]
                }
                batches
            }

            batches <- FRAGMENT(batches[[1]])
            if (!is.null(as.altexp)) {
                originals <- FRAGMENT(originals[[1]])
            }
        }
    } else if (length(batches)==0L) { 
        stop("at least one SingleCellExperiment must be supplied") 
    } else {
        preserve.single <- FALSE
    }

    # Setting up the parallelization environment.
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    # Computing the averages and the size factors.
    if (preserve.single) {
        stats <- .compute_batch_statistics_single(x=batches[[1]], batch=batch, 
            assay.type=assay.type, subset.row=subset.row, BPPARAM=BPPARAM)
    } else {
        stats <- .compute_batch_statistics_list(batches, 
            assay.type=assay.type, subset.row=subset.row, BPPARAM=BPPARAM)
    }

    # Computing rescaling factors.
    sfs <- .rescale_size_factors(stats$averages, stats$size.factors, min.mean=min.mean)
    extra.norm.args <- list(assay.type=assay.type, center.size.factors=FALSE)
    needs.subset <- !normalize.all && !is.null(subset.row)

    if (preserve.single) {
        sce <- batches[[1]]
        if (needs.subset) {
            sce <- sce[subset.row,]
        }

        all.sf <- unlist(sfs, use.names=FALSE)
        reorder <- order(unlist(stats$by.batch))
        sizeFactors(sce) <- all.sf[reorder]

        output <- do.call(logNormCounts, c(list(x=sce), c(norm.args, extra.norm.args)))

        if (!is.null(as.altexp)) {
            altExp(originals[[1]], as.altexp) <- output
            output <- originals[[1]]
        }
    } else {
        output <- vector("list", length(batches))
        names(output) <- names(batches)

        for (i in seq_along(batches)) {
            current <- batches[[i]]
            sizeFactors(current) <- sfs[[i]]
            if (needs.subset) {
                current <- current[subset.row,,drop=FALSE]
            }
            current <- do.call(logNormCounts, c(list(x=current), c(norm.args, extra.norm.args)))
            output[[i]] <- current
        }

        if (!is.null(as.altexp)) {
            for (i in seq_along(output)) {
                altExp(originals[[i]], as.altexp) <- output[[i]]
            }
            output <- originals
        }
    }

    output
}

#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
.compute_batch_statistics_list <- function(batches, assay.type, subset.row, BPPARAM) {
    empty <- vector("list", length(batches))
    names(empty) <- names(batches)
    size.factors <- all.averages <- empty

    for (b in seq_along(batches)) {
        current <- batches[[b]]

        stats <- .compute_batch_statistics(
            assay(current, assay.type), 
            sf=sizeFactors(current),
            subset.row=subset.row, 
            BPPARAM=BPPARAM
        )

        size.factors[[b]] <- stats$size.factors
        all.averages[[b]] <- stats$averages
    }

    list(size.factors=size.factors, averages=all.averages)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom BiocGenerics sizeFactors
.compute_batch_statistics_single <- function(x, batch, assay.type, subset.row, BPPARAM) {
    by.batch <- split(seq_along(batch), batch)
    empty <- vector("list", length(by.batch))
    names(empty) <- names(by.batch)
    size.factors <- all.averages <- empty

    mat <- assay(x, assay.type)
    sf <- sizeFactors(x)

    for (b in seq_along(by.batch)) {
        current <- by.batch[[b]]

        stats <- .compute_batch_statistics(
            mat[,current,drop=FALSE],
            sf=sf[current],
            subset.row=subset.row, 
            BPPARAM=BPPARAM
        )

        size.factors[[b]] <- stats$size.factors
        all.averages[[b]] <- stats$averages
    }

    list(size.factors=size.factors, averages=all.averages, by.batch=by.batch)
}

#' @importFrom scuttle librarySizeFactors calculateAverage
.compute_batch_statistics <- function(mat, sf, subset.row, BPPARAM) {
    if (is.null(sf)) {
        sf <- librarySizeFactors(mat, subset.row=subset.row, BPPARAM=BPPARAM)
    } else {
        sf <- sf / mean(sf)
    }
    ave <- calculateAverage(mat, subset.row=subset.row, size.factors=sf, BPPARAM=BPPARAM)
    list(size.factors=sf, averages=ave)
}

#' @importFrom stats median
.rescale_size_factors <- function(ave.list, sf.list, min.mean) 
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    nbatches <- length(ave.list)
    collected.ratios <- matrix(1, nbatches, nbatches)

    for (first in seq_len(nbatches-1L)) {
        first.ave <- ave.list[[first]]
        first.sum <- sum(first.ave)

        for (second in first + seq_len(nbatches - first)) {
            second.ave <- ave.list[[second]]
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
    rescaling <- collected.ratios[,smallest]

    for (idx in seq_len(nbatches)) {
        sf.list[[idx]] <- sf.list[[idx]] / rescaling[idx]
    }
    
    sf.list
}
