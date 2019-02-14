#' Per-batch scaling normalization
#'
#' Perform scaling normalization within each batch to provide comparable results to the lowest-coverage batch.
#' 
#' @param ... Two or more SingleCellExperiment objects containing counts and size factors.
#' Each object is assumed to represent one batch.
#' @param assay.type A string specifying which assay values contains the counts.
#' @param norm.args A named list of further arguments to pass to \code{\link[scater]{normalize}}.
#' @param min.mean A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.
#' @param subset.row A vector specifying which features to use for correction. 
#' 
#' @details
#' When performing integrative analyses of multiple batches, it is often the case that different batches have large differences in coverage.
#' This function removes systematic differences in coverage across batches to simplify downstream comparisons.
#' It does so by resaling the size factors using median-based normalization on the ratio of the average counts between batches.
#' This is roughly equivalent to the between-cluster normalization described by Lun et al. (2016).
#' 
#' This function will adjust the size factors so that counts in high-coverage batches are scaled \emph{downwards} to match the coverage of the most shallow batch.
#' The \code{\link[scater]{normalize}} function will then add the same pseudo-count to all batches before log-transformation.
#' By scaling downwards, we favour stronger squeezing of log-fold changes from the pseudo-count, mitigating any technical differences in variance between batches.
#' Of course, genuine biological differences will also be shrunk, but this is less of an issue for upregulated genes with large counts.
#' 
#' This function is preferred over running \code{\link{normalize}} directly when computing log-normalized values for use in \code{\link{mnnCorrect}} or \code{\link{fastMNN}}.
#' In most cases, size factors will be computed within each batch;
#' their direct application in \code{\link{normalize}} will not account for scaling differences between batches.
#' In contrast, \code{multiBatchNorm} will rescale the size factors so that they are comparable across batches.
#' 
#' If spike-in transcripts are present, these should be the same across all batches.
#' Spike-in size factors are rescaled separately from those of the endogenous genes, to reflect differences in spike-in quantities across batches.
#' Conversely, spike-in transcripts are not used to compute the rescaling factors for endogenous genes.
#' 
#' Only genes with library size-adjusted average counts greater than \code{min.mean} will be used for computing the rescaling factors.
#' This improves precision and avoids problems with discreteness.
#' Users can also set \code{subset.row} to restrict the set of genes used for computing the rescaling factors.
#' However, this only affects the rescaling of the size factors - normalized values for \emph{all} genes will still be returned.
#' 
#' @return
#' A list of SingleCellExperiment objects with normalized log-expression values in the \code{"logcounts"} assay (depending on values in \code{norm.args}).
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{mnnCorrect}} and \code{\link{fastMNN}} for methods that can benefit from rescaling.
#'
#' \code{\link[scater]{normalize}} for the calculation of log-transformed normalized expression values.
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
#' Lun ATL, Bach K and Marioni JC (2016).
#' Pooling across cells to normalize single-cell RNA sequencing data with many zero counts.
#' \emph{Genome Biol.} 17:75
#'
#' @export
#' @importFrom BiocGenerics normalize
#' @importMethodsFrom scater normalize
#' @importFrom scater centreSizeFactors
#' @importFrom SingleCellExperiment isSpike spikeNames
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), min.mean=1, subset.row=NULL)
# Performs multi-batch normalization, adjusting for differences 
# in scale between SCE objects supplied in '...'.
# 
# written by Aaron Lun
# created 4 June 2018
{
    batches <- list(...)
    checkBatchConsistency(batches)
    if (!checkIfSCE(batches) || length(batches)==0L) {
        stop("at least one SingleCellExperiment object must be supplied")
    }
    checkSpikeConsistency(batches)

    # Centering the size factors, to avoid large rescaling factors due to systematically larger size factors in the input.
    # This would not be reflective of the actual magnitude of the counts.
    batches <- lapply(batches, centreSizeFactors)

    # Adjusting size factors for the non-spike-in transcripts.
    nonspike.subset <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes=FALSE)
    if (is.null(nonspike.subset) || any(nonspike.subset)) { # works logical or integer
        batches <- .batch_rescaler(batches, subset.row=nonspike.subset, exprs_values=assay.type, sf.type=NULL, min.mean=min.mean)
    }

    # Adjusting size factors for each of the spike-ins.
    ref.spike.names <- spikeNames(batches[[1]])
    subset.row <- .row_subset_to_index(batches[[1]], subset.row)

    for (spike in ref.spike.names) {
        ref.spike <- isSpike(batches[[1]], spike)
        for (b in seq_along(batches)) {
            if (!identical(ref.spike, isSpike(batches[[b]], spike))) {
                stop(sprintf("%s spike-in identities differ across batches", spike))
            }
        }

        spike.subset <- intersect(which(ref.spike), subset.row)
        if (length(spike.subset)) { 
            batches <- .batch_rescaler(batches, subset.row=spike.subset, exprs_values=assay.type, sf.type=spike, min.mean=min.mean)
        }
    }

    # Applying the normalization.   
    mapply(FUN=normalize, batches, MoreArgs=c(list(exprs_values=assay.type, centre_size_factors=FALSE), norm.args), SIMPLIFY=FALSE)
}

#' @importFrom scater calcAverage librarySizeFactors
#' @importFrom stats median
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
.batch_rescaler <- function(batches, subset.row, exprs_values, sf.type, min.mean) 
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    collected.ave <- lapply(batches, calcAverage, use_size_factors=TRUE, subset_row=subset.row, exprs_values=exprs_values)
    nbatches <- length(batches)
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

            curratio <- median(second.ave[keep]/first.ave[keep])
            if (!is.finite(curratio) || curratio==0) {
                stop("median ratio of averages between batches is not finite")
            }

            collected.ratios[first,second] <- curratio
            collected.ratios[second,first] <- 1/curratio
        }
    }

    # Rescaling to the lowest coverage.
    smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
    
    for (idx in seq_along(batches)) {
        current <- batches[[idx]]
        cursf <- sizeFactors(current, sf.type)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current, exprs_values=exprs_values, subset_row=subset.row)
            warning(sprintf("no %ssize factors in batch %i, using library sizes instead",
                if (is.null(sf.type)) "" else paste(sf.type, ""), idx))
        }
    
        # Adjusting the size factors for the new reference.
        cursf <- cursf/mean(cursf) / collected.ratios[idx,smallest]
        sizeFactors(current, sf.type) <- cursf
        batches[[idx]] <- current
    }

    return(batches)
}

