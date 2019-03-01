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
#' @param separate.spikes Logical scalar indicating whether spike-in size factors should be rescaled separately from endogenous genes.
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
#' Only genes with library size-adjusted average counts greater than \code{min.mean} will be used for computing the rescaling factors.
#' This improves precision and avoids problems with discreteness.
#' Users can also set \code{subset.row} to restrict the set of genes used for computing the rescaling factors.
#' However, this only affects the rescaling of the size factors - normalized values for \emph{all} genes will still be returned.
#' 
#' @section Handling spike-ins:
#' Spike-in transcripts should be either absent in all batches or, if present, they should be the same across all batches.
#' Rows annotated as spike-in transcripts are not used to compute the rescaling factors for endogenous genes.
#' 
#' By default, the spike-in size factors are rescaled using the same scaling factor for the endogenous genes in the same batch.
#' This preserves the abundances of the spike-in transcripts relative to the endogenous genes, 
#' which is important if the returned SingleCellExperiments are to be used to model technical noise.
#' 
#' If \code{separate.spikes=TRUE}, spike-in size factors are rescaled separately from those of the endogenous genes.
#' This will eliminate differences in spike-in quantities across batches at the cost of losing the ability to compare between endogenous and spike-in transcripts within each batch.
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
#' @importFrom SingleCellExperiment isSpike spikeNames
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), min.mean=1, subset.row=NULL, separate.spikes=FALSE)
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

    # Centering the endogenous size factors, and rescaling all spike-in size factors for the ride.
    for (b in seq_along(batches)) {
        current <- batches[[b]]

        cursf <- sizeFactors(current)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current, exprs_values=assay.type, subset_row=subset.row)
            warning(sprintf("no endogenous size factors in batch %i, using library sizes instead", b))
        }

        centering <- mean(cursf)
        sizeFactors(current) <- cursf / centering

        for (spike in spikeNames(current)) {
            spsf <- sizeFactors(current, spike)
            if (is.null(spsf)) {
                spsf <- librarySizeFactors(current, exprs_values=assay.type, subset_row=isSpike(current, spike))
                warning(sprintf("no %s size factors in batch %i, using library sizes instead", spike, b))
            }

            # If we want to rescale spike-ins separately, they'll need to be centered as well.
            if (separate.spikes) {
                spsf <- spsf/mean(spsf)
            } else {
                spsf <- spsf/centering
            }

            sizeFactors(current, spike) <- spsf
        }

        batches[[b]] <- current
    }

    # Adjusting size factors for the non-spike-in transcripts.
    nonspike.subset <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes=FALSE)
    if (is.null(nonspike.subset) || any(nonspike.subset)) { # works for logical or integer
        endog.rescale <- .compute_batch_rescaling(batches, subset.row=nonspike.subset, assay.type=assay.type, min.mean=min.mean)
        batches <- .rescale_size_factors(batches, endog.rescale)
    } else {
        endog.rescale <- rep(1, length(batches))
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

        if (separate.spikes) {
            spike.subset <- intersect(which(ref.spike), subset.row)
            if (length(spike.subset)) { 
                spike.rescale <- .compute_batch_rescaling(batches, subset.row=spike.subset, 
                    assay.type=assay.type, sf.type=spike, min.mean=min.mean)
            }
        } else {
            spike.rescale <- endog.rescale
        }
        batches <- .rescale_size_factors(batches, spike.rescale, sf.type=spike)
    }

    # Applying the normalization.   
    mapply(FUN=normalize, batches, MoreArgs=c(list(exprs_values=assay.type, centre_size_factors=FALSE), norm.args), SIMPLIFY=FALSE)
}

#' @importFrom scater calcAverage 
#' @importFrom stats median
#' @importFrom BiocGenerics sizeFactors 
.compute_batch_rescaling <- function(batches, subset.row, assay.type, min.mean, sf.type=NULL)
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    nbatches <- length(batches)
    collected.ave <- vector("list", nbatches)

    for (b in seq_along(batches)) {
        collected.ave[[b]] <- calcAverage(assay(batches[[b]], assay.type, withDimnames=FALSE), 
            subset_row=subset.row, use_size_factors=sizeFactors(batches[[b]], sf.type))
     }

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

    # Rescaling to the lowest coverage batch.
    smallest <- which.min(apply(collected.ratios, 2, min, na.rm=TRUE))
    collected.ratios[,smallest]
}

#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importFrom scater librarySizeFactors
.rescale_size_factors <- function(batches, rescaling, assay.type, sf.type=NULL) 
# Applying the chosen rescaling.
{
    for (idx in seq_along(batches)) {
        current <- batches[[idx]]
        cursf <- sizeFactors(current, sf.type)

        # Adjusting the size factors for the new reference.
        cursf <- cursf / rescaling[idx] 
        sizeFactors(current, sf.type) <- cursf
        batches[[idx]] <- current
    }

    batches
}

