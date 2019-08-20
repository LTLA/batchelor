#' Per-batch scaling normalization
#'
#' Perform scaling normalization within each batch to provide comparable results to the lowest-coverage batch.
#' 
#' @param ... Two or more \linkS4class{SingleCellExperiment} objects containing counts and size factors.
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
#' The \code{\link{logNormCounts}} function will then add the same pseudo-count to all batches before log-transformation.
#' By scaling downwards, we favour stronger squeezing of log-fold changes from the pseudo-count, mitigating any technical differences in variance between batches.
#' Of course, genuine biological differences will also be shrunk, but this is less of an issue for upregulated genes with large counts.
#' 
#' For comparison, consider running \code{\link{logNormCounts}} separately in each batch prior to correction.
#' In most cases, size factors will be computed within each batch;
#' batch-specific application in \code{\link{logNormCounts}} will not account for scaling differences between batches.
#' In contrast, \code{multiBatchNorm} will rescale the size factors so that they are comparable across batches.
#' This removes at least one difference between batches to facilitate easier correction.
#' 
#' Only genes with library size-adjusted average counts greater than \code{min.mean} will be used for computing the rescaling factors.
#' This improves precision and avoids problems with discreteness.
#' Users can also set \code{subset.row} to restrict the set of genes used for computing the rescaling factors.
#' However, this only affects the rescaling of the size factors - normalized values for \emph{all} genes will still be returned.
#'
#' @section Note about spike-ins:
#' Rescaling is only performed on endogenous genes in each SingleCellExperiment object.
#' If any spike-in transcripts are present in the \code{\link{altExps}},
#' their abundances will not be rescaled here, and are no longer directly comparable to the rescaled abundances of the genes.
#' This is usually not a major problem as spike-ins are rarely used during the batch correction itself -
#' however, users should not attempt to perform variance modelling with the spike-ins on the output of this function.
#' 
#' @return
#' A list of SingleCellExperiment objects with normalized log-expression values in the \code{"logcounts"} assay (depending on values in \code{norm.args}).
#' 
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{mnnCorrect}} and \code{\link{fastMNN}}, for methods that can benefit from rescaling.
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
#' Lun ATL (2018).
#' Further MNN algorithm development.
#' \url{https://MarioniLab.github.io/FurtherMNN2018/theory/description.html}
#'
#' @export
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
#' @importMethodsFrom scater logNormCounts librarySizeFactors
multiBatchNorm <- function(..., assay.type="counts", norm.args=list(), min.mean=1, subset.row=NULL) {
    batches <- list(...)
    checkBatchConsistency(batches)
    if (length(batches)==0L) { 
        stop("at least one SingleCellExperiment must be supplied") 
    }

    # Centering the endogenous size factors.
    for (b in seq_along(batches)) {
        current <- batches[[b]]

        cursf <- sizeFactors(current)
        if (is.null(cursf)) {
            cursf <- librarySizeFactors(current, exprs_values=assay.type, subset_row=subset.row)
        }

        centering <- mean(cursf)
        sizeFactors(current) <- cursf / centering
        batches[[b]] <- current
    }

    # Adjusting size factors.
    nonspike.subset <- .SCE_subset_genes(subset.row, batches[[1]], get.spikes=FALSE)
    endog.rescale <- .compute_batch_rescaling(batches, subset.row=nonspike.subset, assay.type=assay.type, min.mean=min.mean)
    batches <- .rescale_size_factors(batches, endog.rescale)

    # Applying the normalization.   
    mapply(FUN=logNormCounts, batches, 
        MoreArgs=c(list(exprs_values=assay.type, center_size_factors=FALSE), norm.args), 
        SIMPLIFY=FALSE)
}

#' @importFrom scater calculateAverage 
#' @importFrom stats median
#' @importFrom BiocGenerics sizeFactors sizeFactors<-
.compute_batch_rescaling <- function(batches, subset.row, assay.type, min.mean) 
# Computes the median ratios (a la DESeq normalization),
# finds the smallest ratio and uses that as the reference.
{
    nbatches <- length(batches)
    collected.ave <- lapply(batches, FUN=calculateAverage, exprs_values=assay.type, subset_row=subset.row)

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
.rescale_size_factors <- function(batches, rescaling, assay.type)
# Applying the chosen rescaling.
{
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

