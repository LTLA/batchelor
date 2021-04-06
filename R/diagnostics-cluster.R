#' Cluster-based correction diagnostics
#'
#' A variety of correction diagnostics that make use of clustering information,
#' usually obtained by clustering on cells from all batches in the corrected data.
#' 
#' @param x A factor or vector specifying the assigned cluster for each cell in each batch in the corrected data.
#' Alternatively, a matrix or table containing the number of cells in each cluster (row) and batch (column).
#' @param batch A factor or vector specifying the batch of origin for each cell.
#' Ignored if \code{x} is a matrix or table.
#' @param pseudo.count A numeric scalar containing the pseudo-count to use for the log-transformation.
#' 
#' @return
#' For \code{clusterAbundanceTest}, a named numeric vector of p-values from applying Pearson's chi-squared test on each cluster.
#'
#' For \code{clusterAbundanceVar}, a named numeric vector of variances of log-abundances across batches for each cluster.
#'
#' @details
#' For \code{clusterAbundanceTest}, the null hypothesis for each cluster is that the distribution of cells across batches is proportional to the total number of cells in each batch.
#' We then use \code{\link{chisq.test}} to test for deviations from the expected proportions, possibly indicative of imperfect mixing across batches.
#' This works best for technical replicates where the population composition should be identical across batches.
#' However, the interpretation of the p-value loses its meaning for experiments where there is more biological variability between batches.
#'
#' For \code{clusterAbundanceVar}, we compute log-normalized abundances for each cluster using \code{\link{normalizeCounts}}.
#' We then compute the variance of the log-abundances across batches for each cluster.
#' Large variances indicate that there are strong relative differences in abundance across batches, indicative of either imperfect mixing or genuine batch-specific subpopulations.
#' The idea is to rank clusters by their variance to prioritize them for manual inspection to decide between these two possibilities. 
#' We use a large \code{pseudo.count} by default to avoid spuriously large variances when the counts are low.
#'
#' @author Aaron Lun
#'
#' @examples
#' set.seed(1000)
#' means <- 2^rgamma(1000, 2, 1)
#' A1 <- matrix(rpois(10000, lambda=means), ncol=50) # Batch 1 
#' A2 <- matrix(rpois(10000, lambda=means*runif(1000, 0, 2)), ncol=50) # Batch 2
#'
#' B1 <- log2(A1 + 1)
#' B2 <- log2(A2 + 1)
#' out <- fastMNN(B1, B2) 
#'
#' cluster1 <- kmeans(t(B1), centers=10)$cluster
#' cluster2 <- kmeans(t(B2), centers=10)$cluster
#' merged.cluster <- kmeans(reducedDim(out, "corrected"), centers=10)$cluster
#'
#' # Low p-values indicate unexpected differences in abundance.
#' clusterAbundanceTest(paste("Cluster", merged.cluster), out$batch)
#'
#' # High variances indicate differences in normalized abundance.
#' clusterAbundanceVar(paste("Cluster", merged.cluster), out$batch)
#'
#' @name diagnostics-cluster
NULL

#' @export
#' @rdname diagnostics-cluster
#' @importFrom stats chisq.test
clusterAbundanceTest <- function(x, batch) {
    tab <- .create_abundance_table(x, batch)
    chi.prop <- colSums(tab)/sum(tab)
    chi.results <- apply(tab, 1, FUN=chisq.test, p=chi.prop)
    vapply(chi.results, "[[", i="p.value", 0)
}

.create_abundance_table <- function(x, batch) {
    if (!is.matrix(x)) {
        x <- table(x, batch)        
        x <- unclass(x)
        dimnames(x) <- unname(dimnames(x))
    }
    x
}

#' @export
#' @rdname diagnostics-cluster
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom scuttle normalizeCounts
clusterAbundanceVar <- function(x, batch, pseudo.count=10) {
    tab <- .create_abundance_table(x, batch)
    norm <- normalizeCounts(tab, pseudo.count=pseudo.count)
    out <- rowVars(norm)
    names(out) <- rownames(norm)
    out
}
