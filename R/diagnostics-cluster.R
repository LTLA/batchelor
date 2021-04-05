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
#' @param before Factor or vector specifying the assigned cluster for each cell in the \emph{uncorrected} data (i.e., before batch correction) for a \emph{single} batch.
#' @param after Factor or vector specifying the assigned cluster for each cell in the \emph{corrected} data (i.e., after batch correction), for the same cells in \code{before}.
#' 
#' @return
#' For \code{clusterAbundanceTest}, a named numeric vector of p-values from applying Pearson's chi-squared test on each cluster.
#'
#' For \code{clusterAbundanceVar}, a named numeric vector of variances of log-abundances across batches for each cluster.
#'
#' For \code{compareMergedClusters}, a list containing:
#' \itemize{
#' \item \code{proportions}, a matrix where each row corresponds to one of the \code{after} clusters
#' and contains the (weighted) proportion of its cells that were derived from each of the \code{before} clusters.
#' \item \code{max.prop}, a numeric vector specifying the maximum proportion for each \code{after} cluster (i.e., row of \code{proportions}).
#' \item \code{which.max}, a character vector specifying the \code{before} cluster containing the maximum proportion for each \code{after} cluster.
#' }
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
#' For \code{compareMergedClusters}, the idea is to compare the clustering before and after batch correction for cells in a \emph{single} batch.
#' Correction should not remove heterogeneity within each batch, so each cluster generated from the corrected data (\code{after}) should map to a single cluster in the uncorrected data (\code{before}).
#' We quantify this by looking at each \code{after} cluster and computing the proportion of its cells that were obtained from each \code{before} cluster.
#' (This proportion is computed by weighting each cell according to its \code{before} frequency, to ensure we consider less-frequent clusters.) 
#' For each \code{after} cluster, we then identify the maximum proportion across all \code{before} clusters, i.e., \code{max.prop}.
#' Multiple mappings indicate that the different clusters in \code{before} were merged in the corrected data, representing loss of heterogeneity and manifesting as a low \code{max.prop}.
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
#' # Any low 'max.prop' are candidates for loss of heterogeneity.
#' compareMergedClusters(paste("Before", cluster1),
#'                       paste("After", merged.cluster[out$batch==1]))
#' compareMergedClusters(paste("Before", cluster2),
#'                       paste("After", merged.cluster[out$batch==2]))
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

#' @export
#' @rdname diagnostics-cluster
#' @importFrom DelayedMatrixStats rowMaxs
compareMergedClusters <- function(before, after) { 
    tab <- table(after, before)
    tab <- t(t(tab)/colSums(tab)) # Effectively weighting by 'before' abundance.
    tab <- tab/rowSums(tab)
    list(proportions=tab, max.prop=rowMaxs(tab), which.max=colnames(tab)[max.col(tab, ties.method="first")])
}
