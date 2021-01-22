#' Cluster-based MNN
#'
#' Perform MNN correction based on cluster centroids, 
#' using the corrected centroid coordinates to correct the per-cell expression values 
#' with a variable bandwidth Gaussian kernel.
#'
#' @inheritParams fastMNN
#' @param clusters A list of length equal to \code{...} containing the assigned cluster for each cell in each batch.
#' Alternatively, a BlusterParam object from the \pkg{bluster} package, specifying the clustering to be applied to each batch.
#' @param cluster.d Integer scalar indicating how many PCs should be used in clustering.
#' Only used if \code{clusters} is not a list.
#' If \code{NA}, no PCA is performed and clustering is applied directly to the (cosine-normalized) gene expression values.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.
#' Only used if \code{clusters} is not a list.
#' 
#' @return
#' A \linkS4class{SingleCellExperiment} containing per-cell expression values where each row is a gene and each column is a cell.
#' This has the same format as the output of \code{\link{fastMNN}}
#' but with an additional \code{cluster} field in the \code{\link{colData}} containing the cluster identity of each cell.
#' The \code{\link{metadata}} contains:
#' \itemize{
#' \item \code{merge.info}, a DataFrame with the same format as the output of \code{\link{fastMNN}}.
#' However, the \code{pairs} and \code{lost.var} refer to the cluster centroids, not the cells.
#' \item \code{clusters}, a DataFrame with one row for each cluster across all batches in \code{...}.
#' This can be row-indexed by the values in \code{pairs} to determine the identity of the clusters in each MNN pair.
#' An additional \code{meta} column is provided that describes the meta-cluster to which each cluster belongs.
#' }
#'
#' @details
#' These functions are motivated by the scenario where each batch has been clustered separately
#' and each cluster has already been annotated with some meaningful biological state.
#' We want to identify which biological states match to each other across batches;
#' this is achieved by identifying mutual nearest neighbors based on the cluster centroids with \code{\link{reducedMNN}}.
#' 
#' MNN pairs are identified with \code{k=1} to ensure that each cluster has no more than one match in another batch.
#' This reduces the risk of inadvertently merging together different clusters from the same batch.
#' By comparison, higher values of \code{k} may result in many-to-one mappings between batches
#' such that the correction will implicitly force different clusters together.
#' 
#' Using this guarantee of no-more-than-one mappings across batches,
#' we can form meta-clusters by identifying all components of the resulting MNN graph.
#' Each meta-cluster can be considered to represent some biological state (e.g., cell type),
#' and all of its constituents are the matching clusters within each batch.
#'
#' As an extra courtesy, \code{clusterMNN} will also compute corrected values for each cell.
#' This is done by applying a Gaussian kernel to the correction vectors for the centroids,
#' where the bandwidth is proportional to the distance between that cell and the closest cluster centroid.
#' This yields a smooth correction function that avoids edge effects at cluster boundaries.
#'
#' If \code{clusters} is set to a BlusterParam object (see the \pkg{bluster} package),
#' a PCA is performed in each batch with the specified \code{BSPARAM}.
#' The PCs are then used in clustering with \code{\link[bluster]{clusterRows}} to obtain a list of clusters.
#' This can be used to mimic per-cell batch correction in the absence of \emph{a priori} clusters.
#'
#' @inheritSection fastMNN Dealing with alternative Experiments
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up some data for multiple batches:
#' means <- matrix(rnorm(3000), ncol=3)
#' colnames(means) <- LETTERS[1:3]
#'
#' B1 <- means[,sample(LETTERS[1:3], 500, replace=TRUE)]
#' B1 <- B1 + rnorm(length(B1))
#'
#' B2 <- means[,sample(LETTERS[1:3], 500, replace=TRUE)]
#' B2 <- B2 + rnorm(length(B2)) + rnorm(nrow(B2)) # batch effect.
#'
#' # Applying the correction with some made-up clusters:
#' cluster1 <- kmeans(t(B1), centers=10)$cluster
#' cluster2 <- kmeans(t(B2), centers=10)$cluster
#' out <- clusterMNN(B1, B2, clusters=list(cluster1, cluster2)) 
#'
#' rd <- reducedDim(out, "corrected") 
#' plot(rd[,1], rd[,2], col=out$batch)
#'
#' # Obtaining the clusters internally.
#' out2 <- clusterMNN(B1, B2, clusters=bluster::NNGraphParam())
#' rd2 <- reducedDim(out2, "corrected") 
#' plot(rd2[,1], rd2[,2], col=out$batch)
#' 
#' @references
#' Lun ATL (2019).
#' Cluster-based mutual nearest neighbors correction 
#' \url{https://marionilab.github.io/FurtherMNN2018/theory/clusters.html}
#' 
#' Lun ATL (2019).
#' A discussion of the known failure points of the fastMNN algorithm.
#' \url{https://marionilab.github.io/FurtherMNN2018/theory/failure.html}
#'
#' @seealso
#' \code{\link{reducedMNN}}, which is used internally to perform the correction.
#'
#' @export
#' @importFrom scuttle .bpNotSharedOrUp .unpackLists
#' @importFrom BiocParallel bpstart bpstop SerialParam
#' @importFrom BiocSingular IrlbaParam ExactParam
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom BiocGenerics rbind
#' @importFrom igraph make_graph components
#' @importFrom SummarizedExperiment assay
clusterMNN <- function(..., batch=NULL, restrict=NULL, clusters, cluster.d=50,
    cos.norm=TRUE, merge.order=NULL, auto.merge=FALSE, min.batch.skip=0,
    subset.row=NULL, correct.all=FALSE, assay.type="logcounts", as.altexp=NULL,
    BSPARAM=IrlbaParam(), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    batches <- .unpackLists(...)
    if (!is.null(as.altexp)) {
        batches <- lapply(batches, altExp, e=as.altexp)
    }
    checkBatchConsistency(batches, cells.in.columns=TRUE)
    restrict <- checkRestrictions(batches, restrict, cells.in.columns=TRUE)

    # Extracting information from SCEs.
    is.sce <- checkIfSCE(batches)
    if (any(is.sce)) {
        batches[is.sce] <- lapply(batches[is.sce], assay, i=assay.type)
    }

    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    if (unified <- (length(batches)==1L)) {
        if (is.null(batch)) {
            stop("'batch' must be specified if '...' has only one object")
        }

        divided <- divideIntoBatches(x=batches[[1]], restrict=restrict[[1]], batch=batch, byrow=FALSE)
        restrict <- divided$restrict
        batches <- divided$batches

        if (is.list(clusters)) {
            if (length(clusters)!=1L) {
                stop("'clusters' must be a list of length 1 when '...' contains one element")
            }
            clusters <- split(clusters[[1]], batch)
        }
    } 

    clusters <- .format_clusters(batches=batches, clusters=clusters, 
        restrict=restrict, cluster.d=cluster.d, subset.row=subset.row,
        BPPARAM=BPPARAM, BSPARAM=BSPARAM)

    if (cos.norm) { 
        l2 <- lapply(batches, cosineNorm, mode="l2norm", subset.row=subset.row, BPPARAM=BPPARAM)
        batches <- mapply(batches, l2, FUN=.apply_cosine_norm, SIMPLIFY=FALSE)
    }

    cluster.out <- .compute_centroids(batches=batches, clusters=clusters, restrict=restrict)

    pca <- .full_rank_pca(cluster.out$centers, subset.row=subset.row, correct.all=correct.all)

    merge.out <- reducedMNN(pca, k=1, merge.order=merge.order, auto.merge=auto.merge, 
        BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    prop.out <- .propagate_to_cells(batches, restrict=restrict,
        pca=pca, after=merge.out, 
        BNPARAM=BNPARAM, BPPARAM=BPPARAM,
        correct.all=correct.all, subset.row=subset.row)

    # Formatting the output.
    output <- .convert_to_SCE(prop.out, pca)
    output$cluster <- unlist(cluster.out$clusters)

    metadata(output) <- metadata(merge.out)
    metadata(output)$cluster <- DataFrame(cluster=rownames(merge.out), batch=merge.out$batch)
    
    # Identifying clusters of clusters across batches.
    all.pairs <- unlist(metadata(output)$merge.info$pairs)
    g <- make_graph(rbind(all.pairs$left, all.pairs$right), n=nrow(merge.out), directed=FALSE)
    metadata(output)$cluster$meta <- components(g)$membership

    if (unified) {
        output <- output[,divided$reorder]
    }

    output
}

.full_rank_pca <- function(centers, subset.row, correct.all) 
# NOTE: kept as a separate function for easier testing.
{
    # PCA is strictly for allowing us to impute genes outside of subset.row,
    # and to allow us to return a low-rank matrix of per-cell expression values.
    # There is no signal:noise trade-off because we ask for the maximum 'd'.
    multiBatchPCA(centers, BSPARAM=ExactParam(), deferred=FALSE, 
        get.all.genes=correct.all, subset.row=subset.row, 
        d=sum(vapply(centers, ncol, 0L)) - 1L
    )
}

#' @importFrom Matrix t
#' @importFrom BiocSingular runPCA
.format_clusters <- function(batches, clusters, restrict, subset.row, cluster.d, BSPARAM, BPPARAM) {
    if (is.list(clusters)) {
        if (length(clusters)!=length(batches)) {
            stop("'...' and 'clusters' should be of the same length")
        }

        for (i in seq_along(clusters)) {
            if (ncol(batches[[i]])!=length(clusters[[i]])) {
                stop("corresponding entries of '...' and 'clusters' should have the same number of cells")
            }
        }
    } else {
        if (!is(clusters, "BlusterParam")) {
            stop("'clusters' must be either a list or a BlusterParam object")
        }

        # TODO: respond properly to restriction.
        tmp <- batches
        for (i in seq_along(batches)) {
            current <- batches[[i]]
            if (!is.null(subset.row)) {
                current <- current[subset.row,,drop=FALSE]
            }
            current <- t(current)

            if (!is.na(cluster.d)) {
                pc <- runPCA(current, rank=cluster.d, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
                current <- pc$x
            }
            tmp[[i]] <- bluster::clusterRows(current, BLUSPARAM=clusters)
        } 

        clusters <- tmp
    }

    clusters
}

#' @importFrom scuttle sumCountsAcrossCells
#' @importFrom SummarizedExperiment assay
.compute_centroids <- function(batches, clusters, restrict) {
    centers <- vector("list", length(clusters))
    for (i in seq_along(batches)) {
        B <- batches[[i]]
        C <- clusters[[i]]
        if (!is.null(curres <- restrict[[i]])) {
            C <- C[curres]
            B <- B[,curres,drop=FALSE]
        }
        centers[[i]] <- assay(sumCountsAcrossCells(B, C, average=TRUE))
    }

    names(centers) <- names(clusters) <- names(batches)
    list(centers=centers, clusters=clusters)
}

#' @importFrom stats median
#' @importFrom Matrix colSums crossprod t
#' @importFrom BiocNeighbors queryKNN
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom utils tail
#' @importFrom scuttle .subset2index
.propagate_to_cells <- function(batches, restrict, pca, after, subset.row, correct.all, BNPARAM, BPPARAM) {
    pca.center <- metadata(pca)$centers
    pca.rotation <- metadata(pca)$rotation

    subset.row <- .subset2index(subset.row, batches[[1]], byrow=TRUE)
    if (correct.all) {
        # Undoing the expansion to all vectors.
        pca.rotation <- pca.rotation[subset.row,,drop=FALSE]
        pca.center <- pca.center[subset.row] 
    }
    pca.adj <- drop(pca.center %*% pca.rotation)

    last <- 0L
    renamed <- vector("list", length(batches))

    for (i in seq_along(batches)) {
        curbatch <- batches[[i]][subset.row,,drop=FALSE]
        curbatch <- crossprod(curbatch, pca.rotation)
        curbatch <- t(t(curbatch) - pca.adj)

        centroids <- pca[[i]]
        indices <- last + seq_len(nrow(centroids))
        last <- tail(indices, 1L)

        corrected <- after$corrected[indices,,drop=FALSE]
        delta <- corrected - centroids
        sigma <- median(queryKNN(query=curbatch, X=centroids, k=1,
            subset=restrict[[i]],
            BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.index=FALSE)$distance[,1])
   
        batches[[i]] <- .smooth_gaussian_from_centroids(curbatch, 
            centers=centroids, sigma=sigma, delta=delta)
        renamed[[i]] <- rep(after$batch[indices[1]], nrow(curbatch))
    }

    DataFrame(corrected=I(do.call(rbind, batches)), batch=unlist(renamed))
}

#' @importFrom Matrix rowSums colSums t
.smooth_gaussian_from_centroids <- function(x, centers, sigma, delta) 
# Using a two-pass Gaussian kernel to avoid underflow.
{
    tbatch <- t(x)
    ncenters <- nrow(centers)
    ncells <- nrow(x)

    weights <- matrix(0, ncells, ncenters)
    for (j in seq_len(ncenters)) {
        weights[,j] <- -colSums((tbatch - centers[j,])^2)
    }
    weights <- weights/sigma^2

    top.weight <- weights[cbind(seq_len(ncells), max.col(weights))]
    norm.weights <- exp(weights - top.weight)
    norm.weights <- norm.weights/rowSums(norm.weights)

    # TODO: use LowRankMatrix instead.
    for (j in seq_len(ncenters)) {
        x <- x + outer(norm.weights[,j], delta[j,])
    }

    x
}
