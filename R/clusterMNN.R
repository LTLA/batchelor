#' @export
#' @importFrom scater .bpNotSharedOrUp
#' @importFrom BiocParallel bpstart bpstop
clusterMNN <- function(..., batch=NULL, restrict=NULL, clusters=NULL, d=50, weights=NULL,
    k=1, prop.k=0, ndist=3, merge.order=NULL, auto.merge=FALSE, min.batch.skip=0,
    subset.row=NULL, correct.all=FALSE, assay.type="logcounts",
    BSPARAM=IrlbaParam(deferred=TRUE), BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    # Setting up the parallelization environment.
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    pca <- multiBatchPCA(..., batch=batch, BSPARAM=BSPARAM, BPPARAM=BPPARAM, 
        d=d, weights=weights, subset.row=subset.row, get.all.genes=correct.all)

    .cluster_mnn(pca, batch=batch, restrict=restrict, clusters=clusters,
        k=k, prop.k=prop.k, ndist=ndist, 
        merge.order=merge.order, auto.merge=auto.merge, min.batch.skip=min.batch.skip,
        BNPARAM=BNPARAM, BPPARAM=BPPARAM)
}

#' @export
reducedClusterMNN <- function(..., batch=NULL, restrict=NULL, clusters=NULL, 
    k=1, prop.k=0, ndist=3, merge.order=NULL, auto.merge=FALSE, min.batch.skip=0,
    BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    batches <- list(...)
    checkBatchConsistency(batches, cells.in.columns=FALSE)
    if (length(batches)==1L) {
        pca <- divideIntoBatches(batches[[1]], batch=batch, byrow=FALSE)
    } else {
        pca <- batches
    }

    .cluster_mnn(batches, batch=batch, restrict=restrict, clusters=clusters,
        k=k, prop.k=prop.k, ndist=ndist, 
        merge.order=merge.order, auto.merge=auto.merge, min.batch.skip=min.batch.skip,
        BNPARAM=BNPARAM, BPPARAM=BPPARAM)
}


#' @importFrom utils tail
#' @importFrom BiocNeighbors queryKNN
#' @importFrom S4Vectors DataFrame
.cluster_mnn <- function(batches, batch, restrict, clusters, ..., BNPARAM, BPPARAM) {
    # Checking input clusters.
    if (is.null(clusters)) {
        clusters <- .generate_kmeans_clusters(batches, restrict, BNPARAM, BPPARAM)
    } else if (!is.null(batch)) {
        clusters <- split(clusters, batch)
    } else {
        if (length(clusters)!=length(batches)) {
            stop("'...' and 'clusters' should be of the same length")
        }
        for (i in seq_along(clusters)) {
            if (nrow(batches[[i]])!=length(clusters[[i]])) {
                stop("corresponding entries of '...' and 'clusters' should have the same number of cells")
            }
        }
    }

    # Performing MNN on the cluster centroids, which is relatively easy.
    clusters <- lapply(clusters, as.character)
    centers <- vector("list", length(clusters))
    for (i in seq_along(batches)) {
        if (!is.null(curres <- restrict[[i]])) {
            cls <- clusters[[i]][curres]
            mat <- rowsum(batches[[i]][curres,,drop=FALSE], cls)
        } else {
            cls <- clusters[[i]]
            mat <- rowsum(batches[[i]], cls)
        }
        centers[[i]] <- mat/as.numeric(table(cls)[rownames(mat)])
    }

    output <- do.call(reducedMNN, c(centers, list(..., BNPARAM=BNPARAM, BPPARAM=BPPARAM)))

    # Smoothing the correction applied to the clusters and applying it to the cells.
    # Each cell's smoothing bandwidth is set to the distance from its assigned cluster.
    last <- 0L
    renamed <- vector("list", length(batches))

    for (i in seq_along(batches)) {
        curbatch <- batches[[i]]
        curcenter <- centers[[i]]
        indices <- last + seq_len(nrow(curcenter))
        last <- tail(indices, 1L)

        corrected <- output$corrected[indices,,drop=FALSE]
        delta <- corrected - curcenter
#        closest <- queryKNN(query=curbatch, X=curcenter, k=1, BNPARAM=BNPARAM, 
#            BPPARAM=BPPARAM, get.distance=FALSE)$index
        sigma <- queryKNN(query=curbatch, X=curcenter, k=1, BNPARAM=BNPARAM, 
            BPPARAM=BPPARAM, get.distance=TRUE)$distance[,1]

        adj <- 0
        total <- 0
        tbatch <- t(curbatch)
        for (j in seq_len(nrow(curcenter))) {
            d2 <- colSums((tbatch - curcenter[j,])^2)
            weight <- exp(- d2/sigma^2) # no need to worry about underflow, there is always a weight of exp(-1) somewhere.
            adj <- adj + outer(weight, delta[j,])
            total <- total + weight
        }
    
        batches[[i]] <- batches[[i]] + adj/total
#        batches[[i]] <- batches[[i]] + delta[closest,,drop=FALSE]
        renamed[[i]] <- rep(output$batch[indices[1]], nrow(curbatch))
    }

    DataFrame(corrected=I(do.call(rbind, batches)), batch=unlist(renamed))
}

#' @importFrom stats kmeans
#' @importFrom BiocNeighbors queryKNN
.generate_kmeans_clusters <- function(batches, restrict, BNPARAM, BPPARAM) {
    clusters <- list()
    for (i in seq_along(batches)) {
        curbatch <- batches[[i]]
        if (is.null(curres <- restrict[[i]])) {
            clusters[[i]] <- kmeans(curbatch, centers=sqrt(nrow(curbatch)))$cluster

        } else {
            subbatch <- curbatch[curres,,drop=FALSE]
            k.out <- kmeans(subbatch, centers=sqrt(nrow(subbatch)))
            full.out <- integer(nrow(curbatch))
            full.out[curres] <- k.out$cluster

            leftovers <- !seq_len(nrow(curbatch)) %in% curres
            closest <- queryKNN(X=k.out$centers, query=curbatch[leftovers,,drop=FALSE], k=1, 
                BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)$index
            full.out[leftovers] <- closest
            clusters[[i]] <- full.out
        }
    }
    clusters
}
