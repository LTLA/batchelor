#' @importClassesFrom DelayedArray integer_OR_NULL
setClass("MNN_treenode", slots=c(index="integer", data="ANY", restrict="integer_OR_NULL", origin="integer", extras="list"))

MNN_treenode <- function(index, data, restrict, origin=rep(index, nrow(data)), extras=list()) {
    new("MNN_treenode", index=index, data=data, restrict=restrict, origin=origin, extras=extras)
}

.get_node_index <- function(node) node@index

.get_node_data <- function(node) node@data

.get_node_origin <- function(node) node@origin

.get_node_restrict <- function(node) node@restrict

.get_node_extras <- function(node) node@extras

############################################
# Methods for creating a tree with a pre-defined structure.

.binarize_tree <- function(merge.tree) {
    if (!is.list(merge.tree) && length(merge.tree)==1L) {
        return(merge.tree)
    }

    N <- length(merge.tree)
    if (N==1L) {
        # Get rid of useless internal nodes with only one child.
        .binarize_tree(merge.tree[[1]])
    } else if (N==2L) {
        # Keep going.
        list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
    } else if (N > 2L) {
        # Progressive merge.
        current <- list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
        for (i in 3:N) {
            current <- list(current, .binarize_tree(merge.tree[[i]]))
        }
        current
    } else {
        stop("merge tree contains a node with no children")
    }
}

.fill_tree <- function(merge.tree, batches, restrict) {
    if (!is.list(merge.tree)) {
        val <- batches[[merge.tree]]
        return(MNN_treenode(index=merge.tree, data=val, restrict=restrict[[merge.tree]]))
    }
    if (length(merge.tree)!=2L) {
        stop("merge tree structure should contain two children per node")
    }
    merge.tree[[1]] <- .fill_tree(merge.tree[[1]], batches, restrict)
    merge.tree[[2]] <- .fill_tree(merge.tree[[2]], batches, restrict)
    merge.tree
}

.get_next_merge <- function(merge.tree, path=NULL) {
    if (!is.list(merge.tree[[1]]) && !is.list(merge.tree[[2]])) {
        list(left=merge.tree[[1]], right=merge.tree[[2]], chosen=path)
    } else if (is.list(merge.tree[[2]])) {
        .get_next_merge(merge.tree[[2]], c(path, 2L))
    } else {
        .get_next_merge(merge.tree[[1]], c(path, 1L))
    }
}

.update_tree <- function(merge.tree, path, ...) {
    if (length(path)==0L) {
        return(MNN_treenode(...))
    }
    merge.tree[[path[1]]] <- .update_tree(merge.tree[[path[[1]]]], path[-1], ...)
    merge.tree
}

#' @importFrom utils tail relist
.create_tree_predefined <- function(batches, restrict, merge.order) {
    if (is.null(merge.order)) {
        merge.order <- seq_along(batches)
    }

    if (!is.list(merge.order) && length(merge.order) > 1L) {
        merge.tree <- list(merge.order[1], merge.order[2])
        for (i in tail(merge.order, -2L)) {
            merge.tree <- list(merge.tree, i)
        }
    } else {
        merge.tree <- merge.order
    }

    merge.tree <- .binarize_tree(merge.tree)

    # Checking validity of leaf identities.
    leaves <- unlist(merge.tree)
    if (!is.numeric(leaves)) {
        leaves <- match(leaves, names(batches))
    } else {
        leaves <- as.integer(leaves)
    }
    if (any(is.na(leaves)) || anyDuplicated(leaves) || any(leaves < 1) || any(leaves > length(batches))) {
        stop("invalid leaf nodes specified in 'merge.order'")
    }

    merge.tree <- relist(leaves, merge.tree)
    .fill_tree(merge.tree, batches, restrict)
}

############################################

.restricted_mnn <- function(left.data, left.restrict, right.data, right.restrict, k, prop.k=NULL, ...) {
    if (!is.null(left.restrict)) {
        left.data <- left.data[left.restrict,,drop=FALSE]
    } else {
        left.data <- left.data
    }

    if (!is.null(right.restrict)) {
        right.data <- right.data[right.restrict,,drop=FALSE]
    } else {
        right.data <- right.data
    }

    k1 <- .choose_k(k, prop.k, nrow(left.data))
    k2 <- .choose_k(k, prop.k, nrow(right.data))

    pairs <- findMutualNN(left.data, right.data, k1=k1, k2=k2, ...)
    pairs$first <- .unrestrict_indices(pairs$first, left.restrict)
    pairs$second <- .unrestrict_indices(pairs$second, right.restrict)
    pairs
}

.choose_k <- function(k, prop.k, N) {
    if (is.null(prop.k)) {
        k
    } else {
        min(N, max(k, round(prop.k * N)))
    }
}

############################################
# Tree searching for fastMNN(). This requires some care as the
# search has to orthogonalize each batch before comparing them.

#' @importFrom S4Vectors DataFrame metadata<- metadata
#' @importClassesFrom S4Vectors List
.initialize_auto_search <- function(batches, restrict, ...) {
    remainders <- lapply(seq_along(batches), function(i) {
        MNN_treenode(index=i, data=batches[[i]], restrict=restrict[[i]])
    })
    remainders <- as(remainders, "List")

    collected <- matrix(0L, length(remainders), length(remainders))
    for (i in seq_along(remainders)) {
        below <- i - 1L
        collected[i,seq_len(below)] <- .count_mnn_pairs(remainders[[i]], remainders, upto=below, ...) 
    }

    metadata(remainders)$pairwise <- collected
    remainders
}

#' @importFrom S4Vectors DataFrame
.count_mnn_pairs <- function(left, remainders, upto, ..., orthogonalize=TRUE) {
    left.data <- .get_node_data(left)
    left.restrict <- .get_node_restrict(left)
    left.extras <- .get_node_extras(left)

    n <- integer(upto)
    for (j in seq_len(upto)) {
        right <- remainders[[j]]
        right.data <- .get_node_data(right)
        right.restrict <- .get_node_restrict(right)
        right.extras <- .get_node_extras(right)

        if (orthogonalize) {
            right.data <- .orthogonalize_other(right.data, right.restrict, left.extras)
            left.data <- .orthogonalize_other(left.data, left.restrict, right.extras)
        }

        collected <- .restricted_mnn(left.data, left.restrict, right.data, right.restrict, ...)
        n[j] <- length(collected$first)
    }

    n
}

#' @importFrom S4Vectors metadata
.pick_best_merge <- function(remainders) {
    stats <- metadata(remainders)$pairwise
    chosen <- which(stats==max(stats), arr.ind=TRUE)[1,]
    chosen.left <- chosen[1] 
    chosen.right <- chosen[2] 
    list(left=remainders[[chosen.left]], right=remainders[[chosen.right]], chosen=chosen)
}

#' @importFrom S4Vectors metadata metadata<- List
.update_remainders <- function(remainders, chosen, ..., mnn.args) {
    remainders <- remainders[-chosen]
    val <- MNN_treenode(...)

    if (length(remainders)) {
        old.meta <- metadata(remainders)$pairwise

        # Discarding whatever batches just got merged, and adding the 
        # number of MNN pairs to the newly formed batch.
        old.meta <- old.meta[-chosen,-chosen,drop=FALSE]
        new.stats <- do.call(.count_mnn_pairs, c(list(val, remainders, upto=length(remainders)), mnn.args))
        new.meta <- rbind(old.meta, new.stats)
        new.meta <- cbind(new.meta, 0) # must be square, so that '-chosen' has something to remove.

        remainders <- c(remainders, List(val))
        metadata(remainders)$pairwise <- new.meta
        remainders

    } else {
        val # Return node directly to indicate completion.
    }
}
