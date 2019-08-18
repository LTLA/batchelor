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

############################################

.restricted_mnn <- function(left.data, left.restrict, right.data, right.restrict, ...) {
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

    pairs <- findMutualNN(left.data, right.data, ...)
    pairs$first <- .unrestrict_indices(pairs$first, left.restrict)
    pairs$second <- .unrestrict_indices(pairs$second, right.restrict)
    pairs
}

############################################
# Tree searching for fastMNN(). This requires some care as the
# search has to orthogonalize each batch before comparing them.

.search_for_merge <- function(remainders, ...) {
    left.idx <- right.idx <- n <- list()
    counter <- 1L

    for (i in seq_along(remainders)) {
        left <- remainders[[i]]
        left.data <- .get_node_data(left)
        left.restrict <- .get_node_restrict(left)
        left.extras <- .get_node_extras(left)

        for (j in seq_len(i-1L)) {
            right <- remainders[[j]]
            right.data <- .get_node_data(right)
            right.restrict <- .get_node_restrict(right)
            right.extras <- .get_node_extras(right)

            right.data <- .orthogonalize_other(right.data, right.restrict, left.extras)
            left.data <- .orthogonalize_other(left.data, left.restrict, right.extras)
            collected <- .restricted_mnn(left.data, left.restrict, right.data, right.restrict, ...)

            left.idx[[counter]] <- i
            right.idx[[counter]] <- j
            n[[counter]] <- length(collected$first)
            counter <- counter + 1L
        }
    }
    
    n <- unlist(n)
    chosen <- which.max(n)
    chosen.left <- left.idx[[chosen]]
    chosen.right <- right.idx[[chosen]]
    list(left=remainders[[chosen.left]], right=remainders[[chosen.right]], 
        chosen=c(chosen.left, chosen.right))
}

.update_remainders <- function(remainders, chosen, ...) {
    leftovers <- remainders[-chosen]
    val <- MNN_treenode(...)
    if (length(leftovers)) {
        c(leftovers, list(val))
    } else {
        val # Return node directly to indicate completion.
    }
}
