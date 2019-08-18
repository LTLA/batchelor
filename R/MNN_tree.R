setClass("MNN_treenode", slots=c(index="integer", data="ANY", restrict="integer", origin="integer", extras="list"))

.get_node_index <- function(node) node@index

.get_node_data <- function(node) node@data

.get_node_origin <- function(node) node@origin

.get_node_restrict <- function(node) node@restrict

.get_node_extras <- function(node) node@extras

.fill_tree <- function(merge.tree, batches, restrict) {
    if (!is.list(merge.tree)) {
        val <- batches[[merge.tree]]
        return(new("MNN_treenode", 
            index=merge.tree, 
            data=val,
            restrict=restrict[[merge.tree]],
            origin=rep(merge.tree, nrow(val)),
            extras=list()
        ))
    }
    merge.tree[[1]] <- .fill_tree(merge.tree[[1]], batches)
    merge.tree[[2]] <- .fill_tree(merge.tree[[2]], batches)
    merge.tree
}

.get_next_merge <- function(merge.tree, path=NULL) {
    if (!is.list(merge.tree[[1]]) && !is.list(merge.tree[[2]])) {
        list(left=merge.tree[[1]], right=merge.tree[[2]], path=path)
    } else if (is.list(merge.tree[[2]])) {
        .get_next_merge(merge.tree[[2]], c(path, 2L))
    } else {
        .get_next_merge(merge.tree[[1]], c(path, 1L))
    }
}

.update_tree <- function(merge.tree, path, ...) {
    if (length(path)) {
        return(new("MNN_treenode", ...))
    }
    merge.tree[[path[1]]] <- .update_tree(merge.tree[[path[[1]]]], path[-1], ...)
    merge.tree
}

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
