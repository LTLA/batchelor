.restore_original_order <- function(batch.ordering, ncells.per.batch) 
# Generates a permutation vector to recover the old order of cells
# (assuming that within-batch order is preserved).
{ 
    nbatches <- length(batch.ordering)
    if (nbatches!=length(ncells.per.batch)) {
        stop("length of batch information vectors are not equal")
    }
    if (nbatches==0L) {
        return(integer(0))
    }

    reorder <- vector("list", nbatches)
    last <- 0L
    for (idx in batch.ordering) {
        ncells <- ncells.per.batch[idx]
        reorder[[idx]] <- last + seq_len(ncells)
        last <- last + ncells
    }
    unlist(reorder)
}

.reindex_pairings <- function(pairings, new.order) 
# Restores the MNN pairing indices to the original positions,
# i.e., after cells have been reordered by 'new.order'.
{
    rev.order <- integer(length(new.order))
    rev.order[new.order] <- seq_along(new.order)
    for (x in seq_along(pairings)) {
        current <- pairings[[x]]
        current$left <- rev.order[current$left]
        current$right <- rev.order[current$right]
        pairings[[x]] <- current
    }
    pairings
}

#' @importFrom S4Vectors metadata<- metadata List
.add_batch_names <- function(output, batches) {
    if (!is.null(names(batches))) {
        if (anyDuplicated(names(batches))) {
            stop("names of batches should be unique")
        }
        output$batch <- names(batches)[output$batch]
        L1 <- metadata(output)$merge.info$left 
        R1 <- metadata(output)$merge.info$right 
        L2 <- R2 <- List()
        for (i in seq_along(L1)) {
            L2[[i]] <- names(batches)[L1[[i]]]
            R2[[i]] <- names(batches)[R1[[i]]]
        }
        metadata(output)$merge.info$left <- L2
        metadata(output)$merge.info$right <- R2
        colnames(metadata(output)$merge.info$lost.var) <- names(batches)
    } 
    output
}
