.restore_original_order <- function(batch.ordering, ncells.per.batch) 
# Generates a permutation vector to recover the old order of cells
# (assuming that within-batch order is preserved).
{ 
    reorder <- vector("list", length(batch.ordering))
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
        current$first <- rev.order[current$first]
        current$second <- rev.order[current$second]
        pairings[[x]] <- current
    }
    pairings
}
