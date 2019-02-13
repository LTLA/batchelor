#' @importMethodsFrom BiocGenerics nrow ncol
#' @importFrom BiocGenerics colnames rownames
.rename_output <- function(output, batches, subset.row=NULL, cells.in.columns=TRUE)
# Adds dimension names to the output according to the input 'batches'.
# Also replaces NULL names with empty strings so that the calling function doesn't have
# to worry about inputs where some batches are named and others are not.
{
    GENERATE_NAMES <- function(batches, OTHERDIMFUN, OTHERDIMNAMEFUN) {
        collected <- lapply(batches, OTHERDIMNAMEFUN)
        nulled <- vapply(collected, is.null, FUN.VALUE=TRUE)
        if (any(nulled) && !all(nulled)) {
            collected[nulled] <- lapply(batches[nulled], FUN=function(x) character(OTHERDIMFUN(x)))
        }
        unlist(collected)
    }

    if (!cells.in.columns) {
        cell.names <- GENERATE_NAMES(batches, nrow, rownames)
        rownames(output) <- cell.names
        colnames(output) <- colnames(batches[[1]])
    } else {
        cell.names <- GENERATE_NAMES(batches, ncol, colnames)
        colnames(output) <- cell.names

        feat.names <- rownames(batches[[1]])
        if (!is.null(feat.names) && !is.null(subset.row)) {
            feat.names <- feat.names[.row_subset_to_index(batches[[1]], subset.row)]
        }
        rownames(output) <- feat.names
    }

    output
}

.create_batch_names <- function(batch.labels, ncells.per.batch) 
# Creates batch names if they aren't already available.
{
    if (is.null(batch.labels)) {
        batch.labels <- seq_along(ncells.per.batch)
    }
    batch.ids <- rep(batch.labels, ncells.per.batch)
    list(labels=batch.labels, ids=batch.ids)
}
