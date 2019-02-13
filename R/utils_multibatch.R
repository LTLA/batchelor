.create_batch_names <- function(batch.labels, ncells.per.batch) 
# Creates batch names if they aren't already available.
{
    if (is.null(batch.labels)) {
        batch.labels <- seq_along(ncells.per.batch)
    }
    batch.ids <- rep(batch.labels, ncells.per.batch)
    list(labels=batch.labels, ids=batch.ids)
}
