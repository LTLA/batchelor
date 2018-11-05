#' @importMethodsFrom BiocGenerics nrow ncol
#' @importFrom BiocGenerics colnames rownames
.check_batch_consistency <- function(batches, byrow=TRUE, ignore.null=FALSE) 
# Checking for identical number of rows (and rownames).
# It can also do the same for columns, if we're dealing with PC results.
{
    if (length(batches) < 2L) {
        return(NULL)
    }

    if (byrow) {
        DIMFUN <- nrow
        DIMNAMEFUN <- rownames
        DIM <- "row"
    } else {
        DIMFUN <- ncol
        DIMNAMEFUN <- colnames
        DIM <- "column"
    }

    first <- batches[[1]]
    ref.n <- DIMFUN(first)
    ref.names <- DIMNAMEFUN(first)

    for (b in 2:length(batches)) { 
        current <- batches[[b]]
        if (!identical(DIMFUN(current), ref.n)) {
            stop(sprintf("number of %ss is not the same across batches", DIM))
        }

        cur.names <- DIMNAMEFUN(current)
        if (ignore.null) { 
            if (is.null(cur.names)) { 
                cur.names <- ref.names
            } else if (is.null(ref.names)) {
                ref.names <- cur.names
            }
        }
        if (!identical(cur.names, ref.names)) {
            stop(sprintf("%s names are not the same across batches", DIM))
        }
    }

    return(NULL)
}

#' @importFrom SingleCellExperiment isSpike spikeNames
.check_spike_consistency <- function(batches) 
# Checking for identical spike-in sets and (overall) identities.
{
    if (length(batches) < 2L) {
        return(NULL)
    }

    ref.spike.names <- spikeNames(batches[[1]])
    ref.spike <- isSpike(batches[[1]])
    for (b in seq_along(batches)) {
        if (!identical(ref.spike.names, spikeNames(batches[[b]]))) {
            stop("spike-in sets differ across batches")
        }
        if (!identical(ref.spike, isSpike(batches[[b]]))) {
            stop("spike-in identities differ across batches")
        }
    }
    return(NULL)
}

#' @importFrom methods is
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.check_if_SCEs <- function(batches) 
# Checks that everyone is either an SCE or is not.
# Returns FALSE if 'batches' is empty.
{
    all.sce <- vapply(batches, is, class2="SingleCellExperiment", FUN.VALUE=TRUE)
    if (length(unique(all.sce)) > 1L) {
        stop("cannot mix SingleCellExperiments and other objects")
    }
    any(all.sce) # don't do all.sce[1], avoid errors when length(batches)==0L.
}

.divide_into_batches <- function(x, batch) 
# Splits 'x' by column in 'batch'.
{
    batch <- as.factor(batch)
    if (length(batch)!=ncol(x)) {
        stop("'length(batch)' and 'ncol(x)' are not the same")
    }

    output <- vector("list", nlevels(batch))
    names(output) <- levels(batch)
    reorder <- integer(ncol(x))
    last <- 0L

    for (b in levels(batch)) {
        keep <- batch==b
        current <- x[,keep,drop=FALSE]
        output[[b]] <- current
        reorder[keep] <- last + seq_len(ncol(current))
        last <- last + ncol(current)
    }

    list(batches=output, reorder=reorder) 
}
