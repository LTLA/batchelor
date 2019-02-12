#' Divide into batches
#'
#' Divide a single input object into multiple separate objects according to their batch of origin.
#'
#' @param x A matrix-like object where one dimension corresponds to cells and another represents features.
#' @param batch A factor specifying the batch to which each cell belongs.
#' @param byrow A logical scalar indicating whether rows correspond to cells.
#' @param restrict A subsetting vector specifying which cells should be used for correction.
#'
#' @details
#' This function is intended for internal use and other package developers.
#' It splits a single input object into multiple batches, allowing developers to use the same code for the scenario where \code{batch} is supplied with a single input.
#' 
#' @return
#' A list containing:
#' \itemize{
#' \item \code{batches}, a named list of matrix-like objects where each element corresponds to a level of \code{batch} and contains all cells from that batch.
#' \item \code{reorder}, an integer vector to be applied to the combined \code{batches} to recover the ordering of cells in \code{x}. 
#' \item \code{restricted}, a named list of integer vectors specifying which cells are to be used for correction.
#' Set to \code{NULL} if the input \code{restrict} was also \code{NULL}.
#' }
#'
#' @author Aaron Lun
#'
#' @examples
#' X <- matrix(rnorm(1000), ncol=100)
#' out <- divideIntoBatches(X, sample(3, 100, replace=TRUE))
#' names(out)
#' 
#' # Recovering original order.
#' Y <- do.call(cbind, out$batches)
#' Z <- Y[,out$reorder]
#' all.equal(Z, X) # should be TRUE.
#'
#' @export
divideIntoBatches <- function(x, batch, byrow=FALSE, restrict=NULL) 
{
    if (is.null(batch)) { 
        stop("'batch' must be specified if '...' has only one object")
    }

    batch <- as.factor(batch)
    if (byrow) {
        if (length(batch)!=nrow(x)) {
            stop("'length(batch)' and 'nrow(x)' are not the same")
        }
    } else {
        if (length(batch)!=ncol(x)) {
            stop("'length(batch)' and 'ncol(x)' are not the same")
        }
    }

    output <- vector("list", nlevels(batch))
    names(output) <- levels(batch)
    reorder <- integer(ncol(x))
    last <- 0L
    
    if (!is.null(restrict)) {
        if (byrow) {
            tmp <- .row_subset_to_index(x, restrict)
            restrict <- logical(nrow(x))
        } else {
            tmp <- .col_subset_to_index(x, restrict)
            restrict <- logical(ncol(x))
        }
        restrict[tmp] <- TRUE
        restricted <- output
    } else {
        restricted <- NULL
    }

    for (b in levels(batch)) {
        keep <- batch==b

        if (byrow) {
            current <- x[keep,,drop=FALSE]
            N <- nrow(current)
        } else {
            current <- x[,keep,drop=FALSE]
            N <- ncol(current)
        }

        if (!is.null(restrict)) {
            restricted[[b]] <- which(restrict[keep])
        }

        output[[b]] <- current
        reorder[keep] <- last + seq_len(N)
        last <- last + N
    }

    list(batches=output, reorder=reorder, restricted=restricted) 
}
