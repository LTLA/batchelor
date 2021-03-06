#' Take the intersection of rows across batches
#'
#' Subset multiple batches so that they have the same number and order of rows.
#'
#' @param ... One or more matrix-like objects containing single-cell gene expression matrices.
#' Alternatively, one or more \linkS4class{SingleCellExperiment} objects can be supplied.
#' @param subset.row A vector specifying the subset of genes to retain.
#' Defaults to \code{NULL}, in which case all genes are retained.
#' @param keep.all Logical scalar indicating whether the data should actually be subsetted to \code{subset.row}.
#' If \code{FALSE}, \code{subset.row} is checked but not used.
#'
#' @return A list containing the row-subsetted contents of \code{...}, each of which have the same number and order of rows.
#'
#' @author Aaron Lun
#'
#' @details
#' If any entry of \code{...} contains no row names or if the intersection is empty, an error is raised.
#'
#' If all entries of \code{...} already have the same row names in the same order, this function is a no-op.
#' In such cases, \code{subset.row} can be a character, integer or logical vector.
#'
#' If entries of \code{...} do not have identical row names, \code{subset.row} can only be a character vector.
#' Any other type will lead to an error as the interpretation of the subset is not well defined.
#' 
#' Setting \code{keep.all=TRUE} option gives the same result as having \code{subset.row=NULL} in the first place.
#' However, it is still useful for checking that \code{subset.row} is a character vector.
#' This ensures that downstream applications can safely use \code{subset.row} in conjunction with, e.g., \code{correct.all=TRUE} for \code{\link{fastMNN}}.
#' 
#' @examples
#' X <- rbind(A=c(1,2), B=c(3,4))
#' Y <- rbind(a=c(1,2), B=c(3,4))
#' intersectRows(X, Y) # Only B is retained.
#'
#' # Error is raised when no genes are retained:
#' X <- rbind(A=c(1,2), B=c(3,4))
#' Y <- rbind(a=c(1,2), b=c(3,4))
#' try(intersectRows(X, Y))
#'
#' # Error is raised for non-character subset.row 
#' # when row names are not identical:
#' X <- rbind(A=c(1,2), B=c(3,4), C=c(5,6))
#' Y <- rbind(a=c(1,2), B=c(3,4), C=c(5,6))
#' intersectRows(X, Y)
#' intersectRows(X, Y, subset.row="B")
#' try(intersectRows(X, Y, subset.row=1))
#'
#' # Setting keep.all=TRUE only checks but does not apply subset.row.
#' intersectRows(X, Y, subset.row="B", keep.all=TRUE)
#' try(intersectRows(X, Y, subset.row=1, keep.all=TRUE))
#'
#' @export
#' @importFrom scuttle .unpackLists
intersectRows <- function(..., subset.row=NULL, keep.all=FALSE) {
    all.batches <- .unpackLists(...)

    all.genes <- lapply(all.batches, rownames)
    universe <- Reduce(intersect, all.genes)
    if (length(universe)==0) {
        stop("no genes remaining in the intersection")
    }

    non.same <- FALSE 
    for (i in seq_along(all.batches)) {
        if (!identical(universe, rownames(all.batches[[i]]))) {
            non.same <- TRUE
            all.batches[[i]] <- all.batches[[i]][universe,,drop=FALSE]
        }
    }

    if (!is.null(subset.row)) {
        if (non.same && !is.character(subset.row)) {
            stop("only character 'subset.row' is allowed when '...' have different rownames")
        }
        if (!keep.all) {
            all.batches <- lapply(all.batches, function(x) x[subset.row,,drop=FALSE])
        }
    }

    all.batches
}


