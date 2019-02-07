#' Cosine normalization
#'
#' Perform cosine normalization on the column vectors of an expression matrix.
#'    
#' @param x A gene expression matrix with cells as columns and genes as rows.
#' @param mode A string specifying the output to be returned.
#'
#' @details
#' Cosine normalization removes scaling differences between expression vectors.
#' In the context of batch correction, this is usually applied to remove differences between batches that are normalized separately.
#' For example, \code{\link{fastMNN}} uses this function on the log-expression vectors by default.
#' 
#' Technically, separate normalization introduces scaling differences in the normalized expression, which should manifest as a shift in the log-transformed expression.
#' However, in practice, single-cell data will contain many small counts (where the log function is near-linear) or many zeroes (which remain zero when the pseudo-count is 1).
#' In these applications, scaling differences due to separate normalization are better represented as scaling differences in the log-transformed values.
#' 
#' If applied to the raw count vectors, cosine normalization is similar to library size-related (i.e., L1) normalization.
#' However, we recommend using dedicated methods for computing size factors to normalize raw count data.
#' 
#' While the default is to directly return the cosine-normalized matrix, it may occasionally be desirable to obtain the L2 norm, 
#' e.g., to apply an equivalent normalization to other matrices.
#' This can be achieved by setting \code{mode} accordingly.
#'
#' The function will return a \linkS4class{DelayedMatrix} if \code{x} is a \linkS4class{DelayedMatrix}.
#' This aims to delay the calculation of cosine-normalized values for very large matrices.
#'
#' @return
#' If \code{mode="matrix"}, a double-precision matrix of the same dimensions as \code{X} is returned, containing cosine-normalized values.
#' 
#' If \code{mode="l2norm"}, a double-precision vector is returned containing the L2 norm for each cell.
#' 
#' If \code{mode="all"}, a named list is returned containing the fields \code{"matrix"} and \code{"l2norm"}, which are as described above.
#'
#' @author
#' Aaron Lun
#' 
#' @seealso
#' \code{\link{mnnCorrect}} and \code{\link{fastMNN}}, where this function gets used.
#'
#' @examples
#' A <- matrix(rnorm(1000), nrow=10)
#' str(cosineNorm(A))
#' str(cosineNorm(A, mode="l2norm"))
#' 
#' @export
#' @importClassesFrom DelayedArray DelayedMatrix
#' @importMethodsFrom DelayedArray sweep
#' @importFrom methods is
cosineNorm <- function(x, mode=c("matrix", "all", "l2norm"))
# Computes the cosine norm, with some protection from zero-length norms.
#
# written by Aaron Lun
# 5 July 2018
{
    mode <- match.arg(mode)
    if (is(x, "DelayedMatrix")) {
        l2 <- .Call(cxx_cosine_norm, x, FALSE)[[2]]
        out <- list(matrix=sweep(x, 2, l2, "/", check.margin=FALSE), l2norm=l2)
    } else {
        out <- .Call(cxx_cosine_norm, x, mode!="l2norm")
        names(out) <- c("matrix", "l2norm")
    }

    switch(mode, all=out, matrix=out$matrix, l2norm=out$l2norm)
}
