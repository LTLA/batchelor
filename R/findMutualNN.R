#' Find mutual nearest neighbors
#'
#' Find mutual nearest neighbors (MNN) across two data sets.
#'
#' @param data1 A numeric matrix containing samples (e.g., cells) in the rows and variables/dimensions in the columns.
#' @param data2 A numeric matrix like \code{data1} for another data set with the same variables/dimensions.
#' @param k1 Integer scalar specifying the number of neighbors to search for in \code{data1}.
#' @param k2 Integer scalar specifying the number of neighbors to search for in \code{data2}.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying the neighbour search algorithm to use.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#'
#' @return
#' A list containing the integer vectors \code{first} and \code{second}.
#' Corresponding entries in \code{first} and \code{second} specify a MNN pair of cells from \code{data1} and \code{data2}, respectively.
#'
#' @details
#' The concept of a MNN pair can be explained by considering cells in each of two data sets. 
#' For each cell in data set 1, the set of \code{k2} nearest cells in data set 2 is identified, based on the Euclidean distance in expression space.
#' For each cell in data set 2, the set of \code{k1} nearest cells in data set 1 is similarly identified.
#' Two cells in different batches are considered to be MNNs if each cell is in the other's set.
#'
#' @author
#' Aaron Lun
#'
#' @seealso
#' \code{\link{queryKNN}} for the underlying neighbor search code.
#'
#' \code{\link{fastMNN}} and \code{\link{mnnCorrect}}, which call this function to identify MNNs.
#'
#' @examples
#' B1 <- matrix(rnorm(10000), ncol=50) # Batch 1 
#' B2 <- matrix(rnorm(10000), ncol=50) # Batch 2
#' out <- findMutualNN(B1, B2, k1=20)
#' head(out$first)
#' head(out$second)
#'
#' @export
#' @importFrom BiocNeighbors queryKNN KmknnParam
#' @importFrom BiocParallel SerialParam
findMutualNN <- function(data1, data2, k1, k2=k1, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) 
{
    data1 <- as.matrix(data1)
    data2 <- as.matrix(data2)
    W21 <- queryKNN(data2, query=data1, k=k2, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)
    W12 <- queryKNN(data1, query=data2, k=k1, BNPARAM=BNPARAM, BPPARAM=BPPARAM, get.distance=FALSE)
    out <- find_mutual_nns(W21$index, W12$index)
    names(out) <- c("first", "second")
    return(out)
}
