#' @export
#' @rdname batchCorrect
#' @import methods
setGeneric("batchCorrect", function(..., batch=NULL, subset.row=NULL, correct.all=FALSE, assay.type=NULL, get.spikes=FALSE, PARAM) 
    standardGeneric("batchCorrect"), signature=c("PARAM"))

#' @importFrom BiocNeighbors queryKNN buildIndex KmknnParam
#' @importFrom BiocParallel SerialParam
setGeneric(".advance", function(x, k, BNPARAM=KmknnParam(), BPPARAM=SerialParam()) standardGeneric(".advance"))
