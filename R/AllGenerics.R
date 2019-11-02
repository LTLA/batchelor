#' @export
#' @rdname batchCorrect
#' @import methods
setGeneric("batchCorrect", function(..., batch=NULL, restrict=NULL, subset.row=NULL, correct.all=FALSE, assay.type=NULL, get.spikes=FALSE, PARAM) 
    standardGeneric("batchCorrect"), signature=c("PARAM"))

# Note that assay.type default in each method will override that in the generic,
# so there's no need to reset it explicitly in each method (see ?setMethod). 
# We don't set it here as different methods may expect different assay types, 
# e.g., a correction that only operates on counts rather than log-counts.
