#' BatchelorParam methods
#'
#' Constructors and methods for the batchelor parameter classes.
#'
#' @param ... Named arguments to pass to individual methods upon dispatch.
#' These should not include arguments named in the \code{\link{batchCorrect}} generic.
#'
#' @details
#' BatchelorParam objects are intended to store method-specific parameter settings to pass to the \code{\link{batchCorrect}} generic.
#' These values should refer to data-agnostic parameters; parameters that depend on data (or the data itself) should be specified directly in the \code{\link{batchCorrect}} call.
#'
#' The BatchelorParam classes are all derived from \linkS4class{SimpleList} objects and have the same available methods, e.g., \code{[[}, \code{$}.
#' These can be used to access or modify the object after construction.
#'
#' Note that the BatchelorParam class itself is not useful and should not be constructed directly.
#' Instead, users should use the constructors shown above to create instances of the desired subclass.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{batchCorrect}}, where the BatchelorParam objects are used for dispatch to individual methods.
#'
#' @return
#' The constructors will return a BatchelorParam object of the specified subclass, containing parameter settings for the corresponding batch correction method.
#'
#' @examples
#' # Specifying the number of neighbors, dimensionality.
#' fp <- FastMnnParam(k=20, d=10)
#' fp
#'
#' # List-like behaviour:
#' fp$k
#' fp$k <- 10
#' fp$k
#'
#' @docType methods
#'
#' @rdname BatchelorParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
ClassicMnnParam <- function(...) {
    new("ClassicMnnParam", SimpleList(list(...)))
}

#' @rdname BatchelorParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
FastMnnParam <- function(...) {
    new("FastMnnParam", SimpleList(list(...)))
}

#' @rdname BatchelorParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
RescaleParam <- function(...) {
    new("RescaleParam", SimpleList(list(...)))
}

#' @rdname BatchelorParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
RegressParam <- function(...) {
    new("RegressParam", SimpleList(list(...)))
}

#' @rdname BatchelorParam
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom methods new
NoCorrectParam <- function(...) {
    new("NoCorrectParam", SimpleList(list(...)))
}
