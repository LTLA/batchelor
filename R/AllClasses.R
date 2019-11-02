#' @export
#' @importClassesFrom S4Vectors SimpleList
#' @import methods
#' @rdname BatchelorParam
setClass("BatchelorParam", contains="SimpleList")

#' @export
#' @rdname BatchelorParam
setClass("ClassicMnnParam", contains="BatchelorParam")

#' @export
#' @rdname BatchelorParam
setClass("FastMnnParam", contains="BatchelorParam")

#' @export
#' @rdname BatchelorParam
setClass("RescaleParam", contains="BatchelorParam")

#' @export
#' @rdname BatchelorParam
setClass("RegressParam", contains="BatchelorParam")

#' @export
#' @rdname BatchelorParam
setClass("NoCorrectParam", contains="BatchelorParam")
