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
setClass("NoCorrectParam", contains="BatchelorParam")

#####################################################
# Internal classes for MNN ordering.

setClass("MNN_order", contains="VIRTUAL", 
    slots=c(batches="ANY", restrict="ANY", # storing basic input information.
        reference="ANY", reference.restrict="ANY", reference.indices="integer", # storing reference information.
        current.index="integer", mnn.result="list" # storing values to be returned.
    )
)

setClass("MNN_supplied_order", contains="MNN_order", slots=c(ordering="integer"))

setClass("MNN_auto_order", contains="MNN_order", slots=c(precomputed="list", restricted.batches="ANY"))
