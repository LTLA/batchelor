#' Correct SingleCellExperiment objects
#'
#' Apply a correction to multiple \linkS4class{SingleCellExperiment} objects,
#' while also combining the assay data and column metadata for easy use.
#'
#' @param ... One or more \linkS4class{SingleCellExperiment} objects.
#' @param assay.type A string or integer scalar specifying the assay to use for correction.
#' @inheritParams batchCorrect
#' @param combine.assays Character vector specifying the assays from each entry of \code{...} to combine together without correction.
#' By default, any named assay that is present in all entries of \code{...} are combined.
#' This can be set to \code{character(0)} to avoid combining any assays.
#' @param combine.coldata Character vector specifying the column metadata fields from each entry of \code{...} to combine together.
#' By default, any column metadata field that is present in all entries of \code{...} is combined.
#' This can be set to \code{character(0)} to avoid combining any metadata.
#'
#' @return 
#' A SingleCellExperiment containing the merged expression values in the first assay
#' and a \code{batch} column metadata field specifying the batch of origin for each cell,
#' as described in \code{\link{batchCorrect}}.
#'
#' Additional assays may be present depending on \code{combine.assays}.
#' These contain uncorrected values from each batch that have been simply \code{cbind}ed together.
#' Additional column metadata fields may also be present depending on \code{combine.coldata}.
#'
#' @details
#' This function makes it easy to retain information from the original SingleCellExperiment objects in the post-merge object.
#' Operations like differential expression analyses can be easily performed on the uncorrected expression values,
#' while common annotation can be leveraged in cell-based analyses like clustering.
#'
#' If \code{combine.assays} contains a field that overlaps with the name of the corrected assay from \code{\link{batchCorrect}},
#' a warning will be raised and the corrected assay will be preferentially retained.
#' Similarly, any field named \code{"batch"} in \code{combine.coldata} will be ignored.
#'
#' @author Aaron Lun
#' @examples
#' sce1 <- scater::mockSCE()
#' sce1 <- scater::logNormCounts(sce1)
#' sce2 <- scater::mockSCE()
#' sce2 <- scater::logNormCounts(sce2)
#'
#' f.out <- correctExperiments(sce1, sce2)
#' colData(f.out)
#' assayNames(f.out)
#' 
#' @seealso
#' \code{\link{batchCorrect}}, which does the correction inside this function.
#'
#' \code{\link{noCorrect}}, used to combine uncorrected values for the other assays.
#' @export
#' @importFrom SummarizedExperiment assayNames assay<- assay colData<- colData
correctExperiments <- function(..., batch=NULL, restrict=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", 
    PARAM=FastMnnParam(), combine.assays=NULL, combine.coldata=NULL) 
{
    x <- list(...)
    args <- c(x, list(subset.row=subset.row, correct.all=correct.all, batch=batch))
    merged <- do.call(batchCorrect, c(args, list(restrict=restrict, assay.type=assay.type, PARAM=PARAM)))

    # Adding additional assays.
    if (is.null(combine.assays)) {
        combine.assays <- Reduce(intersect, lapply(x, assayNames))
    }
    merged.assays <- assayNames(merged)
    if (any(combine.assays %in% merged.assays)) {
        warning("ignoring assays with same name as 'batchCorrect' output")
        combine.assays <- setdiff(combine.assays, merged.assays)
    }
    for (nm in combine.assays) {
        nocorrect <- do.call(batchCorrect, c(args, list(assay.type=nm, PARAM=NoCorrectParam())))
        assay(merged, nm) <- assay(nocorrect)
    }

    # Adding additional colData.
    if (is.null(combine.coldata)) {
        combine.coldata <- Reduce(intersect, lapply(x, FUN=function(y) colnames(colData(y))))
    }
    merged.coldata <- colnames(colData(merged))
    if (any(combine.coldata %in% merged.coldata)) {
          warning("ignoring 'colData' fields overlapping 'batchCorrect' output")
          combine.coldata <- setdiff(combine.coldata, merged.coldata)
    }

    if (length(combine.coldata)) {
        collected <- lapply(x, FUN=function(y) colData(y)[,combine.coldata,drop=FALSE])
        colData(merged) <- cbind(colData(merged), do.call(rbind, collected))
    }

    merged
}
