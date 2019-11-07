#' Correct SingleCellExperiment objects
#'
#' Apply a correction to multiple \linkS4class{SingleCellExperiment} objects,
#' while also combining the assay data and column metadata for easy use.
#'
#' @param ... One or more \linkS4class{SingleCellExperiment} objects.
#' If multiple objects are supplied, each object is assumed to contain all and only cells from a single batch.
#' If a single object is supplied, \code{batch} should also be specified.
#' @param assay.type A string or integer scalar specifying the assay to use for correction.
#' @inheritParams batchCorrect
#' @param combine.assays Character vector specifying the assays from each entry of \code{...} to combine together without correction.
#' By default, any named assay that is present in all entries of \code{...} are combined.
#' This can be set to \code{character(0)} to avoid combining any assays.
#' @param combine.coldata Character vector specifying the column metadata fields from each entry of \code{...} to combine together.
#' By default, any column metadata field that is present in all entries of \code{...} is combined.
#' This can be set to \code{character(0)} to avoid combining any metadata.
#' @param include.rowdata Logical scalar indicating whether the function should attempt to include \code{\link{rowRanges}}.
#'
#' @return 
#' A SingleCellExperiment containing the merged expression values in the first assay
#' and a \code{batch} column metadata field specifying the batch of origin for each cell,
#' as described in \code{\link{batchCorrect}}.
#'
#' @details
#' This function makes it easy to retain information from the original SingleCellExperiment objects in the post-merge object.
#' Operations like differential expression analyses can be easily performed on the uncorrected expression values,
#' while common annotation can be leveraged in cell-based analyses like clustering.
#'
#' Additional assays may be added to the merged object, depending on \code{combine.assays}.
#' This will usually contain uncorrected values from each batch that have been simply \code{cbind}ed together.
#' If \code{combine.assays} contains a field that overlaps with the name of the corrected assay from \code{\link{batchCorrect}},
#' a warning will be raised and the corrected assay will be preferentially retained.
#'
#' Any column metadata fields that are shared will also be included in the merged object by default (tunable by setting \code{combine.coldata}).
#' If any existing field is named \code{"batch"}, it will be ignored in favor of that produced by \code{\link{batchCorrect}} and a warning is emitted.
#'
#' Row metadata is only included in the merged object if \code{include.rowdata=TRUE} \emph{and} all row metadata objects are identical across objects in \code{...}.
#' If not, a warning is emitted and no row metadata is attached to the merged object.
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
#' 
#' @export
#' @importFrom SummarizedExperiment assayNames assay<- assay colData<- colData 
#' rowRanges rowRanges<- rowData rowData<-
#' @importFrom S4Vectors mcols mcols<-
correctExperiments <- function(..., batch=NULL, restrict=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", 
    PARAM=FastMnnParam(), combine.assays=NULL, combine.coldata=NULL, include.rowdata=TRUE) 
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

    # Adding additional rowRanges only if they are identical across objects.
    if (include.rowdata) {
        all.rd <- lapply(x, FUN=rowRanges)
        is.id <- TRUE
        for (i in seq_len(length(all.rd)-1L)) {
            if (is.id <- identical(all.rd[[1]], all.rd[[i+1]])) {
                break
            }
        }

        if (!is.id) {
            warning("ignoring non-identical 'rowRanges' fields")
        } else {
            merged.rd <- rowData(merged) 
            replacement <- all.rd[[1]]
            if (!correct.all && !is.null(subset.row)) {
                replacement <- replacement[subset.row,]
            }

            combined <- cbind(mcols(replacement), merged.rd)
            skip <- duplicated(colnames(combined), fromLast=TRUE)
            if (any(skip)) {
                warning("ignoring 'rowData' fields overlapping 'batchCorrect' output")
            }

            mcols(replacement) <- combined[,!skip,drop=FALSE]
            rowRanges(merged) <- replacement
        }
    }

    merged
}
