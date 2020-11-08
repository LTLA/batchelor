#' Correct SingleCellExperiment objects
#'
#' Apply a correction to multiple \linkS4class{SingleCellExperiment} objects,
#' while also combining the assay data and column metadata for easy downstream use.
#' This augments the simpler \code{\link{batchCorrect}} function, which returns only the corrected values.
#'
#' @param ... One or more \linkS4class{SingleCellExperiment} objects.
#' If multiple objects are supplied, each object is assumed to contain all and only cells from a single batch.
#' If a single object is supplied, \code{batch} should also be specified.
#'
#' Alternatively, one or more lists of SingleCellExperiments can be provided;
#' this is flattened so that it is as if the objects inside were passed directly to \code{...}.
#' @param assay.type A string or integer scalar specifying the assay to use for correction.
#' @inheritParams batchCorrect
#' @param combine.assays Character vector specifying the assays from each entry of \code{...} to combine together without correction.
#' By default, any named assay that is present in all entries of \code{...} is combined.
#' This can be set to \code{character(0)} to avoid combining any assays.
#' @param combine.coldata Character vector specifying the column metadata fields from each entry of \code{...} to combine together.
#' By default, any column metadata field that is present in all entries of \code{...} is combined.
#' This can be set to \code{character(0)} to avoid combining any metadata.
#' @param include.rowdata Logical scalar indicating whether the function should attempt to include \code{\link{rowRanges}}.
#' @param add.single Logical scalar indicating whether merged fields should be added to the original SingleCellExperiment.
#' Only relevant when a single object is provided in \code{...}.
#' If \code{TRUE}, \code{combine.assays}, \code{combine.coldata} and \code{include.rowdata} are ignored.
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
#' \itemize{
#' \item All assays shared across the original objects are \code{cbind}ed and added to the merged object.
#' This can be controlled with \code{combine.assays}.
#' Any original assay with the same name as an assay in the output of \code{\link{batchCorrect}} will be ignored with a warning.
#' \item Any column metadata fields that are shared will also be included in the merged object.
#' This can be controlled with \code{combine.coldata}.
#' If any existing field has the same name as any \code{\link{colData}} field produced by \code{\link{batchCorrect}},
#' it will be ignored in favor of the latter.
#' \item Row metadata from \code{...} is included in the merged object if \code{include.rowdata=TRUE}.
#' In such cases, only non-conflicting row data fields are preserved,
#' i.e., fields with different names or identically named fields with the same values between objects in \code{...}.
#' Any conflicting fields are ignored with a warning.
#' \code{rowRanges} are only preserved if they are identical (ignoring the \code{\link{mcols}}) for all objects in \code{...}.
#' }
#'
#' If a single SingleCellExperiment object was supplied in \code{...}, the default behavior is to prepend all \code{\link{assays}}, \code{\link{reducedDims}}, \code{\link{colData}}, \code{\link{rowData}} and \code{\link{metadata}} fields from the merged object into the original (removing any original entries with names that overlap those of the merged object).
#' This is useful as it preserves all (non-overlapping) aspects of the original object, especially the reduced dimensions that cannot, in general, be sensibly combined across multiple objects.
#' Setting \code{add.single=FALSE} will force the creation of a new SingleCellExperiment rather than prepending.
#'
#' @author Aaron Lun
#' @examples
#' sce1 <- scuttle::mockSCE()
#' sce1 <- scuttle::logNormCounts(sce1)
#' sce2 <- scuttle::mockSCE()
#' sce2 <- scuttle::logNormCounts(sce2)
#'
#' f.out <- correctExperiments(sce1, sce2)
#' colData(f.out)
#' assayNames(f.out)
#' 
#' @seealso
#' \code{\link{batchCorrect}}, which does the correction inside this function.
#'
#' \code{\link{noCorrect}}, for another method to combine uncorrected assay values. 
#' 
#' @export
#' @importFrom scuttle .unpackLists
correctExperiments <- function(..., batch=NULL, restrict=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", 
    PARAM=FastMnnParam(), combine.assays=NULL, combine.coldata=NULL, include.rowdata=TRUE, add.single=TRUE) 
{
    x <- .unpackLists(...)
    merged <- batchCorrect(x, subset.row=subset.row, correct.all=correct.all, batch=batch,
        restrict=restrict, assay.type=assay.type, PARAM=PARAM)

    if (length(x)==1L && add.single) {
        .add.single_sce(x[[1]], merged, subset.row=subset.row, correct.all=correct.all)
    } else {
        .create_fresh_combined_sce(x, merged, subset.row=subset.row, correct.all=correct.all,
            combine.assays=combine.assays, combine.coldata=combine.coldata, include.rowdata=include.rowdata)
    }
}

#' @importFrom SummarizedExperiment assayNames assay<- assay 
#' colData<- colData rowRanges rowRanges<- rowData rowData<-
.create_fresh_combined_sce <- function(x, merged, subset.row=NULL, correct.all=FALSE,
    combine.assays=NULL, combine.coldata=NULL, include.rowdata=TRUE) 
{
    # Adding additional assays.
    if (is.null(combine.assays)) {
        combine.assays <- Reduce(intersect, lapply(x, assayNames))
    }
    combine.assays <- .eliminate_overlaps(assayNames(merged), combine.assays, msg="'assays'")

    for (nm in combine.assays) {
        raw.ass <- lapply(x, assay, i=nm)
        if (!is.null(subset.row) && !correct.all) {
            raw.ass <- lapply(raw.ass, "[", i=subset.row, , drop=FALSE)
        }
        assay(merged, nm) <- do.call(cbind, raw.ass)
    }

    # Adding additional colData.
    if (is.null(combine.coldata)) {
        combine.coldata <- Reduce(intersect, lapply(x, FUN=function(y) colnames(colData(y))))
    }
    combine.coldata <- .eliminate_overlaps(colnames(colData(merged)), combine.coldata, msg="'colData' fields")

    if (length(combine.coldata)) {
        collected <- lapply(x, FUN=function(y) colData(y)[,combine.coldata,drop=FALSE])
        colData(merged) <- cbind(colData(merged), do.call(rbind, collected))
    }

    # Adding additional rowRanges.
    if (include.rowdata) {
        all.rrw <- lapply(x, FUN=rowRanges)
        combined <- .accumulate_rowdata(all.rrw)
        combined.rd <- combined$df
        combined.rr <- combined$ranges

        if (!correct.all && !is.null(subset.row)) {
            combined.rd <- combined.rd[subset.row,,drop=FALSE]
            combined.rr <- combined.rr[subset.row]
        }

        merged.rd <- rowData(merged) 
        leftovers.rd <- .eliminate_overlaps(colnames(merged.rd), colnames(combined.rd), msg="'rowData' fields")
        combined.rd <- combined.rd[,leftovers.rd,drop=FALSE]

        # Setting the rowRanges first as this wipes out any existing rowData.
        if (!is.null(combined.rr)) {
            rowRanges(merged) <- combined.rr
        }
        rowData(merged) <- cbind(merged.rd, combined.rd)
    }

    merged
}

.eliminate_overlaps <- function(priority, other, msg="stuff") {
    if (any(other %in% priority)) {
        warning(sprintf("ignoring %s with same name as 'batchCorrect' output", msg))
        other <- setdiff(other, priority)
    }
    other
}

#' @importFrom S4Vectors mcols mcols<- DataFrame make_zero_col_DFrame I
.accumulate_rowdata <- function(all.ranges) {
    all.rd <- lapply(all.ranges, mcols)
    all.ranges <- lapply(all.ranges, FUN=function(x) {
        mcols(x) <- NULL
        x
    })
  
    out.ranges <- all.ranges[[1]]
    for (i in seq_len(length(all.ranges)-1L)) {
        if (!identical(out.ranges, all.ranges[[i+1]])) {
            out.ranges <- NULL
            warning("ignoring non-identical 'rowRanges' fields")
            break
        }
    }

    all.names <- lapply(all.rd, colnames) 
    universe <- Reduce(union, all.names)
    existing <- vector("list", length(universe))
    names(existing) <- universe
    blacklisted <- character(0)

    for (i in seq_along(all.rd)) {
        current <- all.rd[[i]]
        for (j in colnames(current)) {
            values <- current[[j]]
            previous <- existing[[j]]

            if (is.null(previous)) {
                existing[[j]] <- values
            } else if (!identical(previous, values)) {
                warning("ignoring non-identical '", j, "' field in 'rowData'")
                blacklisted <- c(blacklisted, j)
            }
        }
    }

    existing <- existing[setdiff(names(existing), blacklisted)]
    if (length(existing)) {
        existing <- lapply(existing, I)
        out.df <- DataFrame(existing)
    } else {
        out.df <- make_zero_col_DFrame(length(all.ranges[[1]]))
    }

    list(ranges=out.ranges, df=out.df)
}

#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SingleCellExperiment reducedDims<- reducedDims reducedDimNames 
#' @importFrom SummarizedExperiment assays assays<- assayNames rowData rowData<- colData colData<-
.add.single_sce <- function(original, merged, subset.row=NULL, correct.all=FALSE) {
    if (!correct.all && !is.null(subset.row)) {
        original <- original[subset.row,]
    }

    combine.assays <- .eliminate_overlaps(assayNames(merged), assayNames(original), msg="'assays'")
    assays(original) <- c(assays(merged), assays(original)[combine.assays])

    combine.reddim <- .eliminate_overlaps(reducedDimNames(merged), reducedDimNames(original), msg="'reducedDims'")
    reducedDims(original) <- c(reducedDims(merged), reducedDims(original)[combine.reddim])

    combine.coldata <- .eliminate_overlaps(colnames(colData(merged)), colnames(colData(original)), msg="'colData' fields")
    colData(original) <- cbind(colData(merged), colData(original)[,combine.coldata,drop=FALSE])

    combine.rowdata <- .eliminate_overlaps(colnames(rowData(merged)), colnames(rowData(original)), msg="'rowData' fields")
    rowData(original) <- cbind(rowData(merged), rowData(original)[,combine.rowdata,drop=FALSE])

    combine.metadata <- .eliminate_overlaps(names(metadata(merged)), names(metadata(original)), msg="'metadata'")
    metadata(original) <- c(metadata(merged), metadata(original)[combine.metadata])

    original
}
