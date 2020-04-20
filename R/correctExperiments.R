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
#' Any column metadata fields that are shared will also be included in the merged object if \code{combine.coldata=TRUE}.
#' If any existing field has the same name as any \code{\link{colData}} field produced by \code{\link{batchCorrect}},
#' it will be ignored in favor of the latter.
#'
#' Row metadata from \code{...} is included in the merged object if \code{include.rowdata=TRUE}.
#' In such cases, only non-conflicting row data fields are preserved,
#' i.e., fields with different names or identically named fields with the same values between objects in \code{...}.
#' Any conflicting fields are ignored with a warning.
#' \code{rowRanges} are only preserved if they are identical (ignoring the \code{\link{mcols}}) for all objects in \code{...}.
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
#' 
#' @export
#' @importFrom SummarizedExperiment assayNames assay<- assay colData<- colData 
#' rowRanges rowRanges<- rowData rowData<-
#' @importFrom S4Vectors mcols mcols<-
correctExperiments <- function(..., batch=NULL, restrict=NULL, subset.row=NULL, correct.all=FALSE, assay.type="logcounts", 
    PARAM=FastMnnParam(), combine.assays=NULL, combine.coldata=NULL, include.rowdata=TRUE) 
{
    x <- .unpack_batches(...)
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
        shared <- intersect(colnames(merged.rd), colnames(combined.rd))
        if (length(shared)) {
            warning("ignoring 'rowData' fields overlapping 'batchCorrect' output")
            combined.rd <- combined.rd[,setdiff(colnames(combined.rd), shared),drop=FALSE]
        }

        if (!is.null(combined.rr)) {
            rowRanges(merged) <- combined.rr
        }
        rowData(merged) <- cbind(merged.rd, combined.rd)
    }

    merged
}

#' @importFrom S4Vectors mcols mcols<- DataFrame make_zero_col_DFrame
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
    out.df <- if (length(existing)) {
        existing <- lapply(existing, I)
        DataFrame(existing)
    } else {
        make_zero_col_DFrame(length(all.ranges[[1]]))
    }

    list(ranges=out.ranges, df=out.df)
}
