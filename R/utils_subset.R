#' @importFrom S4Vectors normalizeSingleBracketSubscript
.row_subset_to_index <- function(x, index) {
    if (is.null(index)) {
        seq_len(nrow(x))
    } else {
        normalizeSingleBracketSubscript(index, x)
    }
}

.col_subset_to_index <- function(x, index) {
    if (is.null(index)) {
        seq_len(ncol(x))
    } else {
        i <- seq_len(ncol(x))
        names(i) <- colnames(x)
        unname(i[index])
    }
}

#' @importFrom SingleCellExperiment isSpike
.spike_subset <- function(x, get.spikes) 
# Returns a logical vector specifying which rows we should retain,
# if depending on whether or not we want to keep the spike-ins.
{
    if (!get.spikes) {
        suppressWarnings(nokeep <- isSpike(x))
        if (!is.null(nokeep) && any(nokeep)) {
            return(!nokeep)
        }
    } 
    return(NULL)
}

#' @importFrom S4Vectors normalizeSingleBracketSubscript
.SCE_subset_genes <- function(subset.row, x, get.spikes) 
# Convenience function to intersect arbitrary subsetting specification with spike-in information.
# To be mainly used for SingleCellExperiments to further restrict subset.row.
{
    despiked <- .spike_subset(x, get.spikes)
    if (is.null(subset.row)) { 
        subset.row <- despiked
    } else {
        subset.row <- normalizeSingleBracketSubscript(subset.row, x)
        if (!is.null(despiked)) { 
            subset.row <- intersect(subset.row, which(despiked))
        }
    }
    return(subset.row)
}
