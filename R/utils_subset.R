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
