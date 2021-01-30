#' Convert corrected PCs to a SingleCellExperiment
#'
#' Convert low-dimensional corrected PCs to a \linkS4class{SingleCellExperiment} containing corrected expression values.
#' This is a low-level function and most users should not need to call it.
#' 
#' @param corrected.df A \linkS4class{DataFrame} containing a nested matrix of low-dimensional \code{corrected} values and a vector of \code{batch} identities.
#' Typically produced from \code{\link{reducedMNN}}.
#' @param pc.info A list containing PCA statistics, in particular a \code{rotation} matrix.
#' Typically obtained from the \code{\link{metadata}} of the output from \code{\link{multiBatchPCA}}.
#' @param assay.name String specifying the name of the assay to use to store the corrected expression values.
#' @param dimred.name String containing the name fo the \code{\link{reducedDims}} to store the low-dimensional corrected values.
#' Set to \code{NULL} to avoid storing these. 
#'
#' @return A \linkS4class{SingleCellExperiment} containing a \linkS4class{LowRankMatrix} with the corrected per-gene expression values.
#' The \code{\link{colData}} contains the batch identities, the \code{\link{rowData}} contains the rotation matrix,
#' and the \code{\link{reducedDims}} contains the low-dimensional corrected values (if \code{dimred.name} is not \code{NULL}).
#' All additional metadata from \code{corrected.df} and \code{pc.info} is stored in \code{\link{metadata}}.
#' 
#' @author Aaron Lun
#' 
#' @details
#' The corrected expression values are obtained by simply taking the crossproduct of the corrected PCs with the rotation matrix.
#' This reverses the original projection to PC space while retaining the effect of the correction.
#' These values are best used for visualization; 
#' the low-dimensional corrected coordinates are more efficient for per-cell operations like clustering,
#' while the original uncorrected expression values are safer to interpret for per-gene analyses.
#'
#' @seealso
#' \code{\link{reducedMNN}}, to compute \code{corrected.df}; and \code{\link{multiBatchPCA}}, to compute \code{pc.info}.
#' 
#' \code{\link{fastMNN}}, which uses this function to obtain low-rank corrected values.
#'
#' @examples
#' B1 <- matrix(rnorm(10000), nrow=50) # Batch 1
#' B2 <- matrix(rnorm(10000), nrow=50) # Batch 2
#'
#' # Equivalent to fastMNN().
#' cB1 <- cosineNorm(B1)
#' cB2 <- cosineNorm(B2)
#' pcs <- multiBatchPCA(cB1, cB2)
#' mnn.out <- reducedMNN(pcs[[1]], pcs[[2]])
#' 
#' sce <- convertPCsToSCE(mnn.out, metadata(pcs))
#' sce
#'
#' @export
#' @importFrom S4Vectors metadata DataFrame 
#' @importFrom BiocSingular LowRankMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment
convertPCsToSCE <- function(corrected.df, pc.info, assay.name="reconstructed", dimred.name="corrected") {
    rot <- pc.info$rotation
    mat <- LowRankMatrix(rot, corrected.df$corrected)
    assays <- list(mat)
    names(assays) <- assay.name

    rdf <- DataFrame(rotation=I(rot))
    pc.info$rotation <- NULL # no need to store this in the metadata.

    all.meta <- metadata(corrected.df)
    all.meta$pca.info <- pc.info

    dimreds <- list()
    if (!is.null(dimred.name)) {
        dimreds[[dimred.name]] <- corrected.df$corrected
    }

    SingleCellExperiment(assays,
        colData=DataFrame(batch=corrected.df$batch),
        rowData=rdf,
        metadata=all.meta,
        reducedDims=dimreds)
}
