#' Apply function over multiple SingleCellExperiments
#'
#' A generalization of \code{\link{applySCE}} to apply a function to corresponding parts of multiple \link{SingleCellExperiment}s,
#' each of which have one or more alternative Experiments.
#'
#' @inheritParams multiBatchNorm
#' @param FUN Any function that accepts multiple SummarizedExperiment or SingleCellExperiment objects.
#' @param WHICH A character or integer vector containing the names or positions of alternative Experiments to loop over.
#' Defaults to all alternative Experiments that are present in each SingleCellExperiment in \code{...}.
#' @param COMMON.ARGS Further (named) arguments to pass to all calls to \code{FUN}.
#' @param MAIN.ARGS A named list of arguments to pass to calls to \code{FUN} involving the main Experiment(s) only.
#' Alternatively \code{NULL}, in which case the function is \emph{not} applied to the main Experiment.
#' @param ALT.ARGS A named list where each entry is named after an alternative Experiment and contains named arguments to use in \code{FUN} for that Experiment.
#' @param SIMPLIFY Logical scalar indicating whether the output should be simplified to one or more SingleCellExperiments.
#' 
#' @return 
#' In most cases or when \code{SIMPLIFY=FALSE}, a list is returned containing the output of \code{FUN} applied to each \emph{corresponding} Experiment across all \code{...}.
#' If \code{MAIN.ARGS} is not \code{NULL}, the first entry corresponds to the result generated from the main Experiments;
#' all other results are generated according to the entries specified in \code{WHICH} and are named accordingly.
#' 
#' If \code{SIMPLIFY=TRUE} and certain conditions are fulfilled, we can either return:
#' \itemize{
#' \item A single SingleCellExperiment, if all calls to \code{FUN} return a SingleCellExperiment.
#' Here, the results of \code{FUN} on the main/alternative Experiments in \code{...} are mapped to the main or alternative Experiments of the same name in the output.
#' \item A list of SingleCellExperiments, if all calls to \code{FUN} return a list of SingleCellExperiments of the same length.
#' The \code{\link{altExps}} of each output SingleCellExperiment contains the results from the corresponding call to \code{FUN} on the alternative Experiments of the same name in \code{...}.
#' }
#' In both cases, the aim is to mirror the organization of Experiments in each entry of \code{...}.
#'
#' @details
#' This function is a generalization of \code{\link{applySCE}} whereby corresponding Experiments from \code{...} are passed to \code{FUN}.
#' To illustrate, if we passed objects \code{x}, \code{y} and \code{z} in \code{...}:
#' \enumerate{
#' \item We first call \code{FUN} on the set of all main Experiments from \code{...}, obtaining a result equivalent to \code{FUN(x, y, z)} (more on the other arguments later).
#' \item Then we call \code{FUN} on the set of all first alternative Experiments.
#' This is equivalent to \code{FUN(altExp(x), altExp(y), altExp(z))}.
#' \item Then we call \code{FUN} on the set of all second alternative Experiments.
#' This is equivalent to \code{FUN(altExp(x, 2), altExp(y, 2), altExp(z, 2))}.
#' \item And so on.
#' }
#' In effect, much like \code{\link{applySCE}} is analogous to \code{\link{lapply}}, \code{applyMultiSCE} is analogous to \code{\link{mapply}}.
#' This allows users to easily apply the same function to all the Experiments (main and alternative) in a list of \linkS4class{SingleCellExperiment} objects.
#' 
#' Arguments in \code{COMMON.ARGS} (plus some extra arguments, see below) are passed to all calls to \code{FUN}.
#' Arguments in \code{MAIN.ARGS} are only used in the call to \code{FUN} on the main Experiments.
#' Arguments in \code{ALT.ARGS} are passed to the call to \code{FUN} on the alternative Experiments of the same name.
#' For the last two, any arguments therein will override arguments of the same name in \code{COMMON.ARGS}.
#'
#' Arguments in \code{...} that are \emph{not} SingleCellExperiments are actually treated as additional arguments for \code{COMMON.ARGS}.
#' This is purely intended as a user convenience, to avoid the need to write \code{COMMON.ARGS=list()} when specifying these arguments.
#' However, explicitly using \code{COMMON.ARGS} is the safer approach and recommended for developers.
#'
#' By default, looping is performed over alternative Experiments with names that are present across all entries of \code{...}.
#' Values of \code{WHICH} should be unique if any simplification of the output is desired.
#' If \code{MAIN.ARGS=NULL}, the main Experiment is ignored and the function is only applied to the alternative Experiments.
#'
#' The default of \code{SIMPLIFY=TRUE} is aims to make the output easier to manipulate.
#' If \code{FUN} returns a SingleCellExperiment, the outputs across main and alternative Experiments are simplified into a single SingleCellExperiment.
#' If \code{FUN} returns a list of SingleCellExperiments of the same length, the outputs are simplified into one list of SingleCellExperiments.
#' This assumes that \code{WHICH} contains no more than one reference to each alternative Experiment in \code{x}.
#'
#' @seealso
#' \code{\link{applySCE}}, for the simpler version that involves only one SingleCellExperiment object.
#'
#' \code{\link{simplifyToSCE}}, for the conditions required for simplification.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Setting up some objects with alternative Experiments.
#' d1 <- matrix(rnbinom(50000, mu=10, size=1), ncol=100)
#' sce1 <- SingleCellExperiment(list(counts=d1))
#' sizeFactors(sce1) <- runif(ncol(d1))
#' altExp(sce1, "Spike") <- sce1
#' altExp(sce1, "Protein") <- sce1
#' 
#' d2 <- matrix(rnbinom(20000, mu=50, size=1), ncol=40)
#' sce2 <- SingleCellExperiment(list(counts=d2))
#' sizeFactors(sce2) <- runif(ncol(d2))
#' altExp(sce2, "Spike") <- sce2
#' altExp(sce2, "Protein") <- sce2
#'
#' # Applying a function over the main and alternative experiments.
#' normed <- applyMultiSCE(sce1, sce2, FUN=multiBatchNorm)
#' normed
#' altExp(normed[[1]]) # contains log-normalized values
#'
#' regressed <- applyMultiSCE(normed, FUN=regressBatches)
#' regressed
#' altExp(regressed) # contains corrected expression values
#'
#' rescaled <- applyMultiSCE(normed, FUN=rescaleBatches)
#' rescaled
#' altExp(rescaled) # contains corrected expression values
#'
#' # We can also specify 'batch=' directly.
#' combined <- cbind(sce1, sce2)
#' batch <- rep(1:2, c(ncol(sce1), ncol(sce2)))
#'
#' normed <- applyMultiSCE(combined, batch=batch, FUN=multiBatchNorm)
#' normed
#' altExp(normed) # contains log-normalized values
#'
#' regressed <- applyMultiSCE(normed, batch=batch, FUN=regressBatches)
#' regressed
#' altExp(regressed) # contains corrected expression values
#'
#' rescaled <- applyMultiSCE(normed, batch=batch, FUN=rescaleBatches)
#' rescaled
#' altExp(rescaled) # contains corrected expression values
#' 
#' @export
#' @importFrom scuttle .unpackLists
#' @importFrom SingleCellExperiment altExp altExpNames simplifyToSCE
applyMultiSCE <- function(..., FUN, WHICH=NULL, COMMON.ARGS=list(), MAIN.ARGS=list(), ALT.ARGS=list(), SIMPLIFY=TRUE) {
    batches <- .unpackLists(...)
    is.sce <- vapply(batches, is, class="SingleCellExperiment", TRUE)
    COMMON.ARGS <- .dedup_args(COMMON.ARGS, batches[!is.sce])
    batches <- batches[is.sce]

    # Choosing WHICH.
    if (is.null(WHICH)) {
        WHICH <- Reduce(intersect, lapply(batches, altExpNames))
    }

    use.main <- !is.null(MAIN.ARGS)
    output <- vector("list", use.main + length(WHICH))

    ae.names <- WHICH
    if (!is.character(WHICH)) {
        all.ae.names <- lapply(batches, function(x) altExpNames(x)[WHICH])
        all.ae.names <- unique(all.ae.names)
        if (length(all.ae.names) > 1) {
            stop("integer 'WHICH' refers to different 'altExpNames' across '...'")
        }
        ae.names <- all.ae.names[[1]]
    }
    ref.ae.names <- ae.names
    if (use.main) {
        ae.names <- c("", ae.names)
    }
    names(output) <- ae.names

    # Applying on the main Experiment.
    if (use.main) {
        tryCatch({
            output[[1]] <- do.call(FUN, c(batches, .dedup_args(MAIN.ARGS, COMMON.ARGS)))
        }, error=function(err) {
            stop("'FUN' failed on the main Experiments:\n  ", conditionMessage(err))
        })
        ae.names <- c("", ae.names)
    }

    # Applying on the alternative Experiments.
    for (i in seq_along(WHICH)) {
        current <- WHICH[i]
        all.alts <- lapply(batches, altExp, e=current)

        if (is.numeric(current)) {
            cur.args <- ALT.ARGS[[ref.ae.names[current]]]
        } else {
            cur.args <- ALT.ARGS[[current]]
        }

        all.args <- c(all.alts, .dedup_args(cur.args, COMMON.ARGS))
        tryCatch({
            output[[i+use.main]] <- do.call(FUN, all.args)
        }, error=function(err) {
            stop("'FUN' failed on the alternative Experiments '", current, "':\n  ", conditionMessage(err))
        })
    }

    if (SIMPLIFY) {
        already.sce <- vapply(output, is, class="SummarizedExperiment", FUN.VALUE=TRUE)
        which.main <- if (use.main) 1L else NULL

        if (any(already.sce)) {
            attempt <- simplifyToSCE(output, which.main=which.main, warn.level=1)
            if (!is.null(attempt)){ 
                output <- attempt
            }
        } else {
            n.out <- unique(lengths(output))
            if (length(n.out)==1L) {
                success <- TRUE
                attempts <- output[[1]] # preserve names, if any.

                for (i in seq_len(n.out)) {
                    collated <- lapply(output, FUN=`[[`, i=i)
                    attempt <- simplifyToSCE(collated, which.main=which.main, warn.level=1)
                    if (is.null(attempt)) {
                        success <- FALSE
                        break
                    }
                    attempts[[i]] <- attempt
                }

                if (success) {
                    output <- attempts
                }
            } else {
                warning("failed to simplify results with variable numbers of outputs")
            }
        }
    }

    output
}

.dedup_args <- function(arg1, arg2) {
    args <- c(arg1, arg2)
    args[!duplicated(names(args))]
}
