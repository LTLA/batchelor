# General getters.

.get_batches <- function(x) x@batches

.get_restrict <- function(x) x@restrict

.get_reference <- function(x) x@reference

.get_reference_indices <- function(x) x@reference.indices

.get_reference_restrict <- function(x) x@reference.restrict

.get_mnn_result <- function(x) x@mnn.result

.get_current_index <- function(x) x@current.index

# Other utilities 

.do_restrict <- function(x) !is.null(.get_restrict(x))

.coerce_restrict <- function(x, b) {
    out <- .get_restrict(x)[[b]]
    if (!is.null(out)) {
        out
    } else {
        seq_len(nrow(.get_batches(x)[[b]]))
    }
}

.check_clean <- function(x) {
    if (length(.get_current_index(x)) || length(.get_mnn_result(x))) {
        stop("next results have not been assimilated yet")
    }
    invisible(NULL)
}

.check_unclean <- function(x) {
    if (length(.get_current_index(x))!=1L || length(.get_mnn_result(x))!=2L) {
        stop("next results are not available yet")
    }
    invisible(NULL)
}

.unrestrict_indices <- function(index, restrict) {
    if (!is.null(restrict)) index <- restrict[index]
    index
}

.get_current <- function(x) .get_batches(x)[[.get_current_index(x)]]

.get_current_restrict <- function(x) .get_restrict(x)[[.get_current_index(x)]]

#############################################
# This class is designed to hold information across iterations of MNN correction
# for a specified merge order across batches. Checks for correct input are 
# therefore minimal, given that this class is never user-visible.

.get_ordering <- function(x) x@ordering

#' @importFrom methods new
MNN_supplied_order <- function(batches, restrict=NULL, ordering=seq_along(batches)) {
	new("MNN_supplied_order", batches=batches, restrict=restrict, 
        reference=NULL, reference.restrict=NULL, reference.indices=integer(0),
        current.index=integer(0), mnn.result=list(),
        ordering=ordering)
}

#' @importFrom BiocNeighbors queryKNN buildIndex 
setMethod(".advance", "MNN_supplied_order", function(x, k, BNPARAM, BPPARAM) {
    .check_clean(x)
    ordering <- .get_ordering(x)
    position <- length(.get_reference_indices(x))+1L
    batches <- .get_batches(x)

    if (position==1L) { # Initializing results for the first time. 
        ref.idx <- ordering[position]
        x@reference.indices <- ref.idx
        x@reference <- batches[[ref.idx]]
        if (.do_restrict(x)) {
            x@reference.restrict <- .coerce_restrict(x, ref.idx)
        }
        position <- position + 1L
    }

    cur.idx <- ordering[position]
    curdata <- batches[[cur.idx]]
    refdata <- .get_reference(x)

    # Restricting the MNN search, if desired.
    if (.do_restrict(x)) {
        curres <- .get_restrict(x)[[cur.idx]]
        if (!is.null(curres)) {
            curdata <- curdata[curres,,drop=FALSE]
        }
        refdata <- refdata[.get_reference_restrict(x),,drop=FALSE] 
    }

    vals <- findMutualNN(refdata, curdata, k1=k, k2=k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    # Correcting the indices to refer to non-restricted indices.
    if (.do_restrict(x)) {
        vals$first <- .unrestrict_indices(vals$first, .get_reference_restrict(x))
        vals$second <- .unrestrict_indices(vals$second, .get_restrict(x)[[cur.idx]])
    }

    x@current.index <- cur.idx
    x@mnn.result <- vals
    x
})

#############################################
# This class is designed to hold information across iterations of MNN correction
# for an automatically specified merge order across batches. 

.get_precomputed <- function(x) x@precomputed

.get_restricted_batches <- function(x) x@restricted.batches

#' @importFrom methods new
MNN_auto_order <- function(batches, restrict=NULL) {
    restricted.batches <- batches
    if (!is.null(restrict)) {
        for (b in seq_along(restricted.batches)) {
            curres <- restrict[[b]]
            if (!is.null(curres)) {
                restricted.batches[[b]] <- restricted.batches[[b]][curres,,drop=FALSE]
            }
        }
    }

	new("MNN_auto_order", batches=batches, restrict=restrict, 
        reference=NULL, reference.restrict=NULL, reference.indices=integer(0),
        current.index=integer(0), mnn.result=list(),
        precomputed=list(), restricted.batches=restricted.batches)
}

#' @importFrom BiocNeighbors queryKNN buildIndex 
.init_auto_order <- function(x, k, BNPARAM, BPPARAM) {
    restricted.batches <- .get_restricted_batches(x)

    precomputed <- list()
    for (b in seq_along(restricted.batches)) {
        precomputed[[b]] <- buildIndex(restricted.batches[[b]], BNPARAM=BNPARAM)
    }
    x@precomputed <- precomputed

    max.pairs <- list()
    for (First in seq_along(precomputed)) {
        fdata <- restricted.batches[[First]]
        for (Second in seq_len(First-1L)) {
            sdata <- restricted.batches[[Second]]

            W21 <- queryKNN(BNINDEX=precomputed[[Second]], query=fdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
            W12 <- queryKNN(BNINDEX=precomputed[[First]], query=sdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
            out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
            names(out) <- c("first", "second")

            if (length(out$first) > length(max.pairs$first))  {
                max.pairs <- out
                max.first <- First
                max.second <- Second
            }
        }
    }

    x@reference <- .get_batches(x)[[max.first]]
    x@reference.indices <- max.first

    # Correcting the indices to refer to non-restricted indices.
    if (.do_restrict(x)) {
        x@reference.restrict <- .coerce_restrict(x, max.first)
        restrict <- .get_restrict(x)
        max.pairs$first <- .unrestrict_indices(max.pairs$first, restrict[[max.first]])
        max.pairs$second <- .unrestrict_indices(max.pairs$second, restrict[[max.second]])
    }

    x@mnn.result <- max.pairs
    x@current.index <- max.second

    x
}

#' @importFrom BiocNeighbors queryKNN buildIndex 
.extend_auto_order <- function(x, k, BNPARAM, BPPARAM) {
    refdata <- .get_reference(x)
    if (.do_restrict(x)) {
        refdata <- refdata[.get_reference_restrict(x),,drop=FALSE]
    }
    precomp.ref <- buildIndex(refdata, BNPARAM=BNPARAM)
   
    restricted.batches <- .get_restricted_batches(x)
    precomputed <- .get_precomputed(x)
    max.pairs <- list()
    processed <- .get_reference_indices(x)

    for (other in seq_along(restricted.batches)[-processed]) {
        W21 <- queryKNN(BNINDEX=precomputed[[other]], query=refdata, k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        W12 <- queryKNN(BNINDEX=precomp.ref, query=restricted.batches[[other]], k=k, BPPARAM=BPPARAM, get.distance=FALSE)
        out <- .Call(cxx_find_mutual_nns, W21$index, W12$index)
        names(out) <- c("first", "second")
        
        if (length(out$first) > length(max.pairs$first))  {
            max.pairs <- out
            max.other <- other
        }
    }

    # Correcting the indices to refer to non-restricted indices.
    if (.do_restrict(x)) {
        max.pairs$first <- .unrestrict_indices(max.pairs$first, .get_reference_restrict(x))
        max.pairs$second <- .unrestrict_indices(max.pairs$second, .get_restrict(x)[[max.other]])
    }

    x@mnn.result <- max.pairs
    x@current.index <- max.other

    x
}

setMethod(".advance", "MNN_auto_order", function(x, k, BNPARAM, BPPARAM) {
    .check_clean(x)
	if (length(.get_reference_indices(x))==0L) {
        .init_auto_order(x, k, BNPARAM, BPPARAM)
    } else {
        .extend_auto_order(x, k, BNPARAM, BPPARAM)
    }
})

#############################################

.compile <- function(x, corrected) {
    .check_unclean(x)

    # Updating the reference information.
    cur.dex <- .get_current_index(x)
    refN <- nrow(.get_reference(x))
    x@reference <- rbind(.get_reference(x), corrected)
    x@reference.indices <- c(.get_reference_indices(x), cur.dex)

    # Updating the restricted set in the current reference.
    if (.do_restrict(x)) {
        x@reference.restrict <- c(.get_reference_restrict(x), refN + .coerce_restrict(x, cur.dex))
    }

    # Clearing other contents.
    x@current.index <- integer(0)
    x@mnn.result <- list()
    x
}

