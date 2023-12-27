# Checks the application of the mnnCorrect function.
# require(batchelor); require(testthat); source("test-mnn-correct.R")

set.seed(10002)
test_that("Biological subspace is correctly estimated and re-projected", {
    # Estimation yields the PCA rotation vectors.
    A <- matrix(rnorm(10000), ncol=50)
    out <- batchelor:::.get_bio_span(A, 3) 
    ref <- prcomp(t(A))$rotation[,1:3]
    dimnames(ref) <- NULL
    expect_equal(abs(colSums(out/ref)), rep(nrow(out), ncol(out)))

    # Re-projection.
    A <- matrix(rnorm(10000), ncol=50)
    subset <- sample(nrow(A), nrow(A)/2)   
    ref <- batchelor:::.get_bio_span(A[subset,], 3) 
    proj <- batchelor:::.get_bio_span(A, 3, subset.row=subset) 
    expect_equal(ref, proj[subset,])

    B <- rbind(A, A)
    first.half <- seq_len(nrow(A))
    ref <- batchelor:::.get_bio_span(A, 3)
    proj <- batchelor:::.get_bio_span(B, 3, subset.row=first.half) 
    expect_equal(ref, proj[first.half,])
    expect_equal(proj[first.half,], proj[-first.half,])
})

set.seed(10003)
test_that("Batch vectors are correctly calculated", {
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    mnn1 <- 1:10
    mnn2 <- 30:21

    # Constructing a reference function.
    REF <- function(data1, data2, mnn1, mnn2, s2) {
        d <- as.matrix(dist(data2))
        w <- exp(-d^2/s2)

        # Normalizing for differences in MNN density around each cell in
        # dataset 2 involved in a MNN pair. (Inner transposition is not 
        # technically necessary for a symmetric matrix, but whatever.)
        mnn.dens <- rowSums(w[,unique(mnn2)])

        # Expanding for the number of times a cell in dataset 2 is involved in
        # a MNN pair, and dividing the weight appropriately. This is equivalent
        # to averaging across all pairwise vectors for a given MNN-involved
        # cell in dataset 2.
        N <- tabulate(mnn2, nbins=nrow(data2))

        # Applying the two effects above to create a kernel. Each row here
        # represents a cell in dataset 2, and each column represents a 
        # cell in dataset 2 that is in a MNN pair (possibly multiple times).
        kernel <- t(w/(N*mnn.dens))[,mnn2]

        # We want to compute the smoothed average for each cell, so we
        # normalize by the total weight across all MNN-involved cells.
        kernel <- kernel/rowSums(kernel)

        vect <- data1[mnn1,] - data2[mnn2,]
        out <- kernel %*% vect

        dimnames(out) <- NULL
        return(out)
    }

    # Vanilla check
    s2 <- 0.1
    xx <- batchelor:::.compute_correction_vectors(data1, data2, mnn1, mnn2, t(data2), s2)
    ref <- REF(data1, data2, mnn1, mnn2, s2)
    expect_equal(xx, ref)

    # Check with cells in multiple MNN pairs.
    alt.mnn1 <- c(11, 12, 13, mnn1)
    alt.mnn2 <- c(30, 30, 30, mnn2)
    xx <- batchelor:::.compute_correction_vectors(data1, data2, alt.mnn1, alt.mnn2, t(data2), s2)
    ref <- REF(data1, data2, alt.mnn1, alt.mnn2, s2)
    expect_equal(xx, ref)

    # Check with more MNN pairs involved.
    alt.mnn1 <- 1:200
    alt.mnn2 <- 500:301
    xx <- batchelor:::.compute_correction_vectors(data1, data2, alt.mnn1, alt.mnn2, t(data2), s2)
    ref <- REF(data1, data2, alt.mnn1, alt.mnn2, s2)
    expect_equal(xx, ref)

    # Check with a different bandwidth. 
    s2 <- 0.5
    xx <- batchelor:::.compute_correction_vectors(data1, data2, mnn1, mnn2, t(data2), s2)
    ref <- REF(data1, data2, mnn1, mnn2, s2)
    expect_equal(xx, ref)
})

set.seed(100032)
test_that("Variance shift adjustment is correctly performed", {
    data1 <- matrix(rnorm(10000, sd=0.1), nrow=25)
    data2 <- matrix(rnorm(25000, sd=0.1), nrow=25)
    corvect <- matrix(runif(length(data2)), nrow=ncol(data2)) # Yes, this is transposed relative to 'data2'.

    # Constructing a reference function.
    REF <- function(data1, data2, cell.vect, sigma) {
        data1 <- t(data1)
        data2 <- t(data2)
        scaling <- numeric(nrow(cell.vect))

        for (cell in seq_along(scaling)) {
            # For each cell, projecting both data sets onto the normalized correction vector for that cell.
            cur.cor.vect <- cell.vect[cell,]
            l2norm <- sqrt(sum(cur.cor.vect^2))
            cur.cor.vect <- cur.cor.vect/l2norm
            coords2 <- data2 %*% cur.cor.vect
            coords1 <- data1 %*% cur.cor.vect
    
            # Also getting the distance from the correction vector. 
            dist2 <- data2[cell,] - t(data2)
            dist2 <- dist2 - outer(cur.cor.vect, as.numeric(crossprod(dist2, cur.cor.vect)))
            dist2 <- colSums(dist2^2)
            weight2 <- exp(-dist2/sigma)
    
            dist1 <- data2[cell,] - t(data1)
            dist1 <- dist1 - outer(cur.cor.vect, as.numeric(crossprod(dist1, cur.cor.vect)))
            dist1 <- colSums(dist1^2)
            weight1 <- exp(-dist1/sigma)
    
            # Computing the weighted cumulative probability, for quantile-quantile mapping.
            rank2 <- rank(coords2, ties.method="first")
            prob2 <- sum(weight2[rank2 <= rank2[cell]])/sum(weight2)
            ord1 <- order(coords1)
            ecdf1 <- cumsum(weight1[ord1])/sum(weight1)
            
            # Adjusting the length of the correction vector so that the correction will match the quantiles.
            quan1 <- coords1[ord1[min(which(ecdf1 >= prob2))]]
            quan2 <- coords2[cell]
            scaling[cell] <- (quan1 - quan2)/l2norm
        }

        scaling 
    }

    # I dunno, man, I don't have the time to fix this right now.
    skip_on_os("mac", arch="aarch64")

    TEST <- function(data1, data2, cell.vect, sigma) {
        batchelor:::adjust_shift_variance(data1, data2, cell.vect, 
            sigma, seq_len(ncol(data1))-1L, seq_len(ncol(data2))-1L)
    }

    ref <- REF(data1, data2, corvect, 1)
    test <- TEST(data1, data2, corvect, 1)
    expect_equal(ref, test)

    ref <- REF(data1, data2, corvect, 0.1)
    test <- TEST(data1, data2, corvect, 0.1)
    expect_equal(ref, test)

    # Handles subsetting.
    i <- 10:20
    test1 <- batchelor:::.adjust_shift_variance(data1, data2, corvect, 1, subset.row=i)
    test2 <- batchelor:::.adjust_shift_variance(data1[i,], data2[i,], corvect[,i], 1)
    expect_equal(test1[,i], test2)

    # Handles restriction. 
    i1 <- 10:20
    A1 <- cbind(data1, data1[,i1])
    i2 <- 20:10
    A2 <- cbind(data2, data2[,i2])

    test1 <- batchelor:::.adjust_shift_variance(data1, data2, corvect, 1)
    test2 <- batchelor:::.adjust_shift_variance(A1, A2, rbind(corvect, corvect[i2,]), 1, restrict1=i1, restrict2=i2)

    common <- seq_len(ncol(data2))
    expect_identical(test1, test2[common,])
    expect_identical(test1[i2,], test2[-common,])
})

set.seed(10004)
test_that("mnnCorrect behaves consistently with subsetting", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)

    keep <- 1:5 
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep)
    expect_equal(ref, out)

    # Now correcting everything.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    plus.keep <- c(keep, 2:4)
    out <- mnnCorrect(alpha[plus.keep,], bravo[plus.keep,], charlie[plus.keep,], subset.row=keep, correct.all=TRUE)
    expect_equal(ref, out[keep,])
    expect_equal(ref[2:4,], out[length(keep)+2:4-1,])

    # Duplicated genes should have no effect.
    extra <- 2:5
    out <- mnnCorrect(rbind(alpha, alpha[extra,]), rbind(bravo, bravo[extra,]), rbind(charlie, charlie[extra,]), 
        subset.row=1:nrow(alpha), correct.all=TRUE, svd.dim=2)
    ref1 <- out[extra,]
    ref2 <- out[nrow(alpha)+seq_along(extra),]
    expect_equal(ref1, ref2) 
})

set.seed(100041)
test_that("mnnCorrect behaves consistently with subsetting and cosine normalization", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)

    keep <- 1:5 
    default <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])

    # Without cosine normalization of the output.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.out=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.out=FALSE)
    expect_equal(ref, out)

    expect_false(isTRUE(all.equal(assay(default), assay(out))))

    # Without cosine normalization of the input.
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.in=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.in=FALSE)
    expect_equal(ref, out)

    expect_false(isTRUE(all.equal(assay(default), assay(out))))

    # Without any cosine normalization at all.
    keep <- 6:10
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], cos.norm.in=FALSE, cos.norm.out=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, cos.norm.in=FALSE, cos.norm.out=FALSE)
    expect_equal(ref, out)

    default <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    expect_false(isTRUE(all.equal(assay(default), assay(out))))

    # With SVDs.
    keep <- 2:7
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], svd.dim=2)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, svd.dim=2)
    expect_equal(ref, out)

    default <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    expect_false(isTRUE(all.equal(assay(default), assay(out))))

    # Without the variance adjustment.
    keep <- 3:8
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], var.adj=FALSE)
    out <- mnnCorrect(alpha, bravo, charlie, subset.row=keep, var.adj=FALSE)
    expect_equal(ref, out)   

    default <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,])
    expect_false(isTRUE(all.equal(assay(default), assay(out))))
})

set.seed(120000411)
test_that("mnnCorrect uses 'prop.k' correctly", {
    B1 <- matrix(rnorm(10000, 0), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(10000, 1), nrow=100) # Batch 2

    ref <- mnnCorrect(B1, B2) 
    out <- mnnCorrect(B1, B2, k=10, prop.k=20/ncol(B1)) 
    expect_identical(assay(ref), assay(out))

    # max() kicks in.
    ref <- mnnCorrect(B1, B2) 
    out <- mnnCorrect(B1, B2, prop.k=0)
    expect_identical(assay(ref), assay(out))

    # Actually causes a difference in results.
    B2a <- matrix(rnorm(20000, 1), nrow=100)
    ref <- mnnCorrect(B1, B2a)
    out <- mnnCorrect(B1, B2a, prop.k=20/ncol(B1))
    expect_false(identical(assay(ref), assay(out)))
})

set.seed(100042)
test_that("mnnCorrect behaves correctly with an alternative order", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)

    new.order <- c(2, 3, 1)
    out <- mnnCorrect(alpha, bravo, charlie, merge.order=new.order)
    expect_false(is.unsorted(out$batch))

    ref <- mnnCorrect(bravo, charlie, alpha)
    ref$batch <- new.order[as.integer(ref$batch)]
    re.order <- order(ref$batch)
    expect_equal(assay(ref)[,re.order], assay(out))
    expect_equal(ref$batch[re.order], as.integer(out$batch))

    # Checking the pairings in closer detail.
    pairings <- metadata(out)$merge.info$pairs
    ref.pairings <- metadata(ref)$merge.info$pairs
    for (x in seq_along(pairings)) {
        curref <- ref.pairings[[x]]
        curref[,1] <- match(curref[,1], re.order)
        curref[,2] <- match(curref[,2], re.order)
        curout <- pairings[[x]]
        expect_identical(
            curref[order(curref[,1], curref[,2]),],
            curout[order(curout[,1], curout[,2]),]
        )
    }
})

set.seed(100042)
test_that("mnnCorrect behaves correctly with an automatic order", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)
    out <- mnnCorrect(alpha, bravo, charlie, merge.order=c(2,3,1))

    # Works with automatic ordering.
    auto <- mnnCorrect(A=alpha, B=bravo, C=charlie, auto.merge=TRUE)
    expect_identical(auto$batch, LETTERS[out$batch])
    expect_identical(metadata(auto)$merge.info$left[[1]], "C") # charlie and bravo would have most MNNs.
    expect_identical(metadata(auto)$merge.info$right[[1]], "B") 
    expect_identical(metadata(auto)$merge.info$left[[2]], c("C", "B"))
    expect_identical(metadata(auto)$merge.info$right[[2]], "A")

    # Automatic ordering works with options that force same.set=FALSE
    extra <- 5:1
    auto2 <- mnnCorrect(A=rbind(alpha, alpha[extra,]),
        B=rbind(bravo, bravo[extra,]), 
        C=rbind(charlie, charlie[extra,]), 
        auto.merge=TRUE, subset.row=1:10, correct.all=TRUE)

    expect_identical(assay(auto2)[1:nrow(alpha),], assay(auto))
    expect_identical(assay(auto2)[nrow(alpha)+seq_along(extra),], assay(auto)[extra,])
})

set.seed(100043)
test_that("mnnCorrect behaves with SingleCellExperiment inputs", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)
    ref <- mnnCorrect(alpha, bravo, charlie)

    sceA <- SingleCellExperiment(alpha)
    sceB <- SingleCellExperiment(bravo)
    sceC <- SingleCellExperiment(charlie)
    out <- mnnCorrect(sceA, sceB, sceC, assay.type=1)
    expect_equal(ref, out)
})

set.seed(100044)
test_that("mnnCorrect behaves with within-object batches", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000), nrow=100) # Batch 3
    
    combined <- cbind(B1, B2, B3)
    batches <- rep(1:3, c(ncol(B1), ncol(B2), ncol(B3)))
    
    shuffle <- sample(ncol(combined))
    combined <- combined[,shuffle]
    batches <- batches[shuffle]
    
    ref <- mnnCorrect(B1, B2, B3)
    out <- mnnCorrect(combined, batch=batches)
    expect_equal(assay(ref)[,shuffle], assay(out))
    expect_equal(as.character(ref$batch)[shuffle], as.character(out$batch))

    # Checking the pairings in closer detail.
    pairings <- metadata(out)$merge.info$pairs
    ref.pairings <- metadata(ref)$merge.info$pairs
    for (x in seq_along(pairings)) {
        curref <- ref.pairings[[x]]
        curref[,1] <- match(curref[,1], shuffle)
        curref[,2] <- match(curref[,2], shuffle)
        curout <- pairings[[x]]
        expect_identical(
            curref[order(curref[,1], curref[,2]),],
            curout[order(curout[,1], curout[,2]),]
        )
    }
})

set.seed(100045)
test_that("mnnCorrect behaves with restriction", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000, 2), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000, 3), nrow=100) # Batch 3

    i1 <- 20:10
    C1 <- cbind(B1, B1[,i1])
    i2 <- 30:100
    C2 <- cbind(B2, B2[,i2])
    i3 <- 100:20
    C3 <- cbind(B3, B3[,i3])

    keep1 <- seq_len(ncol(B1))
    keep2 <- seq_len(ncol(B2))
    keep3 <- seq_len(ncol(B3))

    # Windows 32-bit fails with var.adj=TRUE, 
    # due to precision issues with an exact comparison
    # when determining the quantile for variance matching.
    N <- if (.Platform$OS.type=="unix") 4L else 1L 

    for (it in seq_len(N)) {
        if (it==1L) {
            args <- list(var.adj=FALSE)
        } else if (it==2L) {
            args <- list()
        } else if (it==3L) {
            args <- list(svd.dim=2L)
        } else {
            # To trigger same.set=FALSE. Do NOT set cos.norm.out=TRUE, as 
            # the large distances result in kernels of 1 (basically only itself)
            # which causes a lot of ties that are broken randomly due to 'sample' below.
            args <- list(subset.row=50:1, correct.all=TRUE) 
        }

        ref <- do.call(mnnCorrect, c(list(B1, B2, B3), args))
        out <- do.call(mnnCorrect, c(list(C1, C2, C3, restrict=list(keep1, keep2, keep3)), args))

        expect_identical(assay(ref)[,ref$batch==1], assay(out)[,out$batch==1][,keep1])
        expect_identical(assay(ref)[,ref$batch==2], assay(out)[,out$batch==2][,keep2])
        expect_identical(assay(ref)[,ref$batch==3], assay(out)[,out$batch==3][,keep3])
    
        expect_identical(assay(ref)[,ref$batch==1][,i1], assay(out)[,out$batch==1][,-keep1])
        expect_identical(assay(ref)[,ref$batch==2][,i2], assay(out)[,out$batch==2][,-keep2])
        expect_identical(assay(ref)[,ref$batch==3][,i3], assay(out)[,out$batch==3][,-keep3])
    
        # Also behaves with a single batch.
        DY <- cbind(C1, C2, C3)
        batch <- rep(1:3, c(ncol(C1), ncol(C2), ncol(C3)))
    
        shuffle <- sample(ncol(DY))
        out2 <- do.call(mnnCorrect, c(
            list(DY[,shuffle], batch=batch[shuffle], 
                restrict=list(shuffle %in% c(keep1, keep2 + ncol(C1), keep3 + ncol(C1) + ncol(C2)))
            ),
            args
        ))
    
        expect_equal(assay(out2), assay(out)[,shuffle])
        expect_identical(out2$batch, as.character(out$batch)[shuffle])
    }
})

set.seed(100046)
test_that("mnnCorrect preserves names in its output", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2

    out <- mnnCorrect(X=B1, Y=B2) 
    expect_identical(out$batch, rep(c("X", "Y"), c(ncol(B1), ncol(B2))))

    rownames(B1) <- rownames(B2) <- sprintf("GENE_%i", seq_len(nrow(B1)))
    out <- mnnCorrect(B1, B2) 
    expect_identical(rownames(out), rownames(B1))

    out <- mnnCorrect(B1, B2, subset.row=1:50) 
    expect_identical(rownames(out), rownames(B1)[1:50])

    colnames(B1) <- sprintf("CELL_%i", seq_len(ncol(B1)))
    colnames(B2) <- sprintf("CELL_%i", seq_len(ncol(B2)))
    out <- mnnCorrect(B1, B2) 
    expect_identical(colnames(out), c(colnames(B1), colnames(B2)))

    # Row names are handled correctly with correct.all=TRUE.
    C1 <- B1
    C2 <- B2
    rownames(C1) <- rownames(C2) <- 1:nrow(C1)

    keep <- 10:100
    sub <- mnnCorrect(C1, C2, subset.row=keep)
    expect_identical(rownames(sub), as.character(keep))

    all <- mnnCorrect(C1, C2, subset.row=keep, correct.all=TRUE)
    expect_identical(rownames(all), rownames(C1))
    expect_identical(assay(all)[keep,], assay(sub))
})

set.seed(10005)
library(Matrix)
test_that("mnnCorrect behaves properly with sparse matrices", {
    alpha <- rsparsematrix(20, 100, density=0.5)
    bravo <- rsparsematrix(20, 100, density=0.5)
    charlie <- rsparsematrix(20, 100, density=0.5)

    out <- mnnCorrect(alpha, bravo, charlie)
    ref <- mnnCorrect(as.matrix(alpha), as.matrix(bravo), as.matrix(charlie))
    expect_equivalent(ref, out)
})
