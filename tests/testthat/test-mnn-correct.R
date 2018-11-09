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

        mnn.dens <- rowSums(w[,unique(mnn2)])
        N <- tabulate(mnn2, nbins=nrow(data2))
        kernel <- t(w/(N*mnn.dens))[,mnn2]

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
    data1 <- matrix(rnorm(10000, sd=0.1), ncol=25)
    data2 <- matrix(rnorm(25000, sd=0.1), ncol=25)
    corvect <- matrix(runif(length(data2)), nrow=nrow(data2))

    # Constructing a reference function.
    REF <- function(data1, data2, cell.vect, sigma) {
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
        return(scaling)
    }

    ref <- REF(data1, data2, corvect, 1)
    test <- .Call(batchelor:::cxx_adjust_shift_variance, t(data1), t(data2), corvect, 1)
    expect_equal(ref, test)

    ref <- REF(data1, data2, corvect, 0.1)
    test <- .Call(batchelor:::cxx_adjust_shift_variance, t(data1), t(data2), corvect, 0.1)
    expect_equal(ref, test)
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
    out <- mnnCorrect(rbind(alpha, alpha), rbind(bravo, bravo), rbind(charlie, charlie), subset.row=1:nrow(alpha), correct.all=TRUE, svd.dim=2)
    ref1 <- out[1:nrow(alpha),]
    ref2 <- out[nrow(alpha)+1:nrow(alpha),]
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
    
set.seed(100042)
test_that("mnnCorrect behaves correctly with the order", {
    alpha <- matrix(rnorm(1000), ncol=100)
    bravo <- matrix(rnorm(2000), ncol=200)
    charlie <- matrix(rnorm(3000), ncol=300)

    new.order <- c(2, 3, 1)
    out <- mnnCorrect(alpha, bravo, charlie, order=new.order)
    expect_false(is.unsorted(out$batch))

    ref <- mnnCorrect(bravo, charlie, alpha)
    ref$batch <- new.order[as.integer(ref$batch)]
    re.order <- order(ref$batch)
    expect_equal(assay(ref)[,re.order], assay(out))
    expect_equal(ref$batch[re.order], as.integer(out$batch))

    # Checking the pairings in closer detail.
    pairings <- metadata(out)$pairs
    ref.pairings <- metadata(ref)$pairs
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

    # Behaves with spikes as input.
    isp <- sample(nrow(alpha), nrow(alpha)/2)
    isSpike(sceA, "ERCC") <- isp
    isSpike(sceB, "ERCC") <- isp
    isSpike(sceC, "ERCC") <- isp
    out <- mnnCorrect(sceA, sceB, sceC, assay.type=1)
    ref <- mnnCorrect(alpha[-isp,], bravo[-isp,], charlie[-isp,])
    expect_equal(ref, out)

    # Spikes and subsetting interact correctly
    i <- rbinom(nrow(alpha), 1, 0.5)==1L
    out <- mnnCorrect(sceA, sceB, sceC, assay.type=1, subset.row=i)
    keep <- setdiff(which(i), isp)
    ref <- mnnCorrect(alpha[keep,], bravo[keep,], charlie[keep,], assay.type=1)
    expect_equal(ref, out)
})

set.seed(100043)
test_that("mnnCorrect behaves with within-object batches", {
    B1 <- matrix(rnorm(10000), nrow=100) # Batch 1 
    B2 <- matrix(rnorm(20000), nrow=100) # Batch 2
    B3 <- matrix(rnorm(15000), nrow=100) # Batch 2
    
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
    pairings <- metadata(out)$pairs
    ref.pairings <- metadata(ref)$pairs
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
