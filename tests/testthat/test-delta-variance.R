# This checks that the delta variance calculations are performed properly. 
# library(testthat); library(batchelor); source("test-delta-variance.R")

set.seed(100)
means <- 2^rexp(200) * 10
B1 <- matrix(rpois(10000, means), ncol=50) # Batch 1 
B2 <- matrix(rpois(10000, means), ncol=50) # Batch 2

B1[1:2,1:25] <- B1[1:2,1:25] * 50 
B2[2,1:25] <- B2[2,1:25] * 50 

library(scuttle)
B1 <- normalizeCounts(B1)
B2 <- normalizeCounts(B2)
out <- fastMNN(B1, B2)

test_that("mnnDeltaVariance works as expected in basic cases", {
    pp <- metadata(out)$merge.info$pairs[[1]]
    df <- mnnDeltaVariance(B1, B2, pairs=pp, cos.norm=FALSE)
    expect_identical(which.max(df$adjusted), 1L)

    expect_equal(df$mean, rowMeans(B1[,pp[,1]] + B2[,pp[,2] - ncol(B1)]) / 2)
    expect_equal(df$total, rowVars(B1[,pp[,1]] - B2[,pp[,2] - ncol(B1)]))
    expect_equal(df$trend, scran::fitTrendVar(df$mean, df$total)$trend(df$mean))

    df2 <- mnnDeltaVariance(B1, B2, pairs=metadata(out)$merge.info$pairs[1], cos.norm=FALSE)
    expect_equal(df, df2)

    # Subsetting works as expected.
    df <- mnnDeltaVariance(B1[1:100,], B2[1:100,], pairs=pp, cos.norm=FALSE)
    df2 <- mnnDeltaVariance(B1, B2, pairs=pp, subset.row=1:100, cos.norm=FALSE)
    expect_equal(df, df2)

    i <- c(1:100, 1:10)
    df2 <- mnnDeltaVariance(B1[i,], B2[i,], pairs=pp, subset.row=1:100, cos.norm=FALSE, compute.all=TRUE)
    expect_equal(df, df2[1:100,])
    expect_equal(df2[101:110,], df2[1:10,])
})

test_that("mnnDeltaVariance works with the cosine normalization", {
    pp <- metadata(out)$merge.info$pairs[[1]]
    df <- mnnDeltaVariance(B1, B2, pairs=pp, cos.norm=TRUE)
    expect_identical(which.max(df$adjusted), 1L)

    df2 <- mnnDeltaVariance(cosineNorm(B1), cosineNorm(B2), pairs=pp, cos.norm=TRUE)
    expect_true(mad(df$total/df2$total) < 1e-8)

    # Interacts sanely with the subsetting.
    df <- mnnDeltaVariance(B1[1:100,], B2[1:100,], pairs=pp, cos.norm=TRUE)
    i <- c(1:100, 1:10)
    df2 <- mnnDeltaVariance(B1[i,], B2[i,], pairs=pp, subset.row=1:100, cos.norm=TRUE, compute.all=TRUE)
    expect_equal(df$adjusted, df2$adjusted[1:100])
    expect_equal(df2$adjusted[1:10], df2$adjusted[101:110])
})

test_that("mnnDeltaVariance works as expected with blocking", {
    pp <- metadata(out)$merge.info$pairs
    df <- mnnDeltaVariance(B1, B2, pairs=pp)
    df2 <- mnnDeltaVariance(B1, B2, pairs=c(pp, pp))
    expect_equal(df$adjusted, df2$adjusted)
    expect_equal(df, df2$per.step[,1])

    # A more interesting example.
    B3 <- matrix(rpois(10000, means), ncol=50) # Batch 2
    B3[3,26:50] <- B1[3,26:50] * 50 
    B3 <- normalizeCounts(B3)
    out <- fastMNN(B1, B2, B3)

    pp <- metadata(out)$merge.info$pairs
    df3 <- mnnDeltaVariance(B1, B2, B3, pairs=pp)
    expect_identical(ncol(df3$per.step), 2L)
    expect_false(isTRUE(all.equal(df2$adjusted, df3$adjusted)))
})
