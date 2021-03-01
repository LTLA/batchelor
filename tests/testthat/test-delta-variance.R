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
    df <- mnnDeltaVariance(B1, B2, pairs=pp)
    expect_identical(which.max(df$adjusted), 1L)

    expect_equal(df$mean, rowMeans(B1[,pp[,1]] + B2[,pp[,2] - ncol(B1)]) / 2)
    expect_equal(df$total, rowVars(B1[,pp[,1]] - B2[,pp[,2] - ncol(B1)]))
    expect_equal(df$trend, scran::fitTrendVar(df$mean, df$total)$trend(df$mean))

    df2 <- mnnDeltaVariance(B1, B2, pairs=metadata(out)$merge.info$pairs[1])
    expect_equal(df, df2)
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
