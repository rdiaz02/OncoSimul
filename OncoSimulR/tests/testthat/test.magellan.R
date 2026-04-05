## Minimal exercising and testing of some
## functionality from MAGELLAN

inittime <- Sys.time()
cat(paste("\n Starting test.magellan at", date(), "\n"))

test_that("Abort in NK", {
    expect_error(rfitness(4, K = 5, model = "NK"),
                 "It makes no sense to have K >= g", fixed = TRUE)
    expect_error(rfitness(6, K = 6, model = "NK"),
                 "It makes no sense to have K >= g", fixed = TRUE)
    }
)


test_that("Call Magellan stats on NK", {
    rnk1 <- rfitness(6, K = 1, model = "NK")
    expect_true(is.numeric(Magellan_stats(rnk1)))
    
    rnk2 <- rfitness(6, K = 4, model = "NK")
    expect_true(is.numeric(Magellan_stats(rnk2)))
    })


test_that("Call Magellan stats on RMF", {
    rmf1 <- rfitness(6)
    expect_true(is.numeric(Magellan_stats(rmf1)))
    })



test_that("Minimal Magellan tests", {
  r1 <- rfitness(2)

  r1_noepi <- r1
  r1_noepi[, 3] <- c(1, 2, 3, 4)

  r1_sign <- r1
  r1_sign[, 3] <- c(1, 2, 0.1, 4)

  r1_rsign <- r1
  r1_rsign[, 3] <- c(1, 3, 4, 2)

  expect_equal(Magellan_stats(r1_noepi)[c("sign", "rsign")],
               c(sign = 0, rsign = 0))
  expect_equal(Magellan_stats(r1_sign)[c("sign", "rsign")],
               c(sign = 1, rsign = 0))
  expect_equal(Magellan_stats(r1_rsign)[c("sign", "rsign")],
               c(sign = 0, rsign = 1))

  r2 <- rfitness(3)

  r2[, 4] <- c(1, 3, 1.5, 0.1, 2, 4, 1.2, 1)
  ## sign:
  ##    WT, A, B, AB
  ##    WT, A, C, AC
  ##    A, AB, AC, ABC
  ##    B, AB, BC, ABC
  ## rsign:
  ##    C, AC, BC, ABC

  expect_equal(Magellan_stats(r2)[c("sign", "rsign")],
               c(sign = 4/6, rsign = 1/6))
})


cat(paste("\n Ending test.magellan at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
