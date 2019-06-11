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
    Magellan_stats(rnk1)
    
    rnk2 <- rfitness(6, K = 4, model = "NK")
    Magellan_stats(rnk2)
    })


test_that("Call Magellan stats on RMF", {
    rmf1 <- rfitness(6)
    Magellan_stats(rmf1)
    
    })


cat(paste("\n Ending test.magellan at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
