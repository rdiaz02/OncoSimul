inittime <- Sys.time()
cat(paste("\n Dummy empty test ", date(), "\n"))

## The default reporter will not show the number
## of tests run unless there is a failure or a warning
## But I want to see how many pass (and skip) tests I get.
## So: create one dummy skip test and one warning to force full reporting.
test_that("Dummy empty (skip) test", {
    six <- 2 * 3
})

test_that("Dummy warning to force full reporting" , {
    warning("Dummy warning")
})
