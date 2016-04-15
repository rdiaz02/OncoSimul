## Evaluate whole block, as i is local to function. The third_tryE and
## fourth_tryE based on the help for try_again.

## Fails as it should
third_tryE <- local({
       i <- 3
       function() {
         i <<- i - 1
         expect_false(i > 0)
       }
})
fourth_tryE <- local({
       i <- 4
       function() {
         i <<- i - 1
         expect_false(i > 0)
       }
})
try_again(3, test_that("something", {
    third_tryE()
    fourth_tryE()
}))


## Works
third_tryE <- local({
       i <- 3
       function() {
         i <<- i - 1
         expect_false(i > 0)
       }
})
fourth_tryE <- local({
       i <- 4
       function() {
         i <<- i - 1
         expect_false(i > 0)
       }
})
try_again(4, test_that("something", {
    third_tryE()
    fourth_tryE()
}))



## plat with number of tries to get an idea
try_again(3, test_that("runif", {
    x <- runif(1)
    cat("\n x is ", x, "\n")
    expect_true(x < 0.5)
    expect_true(x > 0.25)
}))


