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



## play with number of tries to get an idea
try_again(3, test_that("runif", {
    x <- runif(1)
    cat("\n x is ", x, "\n")
    expect_true(x < 0.5)
    expect_true(x > 0.25)
}))


try_again2 <- function (times, code, message_times = TRUE) 
{
    init_times <- times
    while (times > 0) {
        e <- tryCatch(withCallingHandlers({
            code
            NULL
        }, warning = function(e) {
            if (identical(e$message, "restarting interrupted promise evaluation")) {
                invokeRestart("muffleWarning")
            }
        }), expectation_failure = function(e) {
            e
        }, error = function(e) {
            e
        })
        if (is.null(e)) {
            if(message_times)
                message("\n times run: ", init_times - times + 1)
            return(invisible(TRUE))
        }
        times <- times - 1L
    }
    stop(e)
}


try_again2(3, test_that("runif", {
    x <- runif(1)
    cat("\n x is ", x, "\n")
    expect_true(x < 0.5)
    expect_true(x > 0.25)
}))
