## A minimal thing, to make sure no screw ups

epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "d:c"), s = c(0.2, 0.3))

oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)
oeffects1b <- data.frame(ids = c("d>a", "c > d"),
                         s = c(0.4, -0.3))

oeffects1c <- data.frame(ids = c("a>d", "c > d"),
                         s = c(0.4, -0.3))



test_that("epist to long identical", {
    expect_identical(to.long.epist.order(epistm1, ":"),
                     to.long.epist.order(epistm1b, ":"))
})

test_that("order effects to long identical", {
    expect_identical(to.long.epist.order(oeffects1, ">"),
                     to.long.epist.order(oeffects1b, ">"))
})

test_that("order effects to long different", {
    expect_false(identical(to.long.epist.order(oeffects1, ">"),
                     to.long.epist.order(oeffects1c, ">")))
})


