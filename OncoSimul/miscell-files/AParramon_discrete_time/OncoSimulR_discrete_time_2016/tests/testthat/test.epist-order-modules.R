cat(paste("\n Starting epist-order-modules at", date()))

## A minimal thing, to make sure no screw ups

## We no longer accept data frames. Those checks disabled

epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
## epistm1b <- data.frame(ids = c("a:d", "d:c"), s = c(0.2, 0.3))

oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)
## oeffects1b <- data.frame(ids = c("d>a", "c > d"),
##                          s = c(0.4, -0.3))

## oeffects1c <- data.frame(ids = c("a>d", "c > d"),
##                          s = c(0.4, -0.3))



## test_that("epist to long identical", {
##     expect_identical(to.long.epist.order(epistm1, ":"),
##                      to.long.epist.order(epistm1b, ":"))
## })

## test_that("order effects to long identical", {
##     expect_identical(to.long.epist.order(oeffects1, ">"),
##                      to.long.epist.order(oeffects1b, ">"))
## })

## test_that("order effects to long different", {
##     expect_false(identical(to.long.epist.order(oeffects1, ">"),
##                      to.long.epist.order(oeffects1c, ">")))
## })




gM <- c("Root" = "Root", "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")
gM2 <- c("Root" = "Root", "M" = "1", "B" = "3, 4", "A" = "6")
gM3 <- c("Root" = "Root", "b" = "1, 2", "a" = "3, 4, 5", "c" = "6")

## next two check ordered OK as they are matched against same output
gM4 <- c("Root" = "Root", "M" = "a , b", "B" = "uVw , vzt", "A" = "mXy")
gM5 <- c("Root" = "Root",  "A" = "mXy", "M" = "a , b", "B" = "uVw , vzt")




out1 <- structure(list(Gene = c("Root", "1", "2", "3", "4", "5", "6"), 
                       Module = c("Root", "a", "a", "b", "b", "b", "c"),
                       GeneNumID = 0:6, 
                       ModuleNumID = c(0L, 1L, 1L, 2L, 2L, 2L, 3L)),
                  .Names = c("Gene", 
                      "Module", "GeneNumID", "ModuleNumID"), row.names = c(NA, 7L),
                  class = "data.frame")

out2 <- structure(list(Gene = c("Root", "1", "3", "4", "6"),
                       Module = c("Root", 
                           "M", "B", "B", "A"),
                       GeneNumID = 0:4, ModuleNumID = c(0L, 1L, 
                                            2L, 2L, 3L)),
                  .Names = c("Gene", "Module", "GeneNumID", "ModuleNumID"
                             ), row.names = c(NA, 5L), class = "data.frame")

out3 <- structure(list(Gene = c("Root", "1", "2", "3", "4", "5", "6"), 
                       Module = c("Root", "b", "b", "a", "a", "a", "c"),
                       GeneNumID = 0:6, 
                       ModuleNumID = c(0L, 1L, 1L, 2L, 2L, 2L, 3L)),
                  .Names = c("Gene", 
                      "Module", "GeneNumID", "ModuleNumID"),
                  row.names = c(NA, 7L), class = "data.frame")

out4 <- structure(list(Gene = c("Root", "a", "b", "mXy", "uVw", "vzt"
), Module = c("Root", "M", "M", "A", "B", "B"), GeneNumID = 0:5, 
    ModuleNumID = c(0L, 1L, 1L, 3L, 2L, 2L)), .Names = c("Gene", 
"Module", "GeneNumID", "ModuleNumID"), row.names = c(1L, 2L, 
3L, 4L, 5L, 6L), class = "data.frame")


out5 <- structure(list(Gene = c("Root", "a", "b", "mXy", "uVw", "vzt"
), Module = c("Root", "M", "M", "A", "B", "B"), GeneNumID = 0:5, 
    ModuleNumID = c(0L, 2L, 2L, 1L, 3L, 3L)), .Names = c("Gene", 
"Module", "GeneNumID", "ModuleNumID"), row.names = c(1L, 2L, 
3L, 4L, 5L, 6L), class = "data.frame")


test_that("gm to geneModule, ex1", {
    expect_identical(OncoSimulR:::gm.to.geneModuleL(gM, FALSE),
                     out1)
})

test_that("gm to geneModule, ex2", {
    expect_identical(OncoSimulR:::gm.to.geneModuleL(gM2, FALSE),
                     out2)
})

test_that("gm to geneModule, ex3", {
    expect_identical(OncoSimulR:::gm.to.geneModuleL(gM3, FALSE),
                     out3)
})

test_that("gm to geneModule, ex4", {
    expect_identical(OncoSimulR:::gm.to.geneModuleL(gM4, FALSE),
                     out4)
})

test_that("gm to geneModule, ex5", {
    expect_identical(OncoSimulR:::gm.to.geneModuleL(gM5, FALSE),
                     out5)
})


gb <- c("a" = "1, 2", "b" = "3, 4, 5", "c" = "6")
gb1 <- c("a" = "1, 2", "Root" = "Root", "b" = "3, 4, 5", "c" = "6")

## this is no longer an error as we can deal with no Root
## test_that("error if no root in gm", {
##              expect_error(OncoSimulR:::gm.to.geneModuleL(gb, FALSE))
##           })


test_that("error if root out of place in gm", {
             expect_error(OncoSimulR:::gm.to.geneModuleL(gb1, FALSE))
          })


cat(paste("\n Ending epist-order-modules at", date()))
