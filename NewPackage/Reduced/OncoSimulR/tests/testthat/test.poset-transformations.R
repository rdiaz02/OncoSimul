## A preliminary set of tests. More in the current unreleased code.
## library(OncoSimulR); library(testthat)





## This code is in devel
## m0 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m0) <- rownames(m0) <- c(0, 2, 3, 5)
## m0[1, 4] <- 1L

## test_that("adjmat illegal even in this conversion",
##           expect_error(OncoSimulR:::OTtoPoset(m0)))


## m0 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m0) <- rownames(m0) <- c(0, 2, 3, 5)
## m0[1,  ] <- 1L
## test_that("does not conform to having root be called Root",
##           expect_error(OncoSimulR:::OTtoPoset(m0)))

## m0 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m0) <- rownames(m0) <- c("Root", 2, 3, 5)
## m0[1, ] <- 1L
## test_that("but an edge to Root",
##           expect_error(OncoSimulR:::OTtoPoset(m0)))

## m0 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m0) <- rownames(m0) <- c("Root", 2, 3, 5)
## m0[1, 2:4] <- 1L
## test_that("OT to the smallest Poset",
##           expect_equal(matrix(nrow=0, ncol=2), OncoSimulR:::OTtoPoset(m0)))

## new devel code
## test_that("adjmat illegal because empty nodes",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m0, root = TRUE)))

## test_that("adjmat illegal because empty nodes, 2",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m0, root = FALSE)))

## m1 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m1) <- rownames(m1) <- c(0, 2, 3, 5)
## m1[, 4] <- 1L

## test_that("adjmat illegal because lack of incoming",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = TRUE)))
## test_that("adjmat illegal because lack of incoming, 2",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = FALSE)))
## test_that("adjmat illegal because non-int column names",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1[-1, -1], root = FALSE)))
## test_that("adjmat illegal because of no root",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1[-1, -1], root = TRUE)))

## m2 <- matrix(0L, ncol = 4, nrow = 4)
## colnames(m2) <- rownames(m2) <- c(0:3)
## m2[, 4] <- 1L; m2[1, 3] <- 1L

## test_that("adjmat illegal because lack of incoming",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m2, root = TRUE)))
## test_that("adjmat illegal because lack of incoming, 2",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m2, root = FALSE)))


m3 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m3) <- rownames(m3) <- c(0:3)
m3[1, 2:4] <- 1L

## test_that("inconsistency in root options",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = FALSE)))

rt3 <- cbind(1:3, 0L); class(rt3) <- "restrictionTable"
test_that("simple correct adjmat -> rT", {
    expect_equal(rt3, OncoSimulR:::adjmat.to.restrictTable(m3, root = TRUE))})


p3 <- cbind(0L, 3L)
test_that("simple correct poset -> graph", {
    expect_equal(m3,
                 OncoSimulR:::posetToGraph(p3, names = 0:3, addroot = TRUE, type = "adjmat"))})



m3 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m3) <- rownames(m3) <- c(0:3)
m3[1, ] <- 1L
## test_that("adjmat illegal because incoming to Root",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = TRUE)))
## test_that("adjmat illegal because incoming to Root",
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = FALSE)))




p5 <- cbind(c(0L), c(5L))
m5 <- matrix(0L, ncol = 6, nrow = 6); colnames(m5) <- rownames(m5) <- 0:5
m5[1, 2:6] <- 1L
test_that("simple correct poset -> graph, 5 nodes, root", {
    expect_equal(m5,
                 OncoSimulR:::posetToGraph(p5, names = 0:5,
                                           addroot = TRUE, type = "adjmat"))})
test_that("simple correct poset -> graph, 5 nodes, no root", {
    expect_equal(m5[-1, -1],
                 OncoSimulR:::posetToGraph(p5, names = 1:5,
                                           addroot = FALSE, type = "adjmat"))})


test_that("to rT from poset, through adjmat with and w.o. root",
          {
              expect_equal(
                  OncoSimulR:::adjmat.to.restrictTable(
                      OncoSimulR:::posetToGraph(p5, names = 1:5, addroot = FALSE, type = "adjmat"),
                      root = FALSE),
                  OncoSimulR:::adjmat.to.restrictTable(
                      OncoSimulR:::posetToGraph(p5, names = 0:5, addroot = TRUE, type = "adjmat"),
                      root = TRUE))
          })




## can convert if OT, not o.w.
## m1 <- structure(c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L), .Dim = c(3L, 
## 3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
##                                             )))
## test_that("from adjmat to rT, if OT", {
##           expect_identical(cbind(2L, 4L), OncoSimulR:::OTtoPoset(m1))})
## test_that(".. but fails if not considered OT", {
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = TRUE))})
## test_that(".. but fails if not considered OT", {
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = FALSE))})



## m1b <- structure(c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L), .Dim = c(3L, 
## 3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
## )))
## test_that("Change pos of elements: from adjmat to rT, if OT", {
##           expect_identical(cbind(4L, 2L), OncoSimulR:::OTtoPoset(m1b))})
## test_that(".. but fails if not considered OT", {
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1b, root = TRUE))})
## test_that(".. but fails if not considered OT", {
##           expect_error(OncoSimulR:::adjmat.to.restrictTable(m1b, root = FALSE))})



## posets with some nodes not explicit
pm1 <- structure(c(0L, 2L, 2L, 4L), .Dim = c(2L, 2L))
pm1.nr <- structure(c(2L, 4L), .Dim = 1:2)
test_that("poset to rT, with some nodes missing",{
          expect_identical(OncoSimulR:::poset.to.restrictTable(pm1.nr),
                           OncoSimulR:::poset.to.restrictTable(pm1))})


pm1b <- structure(c(4L, 0L, 2L, 4L), .Dim = c(2L, 2L))
pm1b.nr <- structure(c(4L, 2L), .Dim = 1:2)
test_that("poset to rT, with some nodes missing, nodes exchanged",{
          expect_identical(OncoSimulR:::poset.to.restrictTable(pm1b.nr),
                           OncoSimulR:::poset.to.restrictTable(pm1b))})





