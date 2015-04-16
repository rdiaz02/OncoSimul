## A set of tests to verify transformations between graph formats are
## OK. Recall that restrictionTable acts like a "sink": nothing is
## transformed from a restctiionTable to anything else.


## For really exhaustive, set this to a big number
## numTests <- 5000
numTests <- 100

## FIXME: still need graph.to.poset and recheck others

createAndConvert <- function(rangeNodes = 4:30,
                            rangeParents = 2:5,
                             verbose = FALSE) {
    tp <- within(list( 
        nodes = sample(rangeNodes, 1),
        parents = sample(rangeParents, 1)
    ),
                 h <- sample(nodes, 1)
                 )

    if(verbose)
        print(tp)

    g1 <- OncoSimulR:::simOGraph(n = tp$nodes,
                    h = tp$h,
                    nparents = tp$parents,
                    conjunction = TRUE,
                    multilevelParent = TRUE,
                    removeIndirect = TRUE,
                    rootName = "0"
                    )

    g1.To.p1 <- OncoSimulR:::intAdjMatToPoset(g1)
    g1.To.p1.To.g1 <- OncoSimulR:::posetToGraph(g1.To.p1,
                                                names = 0:tp$nodes,
                                                addroot = TRUE,
                                                type = "adjmat")


    g1.To.p1.To.g1.No.Root <- OncoSimulR:::posetToGraph(g1.To.p1,
                                                        names = 1:tp$nodes,
                                                        addroot = FALSE,
                                                        type = "adjmat")

    
    g1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(g1, root = TRUE)
    g1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(g1[-1, -1], root = FALSE)
    g1.To.p1.To.rt <- OncoSimulR:::poset.to.restrictTable(g1.To.p1)

    return(list(g1 = g1,
                g1.To.p1 = g1.To.p1,
                g1.To.p1.To.g1 = g1.To.p1.To.g1,
                g1.To.p1.To.g1.No.Root = g1.To.p1.To.g1.No.Root,
                g1.To.rt = g1.To.rt,
                g1.To.rt.2 = g1.To.rt.2,
                g1.To.p1.To.rt = g1.To.p1.To.rt
                )
       )
}



masterTestCall <- function(rangeNodes = 4:30,
                           rangeParents = 2:5,
                           verbose = FALSE) {

    gg <- createAndConvert(rangeNodes = rangeNodes,
                           rangeParents = rangeParents,
                           verbose = verbose)


    test_that("graph identical to graph->poset->graph", {
                  expect_identical(gg$g1, gg$g1.To.p1.To.g1)
           })
    test_that("graph identical to graph->poset->graph, without Root", {
        expect_identical(gg$g1[-1, -1], gg$g1.To.p1.To.g1.No.Root)
    })
    test_that("restriction table: identical graph->rT and graph->poset->rT", {
        expect_identical(gg$g1.To.rt, gg$g1.To.p1.To.rt)
    })
    test_that("restriction table: identical graph->rT and graph->poset->rT, No root", {
        expect_identical(gg$g1.To.rt.2, gg$g1.To.p1.To.rt)
    })
    ## cat("A full round of tests completed OK.\n")
    return("OK")
}

tmp <- replicate(numTests, masterTestCall() )




## testGraphTransf <- function(rangeNodes = 4:30,
##                             rangeParents = 2:5,
##                             verbose = FALSE) {
##     ## tp <- within(list( 
##     ##     nodes = sample(rangeNodes, 1),
##     ##     parents = sample(rangeParents, 1)
##     ## ),
##     ##              h <- sample(nodes, 1)
##     ##              )

##     ## if(verbose)
##     ##     print(tp)

##     ## g1 <- simOGraph(n = tp$nodes,
##     ##                 h = tp$h,
##     ##                 nparents = tp$parents,
##     ##                 conjunction = TRUE,
##     ##                 multilevelParent = TRUE,
##     ##                 removeIndirect = TRUE,
##     ##                 rootName = "0"
##     ##                 )

##     ## g1.To.p1 <- intAdjMatToPoset(g1)
##     ## g1.To.p1.To.g1 <- OncoSimulR:::posetToGraph(g1.To.p1,
##     ##                                             names = 0:tp$nodes,
##     ##                                             addroot = TRUE,
##     ##                                             type = "adjmat")


##     ## g1.To.p1.To.g1.No.Root <- OncoSimulR:::posetToGraph(g1.To.p1,
##     ##                                                     names = 1:tp$nodes,
##     ##                                                     addroot = FALSE,
##     ##                                                     type = "adjmat")

    
##     ## g1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(g1, root = TRUE)
##     ## g1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(g1[-1, -1], root = FALSE)
##     ## p1.To.rt <- OncoSimulR:::poset.to.restrictTable(g1.To.p1)

##     stopifnot(identical(g1, g1.To.p1.To.g1))
##     stopifnot(identical(g1[-1, -1], g1.To.p1.To.g1.No.Root))
##     stopifnot(identical(g1.To.rt, p1.To.rt))
##     stopifnot(identical(g1.To.rt.2, p1.To.rt))
##     return(TRUE)
## }

## resT <- replicate(numTests, testGraphTransf() )

## stopifnot(length(resT) == numTests)
## stopifnot(all(resT))

