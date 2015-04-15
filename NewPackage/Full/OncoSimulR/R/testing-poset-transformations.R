## A set of tests to verify transformations between graph formats are
## OK. Recall that restrictionTable acts like a "sink": nothing is
## transformed from a restctiionTable to anything else.


source("generate-random-trees.R")
library(OncoSimulR)

numTests <- 5000

testGraphTransf <- function(rangeNodes = 4:30,
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

    g1 <- simOGraph(n = tp$nodes,
                    h = tp$h,
                    nparents = tp$parents,
                    conjunction = TRUE,
                    multilevelParent = TRUE,
                    removeIndirect = TRUE,
                    rootName = "0"
                    )

    g1.To.p1 <- intAdjMatToPoset(g1)
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
    p1.To.rt <- OncoSimulR:::poset.to.restrictTable(g1.To.p1)

    stopifnot(identical(g1, g1.To.p1.To.g1))
    stopifnot(identical(g1[-1, -1], g1.To.p1.To.g1.No.Root))
    stopifnot(identical(g1.To.rt, p1.To.rt))
    stopifnot(identical(g1.To.rt.2, p1.To.rt))
    return(TRUE)
}

resT <- replicate(numTests, testGraphTransf() )

stopifnot(length(resT) == numTests)
stopifnot(all(resT))
