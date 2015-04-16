## A set of tests to verify transformations between graph formats are
## OK. Recall that restrictionTable acts like a "sink": nothing is
## transformed from a restctiionTable to anything else.


## For really exhaustive, set this to a big number
## numTests <- 5000
numTests <- 100

## FIXME: still need graph.to.poset and recheck others

createAndConvert <- function(rangeNodes = 4:30,
                             rangeParents = 2:5,
                             verbose = FALSE,
                             names = 0) {
    tp <- within(list( 
        nodes = sample(rangeNodes, 1),
        parents = sample(rangeParents, 1)
    ),
                 h <- sample(nodes, 1)
                 )

    if(verbose)
        print(tp)

    
    am1 <- OncoSimulR:::simOGraph(n = tp$nodes,
                    h = tp$h,
                    nparents = tp$parents,
                    conjunction = TRUE,
                    multilevelParent = TRUE,
                    removeIndirect = TRUE,
                    rootName = "0"
                    )


    if(names == 1) {
        rn <- replicate(nrow(am1) - 1,
                        paste(sample(letters, 10, replace = TRUE),
                              collapse = ""))
        rn <- c("Root", rn)
    } else if(names == 2) {
        rn <- replicate(nrow(am1) - 1,
                        paste(sample(letters, 10, replace = TRUE),
                              collapse = ""))
        rn <- c("0", rn)
    } else if(names == 3) {
        rn <- paste0("G", seq.int(nrow(am1) - 1))
        rn <- c("Root", rn)
    } else if(names == 4) {
        rn <- paste0("G", seq.int(nrow(am1) - 1))
        rn <- c("0", rn)
    }
    if(names != 0)
        rownames(am1) <- colnames(am1) <- rn
    
    am1.To.p1 <- OncoSimulR:::intAdjMatToPoset(am1)
    am1.To.p1.To.am1 <- OncoSimulR:::posetToGraph(am1.To.p1,
                                                  names = NULL,
                                                  ## names = 0:tp$nodes,
                                                  addroot = TRUE,
                                                  type = "adjmat")


    am1.To.p1.To.am1.No.Root <- OncoSimulR:::posetToGraph(am1.To.p1,
                                                          names = NULL,
                                                          ## names = 1:tp$nodes,
                                                          addroot = FALSE,
                                                          type = "adjmat")

    
    am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1, root = TRUE)
    am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)
    am1.To.p1.To.rt <- OncoSimulR:::poset.to.restrictTable(am1.To.p1)


    gf1 <- as(am1, "graphNEL")

    gf1.To.p1 <- OncoSimulR:::graphToPoset(gf1)

    return(list(am1 = am1,
                am1.To.p1 = am1.To.p1,
                am1.To.p1.To.am1 = am1.To.p1.To.am1,
                am1.To.p1.To.am1.No.Root = am1.To.p1.To.am1.No.Root,
                am1.To.rt = am1.To.rt,
                am1.To.rt.2 = am1.To.rt.2,
                am1.To.p1.To.rt = am1.To.p1.To.rt
                )
       )
}






## posetToGraph as graphnel output

## am -> graphNel  and graphNel -> am
## graphNel -> poset and poset -> graphNEL


##  am <-> poset
##  am  <-> graphNel
##  graphNel <-> poset

## Nice full circle: am -> graphNel -> poset -> am
## am -> poset -> graphNEL
## etc, starting at each point


masterTestCall <- function(rangeNodes = 4:30,
                           rangeParents = 2:5,
                           verbose = FALSE,
                           names = 0) {

    out <- createAndConvert(rangeNodes = rangeNodes,
                           rangeParents = rangeParents,
                            verbose = verbose,
                            names = names)


    test_that("adjmat identical to adjmat->poset->adjmat", {
                  expect_identical(out$am1, out$am1.To.p1.To.am1)
           })
    test_that("adjmat identical to adjmat->poset->adjmat, without Root", {
        expect_identical(out$am1[-1, -1], out$am1.To.p1.To.am1.No.Root)
    })
    test_that("restriction table: identical adjmat->rT and adjmat->poset->rT", {
        expect_identical(out$am1.To.rt, out$am1.To.p1.To.rt)
    })
    test_that("restriction table: identical adjmat->rT and adjmat->poset->rT, No root", {
        expect_identical(out$am1.To.rt.2, out$am1.To.p1.To.rt)
    })
    ## cat("A full round of tests completed OK.\n")
    return("OK")
}

for(nn in 0:4) 
    tmp <- replicate(numTests, masterTestCall(names = nn) )




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

