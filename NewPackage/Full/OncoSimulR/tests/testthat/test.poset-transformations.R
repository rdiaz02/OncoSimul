## posetToAdjMat and posetToGraph: always keeproot.
## so modify plotPoset and see what run.oncotree does.


## A set of tests to verify transformations between graph formats are
## OK. Recall that restrictionTable acts like a "sink": nothing is
## transformed from a restctiionTable to anything else.

## Recall that restrcitionTable only uses integers for now.

## We have, first, some specific tests and then a general testing call
## with random trees. The second might not catch specific issues, such as
## not all nodes being in the poset, etc.


m1 <- structure(c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L), .Dim = c(3L, 
3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
                                             )))

m1b <- structure(c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L), .Dim = c(3L, 
3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
)))


pm1 <- structure(c(0L, 2L, 2L, 4L), .Dim = c(2L, 2L))
pm1.nr <- structure(c(2L, 4L), .Dim = 1:2)

pm1b <- structure(c(4L, 0L, 2L, 4L), .Dim = c(2L, 2L))
pm1b.nr <- structure(c(4L, 2L), .Dim = 1:2)

test_that("adjmat to Poset, example 1", {
    expect_identical(OncoSimulR:::adjMatToPoset(m1, dropRoot = FALSE),
                     pm1)
})

test_that("adjmat to Poset, example 1, no root", {
    expect_identical(OncoSimulR:::adjMatToPoset(m1, dropRoot = TRUE),
                     pm1.nr)
})

test_that("adjmat to Poset, example 2", {
    expect_identical(OncoSimulR:::adjMatToPoset(m1b, dropRoot = FALSE),
                     pm1b)
})

test_that("adjmat to Poset, example 2, no root", {
    expect_identical(OncoSimulR:::adjMatToPoset(m1b, dropRoot = TRUE),
                     pm1b.nr)
})

## adjmat -> poset -> adjMat



## For really exhaustive, set this to a big number
## numTests <- 5000
numTests <- 100

## FIXME: still need graph.to.poset and recheck others

createAndConvert <- function(rangeNodes = 4:30,
                             rangeParents = 2:5,
                             verbose = TRUE,
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

    am11 <<- am1
    
    ## FIXME: for now, we only excercise names = 0.
    ## Posets only deal with integer codes, and same for restrictTable.
    
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

    
    
    am1.nr <- am1[-1, -1]
    gf1 <- as(am1, "graphNEL")
    p1 <- OncoSimulR:::adjMatToPoset(am1) ## am1.To.p1

    p1.nr <- OncoSimulR:::adjMatToPoset(am1, dropRoot = TRUE)
    
    p1.am1 <- OncoSimulR:::posetToGraph(p1,
                                        names = 0:tp$nodes,
                                        addroot = TRUE,
                                        type = "adjmat")

    p1.am1.p1 <- OncoSimulR:::adjMatToPoset(p1.am1)

    
    p1.am1.nr <- OncoSimulR:::posetToGraph(p1,
                                           names = 1:tp$nodes,
                                           addroot = FALSE,
                                           type = "adjmat")
    

    p1.nr.am1 <- OncoSimulR:::posetToGraph(p1.nr,
                                           names = 0:tp$nodes,
                                           addroot = TRUE,
                                           type = "adjmat")

    p1.nr.am1.nr <- OncoSimulR:::posetToGraph(p1.nr,
                                              names = 1:tp$nodes,
                                              addroot = FALSE,
                                              type = "adjmat")


    am1.gf1.am1 <- as(as(am1, "graphNEL"), "matrix")
    gf1.am1.gf1 <- as(as(gf1, "matrix"), "graphNEL")

    gf1.nr <-  as(am1[-1, -1], "graphNEL")
    
    NR.am1.gf1.am1 <- as(as(am1.nr, "graphNEL"), "matrix")
    NR.gf1.am1.gf1 <- as(as(gf1.nr, "matrix"), "graphNEL")

    
    gf1.p1 <- OncoSimulR:::graphToPoset(gf1)
    p1.gf1 <- OncoSimulR:::posetToGraph(p1,
                           names = 0:tp$nodes,
                           addroot = TRUE,
                           type = "graphNEL")
    p1.gf1.p1 <- OncoSimulR:::graphToPoset(p1.gf1)

    p1.gf1.nr <- OncoSimulR:::posetToGraph(p1,
                              names = 1:tp$nodes,
                              addroot = FALSE,
                              type = "graphNEL")

    ## these cannot be done now, as adj matrices without 0 are disabled
    ## gf1.nr.p1 <- OncoSimulR:::graphToPoset(gf1.nr)
    ## gf1.nr.p1.B <- OncoSimulR:::graphToPoset(p1.gf1.nr)

    p1.nr.gf1 <- OncoSimulR:::posetToGraph(p1.nr,
                              names = 0:tp$nodes,
                              addroot = TRUE,
                              type = "graphNEL")
    p1.nr.gf1.nr <- OncoSimulR:::posetToGraph(p1.nr,
                              names = 1:tp$nodes,
                              addroot = FALSE,
                              type = "graphNEL")


    am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1, root = TRUE)
    am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)
    p1.To.rt <- OncoSimulR:::poset.to.restrictTable(p1)
    p1.To.rt.2 <- OncoSimulR:::poset.to.restrictTable(p1.nr)
   
   
    ## p1.gf1.p1 <- graphToPoset(posetTograph(p1, keeproot = TRUE))
    ## gf1.p1.gf1 <- posetToGraph(graphToPoset(gf1), keeproot = TRUE)

    

    
    
    ## am1.To.p1.To.am1 <- OncoSimulR:::posetToGraph(am1.To.p1,
    ##                                               keeproot = TRUE,
    ##                                               type = "adjmat")

    ## posets must have root
    ## am1.To.p1.To.am1.No.Root <- OncoSimulR:::posetToGraph(am1.To.p1,
    ##                                                       keeproot = FALSE,
    ##                                                       type = "adjmat")

    
    ## if(names == 0) {
    ##     am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1)
    ##     ## am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)
    ##     p1.To.rt <- OncoSimulR:::poset.to.restrictTable(p1)
    ## } else {
    ##     am1.To.rt <- NA
    ##     ## am1.To.rt.2 <- NA
    ##     p1.To.rt <- NA
    ## }

    return(list(
        am1 = am1,
        am1.nr = am1.nr,
        gf1 = gf1,
        p1 = p1,
        p1.nr = p1.nr,
        p1.am1 = p1.am1,
        p1.am1.p1 = p1.am1.p1,
        p1.am1.nr = p1.am1.nr,
        p1.nr.am1 = p1.nr.am1,
        p1.nr.am1.nr = p1.nr.am1.nr,
        am1.gf1.am1 = am1.gf1.am1,
        gf1.am1.gf1 = gf1.am1.gf1,
        gf1.nr = gf1.nr,
        NR.am1.gf1.am1 = NR.am1.gf1.am1,
        NR.gf1.am1.gf1 = NR.gf1.am1.gf1,
        gf1.p1 = gf1.p1,
        p1.gf1 = p1.gf1,
        p1.gf1.p1 = p1.gf1.p1,
        p1.gf1.nr = p1.gf1.nr,
        ## gf1.nr.p1 = gf1.nr.p1,
        ## gf1.nr.p1.B = gf1.nr.p1.B,
        p1.nr.gf1 = p1.nr.gf1,
        p1.nr.gf1.nr = p1.nr.gf1.nr,
        am1.To.rt   = am1.To.rt, 
        am1.To.rt.2 = am1.To.rt.2, 
        p1.To.rt    = p1.To.rt, 
        p1.To.rt.2  = p1.To.rt.2        
    ))  
    ## return(list(am1 = am1,
    ##             gf1 = gf1,
    ##             p1 = p1,
    ##             p1.am1 = p1.am1,
    ##             p1.am1.p1 = p1.am1.p1,
    ##             am1.gf1.am1 = am1.gf1.am1,
    ##             gf1.am1.gf1 = gf1.am1.gf1,
    ##             p1.gf1.p1 = p1.gf1.p1,
    ##             gf1.p1.gf1 = gf1.p1.gf1,
    ##             am1.nr = am1.nr,
    ##             gf1.nr = gf1.nr,
    ##             NR.am1.gf1.am1 = NR.am1.gf1.am1,
    ##             NR.gf1.am1.gf1 = NR.gf1.am1.gf1))
}





## FIXME
## For now, any conversion to rT is from things with integer names
## between other things with arbitrary names??




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
                           verbose = TRUE,
                           names = 0) {

    out <- createAndConvert(rangeNodes = rangeNodes,
                           rangeParents = rangeParents,
                            verbose = verbose,
                            names = names)


    test_that("adjmat->poset->adjmat", {
        expect_identical(out$am1, out$p1.am1)
    })
    ## the next one is currently redundant
    test_that("poset->adjmat->poset", {
        expect_identical(out$p1, out$p1.am1.p1)
    })

    test_that("adjmat->poset->adjmat, no root", {
        expect_identical(out$am1.nr, out$p1.am1.nr)
    })


    test_that("sizes return adj mats w/w.o nr", {
        expect_true( nrow(out$p1.am1) == (nrow(out$p1.am1.nr) + 1) ) 
    })

    test_that("different posets, p1, p1.nr", {
        expect_false(identical(out$p1 , out$p1.nr))
    })

    test_that("p1 larger than p1.nr", {
        expect_true( nrow(out$p1) > nrow(out$p1.nr) )
    })
   
    
    test_that("adjmat->poset.nr ->adjmat", {
        expect_identical(out$am1, out$p1.nr.am1)
    })

    test_that("adjmat->poset.nr ->adjmat, no root", {
        expect_identical(out$am1.nr, out$p1.nr.am1.nr)
    })

    ## next two are redundant
    test_that("adjmat->graph->adjmat", {
        ## equality, not identical, since storage.mode differs
        expect_equal(out$am1, out$am1.gf1.am1)
    })
    test_that("graph->adjmat->graph", {
        expect_equal(out$gf1, out$gf1.am1.gf1)
    })
    
    ## next two are redundant
    test_that("adjmat->graph->adjmat, no root", {
        expect_equal(out$am1.nr, out$NR.am1.gf1.am1)
    })
    test_that("graph->adjmat->graph, no root", {
        expect_equal(out$gf1.nr, out$NR.gf1.am1.gf1)
    })


    test_that("poset = (graph->poset)", {
        expect_identical(out$p1, out$gf1.p1)
    })
    
    test_that("graph =   (poset -> graph) ", {
        expect_identical(out$gf1, out$p1.gf1)
    })

    ## redundant
    test_that("poset -> graph->poset ", {
        expect_identical(out$p1, out$p1.gf1.p1)
    })

    test_that("poset -> graph.nr  = gf1.nr ", {
        expect_identical(out$gf1.nr, out$p1.gf1.nr)
    })

    ## No adjanceny matrices without 0 allowed
    ## test_that("graph.nr -> poset  =  am.nr -> poset ", {
    ##     expect_identical(out$gf1.nr.p1, out$p1.nr)
    ## })

    ## test_that("graph.nr -> poset =  poset-> graf.nr -> poset ", {
    ##     expect_identical(out$gf1.nr.p1, out$gf1.nr.p1.B)
    ## })

    test_that("gf1 =  poset.nr ->graph", {
        expect_identical(out$gf1, out$p1.nr.gf1)
    })

    test_that("gf1.nr =  poset.nr ->graph nr", {
        expect_identical(out$gf1.nr, out$p1.nr.gf1.nr)
    })

    test_that(" from p1.nr to gf1 w and w/o root differ", {
        expect_false(identical(out$p1.nr.gf1, out$p1.nr.gf1.nr))
    })

    ## test_that("graph->poset->graph", {
    ##     expect_identical(out$gf1, out$gf1.p1.gf1)
    ## })
  
    
    ## if(names == 0) {
        test_that("restriction table: identical from adjmat and poset", {
            expect_identical(out$am1.To.rt, out$p1.To.rt)
        })
        test_that("restriction table: identical from adjmat w/wo root", {
            expect_identical(out$am1.To.rt, out$am1.To.rt.2)
        })
        test_that("restriction table: identical from poset w/wo root", {
            expect_identical(out$p1.To.rt, out$p1.To.rt.2)
        })
  ##   }
    ## cat("A full round of tests completed OK.\n")
    return("OK")
}







tmp <- replicate(numTests, masterTestCall() )


## for(nn in 0:4) 
##     tmp <- replicate(numTests, masterTestCall(names = nn) )






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

