## A set of tests to verify transformations between graph formats are
## OK. Recall that restrictionTable acts like a "sink": nothing is
## transformed from a restctiionTable to anything else.

## We have, first, some specific tests and then a general testing call
## with random trees. The second might not catch specific issues, such as
## not all nodes being in the poset, etc.


## include code that tests convertRestrictTable. Nope, not now.

## Verify

## poset to rT. From the different kinds of posets?

## all the code for new simulations is based upon going from adjmat to rT.
## verify that.

cat(paste("\n Starting poset-transformations tests", date(), "\n"))

## RNGkind("Mersenne-Twister")

test_that("posetToGraph stop in incorrect entry type", {
    expect_error(OncoSimulR:::posetToGraph(1:5, letters[1:5]),
                 "If poset is not a matrix, it must be a vector of length 1")
    expect_error(OncoSimulR:::posetToGraph(1:5, letters[1:2]),
                 "If poset is not a matrix, it must be a vector of length 1")
    expect_error(OncoSimulR:::posetToGraph(data.frame(a = 1:5, b = 1:5),
                                           letters[1:5]),
                 "poset contains non-integers")
})



m0 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m0) <- rownames(m0) <- c(0, 2, 3, 5)
m0[1, 4] <- 1L

test_that("adjmat illegal even in this conversion",
          expect_error(OncoSimulR:::OTtoPoset(m0)))


m0 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m0) <- rownames(m0) <- c(0, 2, 3, 5)
m0[1,  ] <- 1L
test_that("does not conform to having root be called Root",
          expect_error(OncoSimulR:::OTtoPoset(m0)))

m0 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m0) <- rownames(m0) <- c("Root", 2, 3, 5)
m0[1, ] <- 1L
test_that("but an edge to Root",
          expect_error(OncoSimulR:::OTtoPoset(m0)))

m0 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m0) <- rownames(m0) <- c("Root", 2, 3, 5)
m0[1, 2:4] <- 1L
test_that("OT to the smallest Poset",
          expect_equal(matrix(nrow=0, ncol=2), OncoSimulR:::OTtoPoset(m0)))


test_that("adjmat illegal because empty nodes",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m0, root = TRUE)))

test_that("adjmat illegal because empty nodes, 2",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m0, root = FALSE)))

m1 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m1) <- rownames(m1) <- c(0, 2, 3, 5)
m1[, 4] <- 1L

test_that("adjmat illegal because lack of incoming",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = TRUE)))
test_that("adjmat illegal because lack of incoming, 2",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = FALSE)))
test_that("adjmat illegal because non-int column names",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1[-1, -1], root = FALSE)))
test_that("adjmat illegal because of no root",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1[-1, -1], root = TRUE)))

m2 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m2) <- rownames(m2) <- c(0:3)
m2[, 4] <- 1L; m2[1, 3] <- 1L

test_that("adjmat illegal because lack of incoming",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m2, root = TRUE)))
test_that("adjmat illegal because lack of incoming, 2",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m2, root = FALSE)))


m3 <- matrix(0L, ncol = 4, nrow = 4)
colnames(m3) <- rownames(m3) <- c(0:3)
m3[1, 2:4] <- 1L

test_that("inconsistency in root options",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = FALSE)))

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
test_that("adjmat illegal because incoming to Root",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = TRUE)))
test_that("adjmat illegal because incoming to Root",
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m3, root = FALSE)))




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
m1 <- structure(c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L), .Dim = c(3L, 
3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
                                             )))
test_that("from adjmat to rT, if OT", {
          expect_identical(cbind(2L, 4L), OncoSimulR:::OTtoPoset(m1))})
test_that(".. but fails if not considered OT", {
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = TRUE))})
test_that(".. but fails if not considered OT", {
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1, root = FALSE))})



m1b <- structure(c(0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L), .Dim = c(3L, 
3L), .Dimnames = list(c("Root", "2", "4"), c("Root", "2", "4"
)))
test_that("Change pos of elements: from adjmat to rT, if OT", {
          expect_identical(cbind(4L, 2L), OncoSimulR:::OTtoPoset(m1b))})
test_that(".. but fails if not considered OT", {
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1b, root = TRUE))})
test_that(".. but fails if not considered OT", {
          expect_error(OncoSimulR:::adjmat.to.restrictTable(m1b, root = FALSE))})



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





numTests <- 40 ## in long, we use 100

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

    
    am1 <- OncoSimulR::simOGraph(n = tp$nodes,
                    h = tp$h,
                    nparents = tp$parents,
                    conjunction = TRUE,
                    multilevelParent = TRUE,
                    removeDirectIndirect = TRUE,
                    rootName = "0"
                    )
    
    am1.nr <- am1[-1, -1]
    gf1 <- as(am1, "graphNEL")

    am1.gf1.am1 <- as(as(am1, "graphNEL"), "matrix")
    gf1.am1.gf1 <- as(as(gf1, "matrix"), "graphNEL")

    gf1.nr <-  as(am1[-1, -1], "graphNEL")
    
    NR.am1.gf1.am1 <- as(as(am1.nr, "graphNEL"), "matrix")
    NR.gf1.am1.gf1 <- as(as(gf1.nr, "matrix"), "graphNEL")

    am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1, root = TRUE)
    am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)

    return(list(
        am1 = am1,
        am1.nr = am1.nr,
        gf1 = gf1,
        am1.gf1.am1 = am1.gf1.am1,
        gf1.am1.gf1 = gf1.am1.gf1,
        gf1.nr = gf1.nr,
        NR.am1.gf1.am1 = NR.am1.gf1.am1,
        NR.gf1.am1.gf1 = NR.gf1.am1.gf1,
        am1.To.rt   = am1.To.rt, 
        am1.To.rt.2 = am1.To.rt.2 
    ))  
}


masterTestCall <- function(rangeNodes = 4:30,
                           rangeParents = 2:5,
                           verbose = FALSE) {

    out <- createAndConvert(rangeNodes = rangeNodes,
                           rangeParents = rangeParents,
                            verbose = verbose)

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
    
    
    test_that("restriction table: identical from adjmat w/wo root", {
        expect_identical(out$am1.To.rt, out$am1.To.rt.2)
    })
    return("OK")
} 

tmp <- replicate(numTests, masterTestCall() )



## verify the simulator does generate correct adjacency matrices
numSimul <- 50  ## in long, we use a 100

checkAdjMatOGraph <- function(rangeNodes = 4:30,
                              rangeParents = 2:5,
                              multilevelP = TRUE,
                              verbose = FALSE) {
    tp <- within(list( 
        nodes = sample(rangeNodes, 1),
        parents = sample(rangeParents, 1)
    ),
                 h <- sample(nodes, 1)
                 )
    
    if(verbose)
        print(tp)
    
    
    am1 <- simOGraph(n = tp$nodes,
                                  h = tp$h,
                    nparents = tp$parents,
                    conjunction = TRUE,
                    multilevelParent = multilevelP,
                    removeDirectIndirect = TRUE,
                    rootName = "0"
                                  )
    am1 <<- am1
    am1nr <- am1[-1]
    am1nrnr <- am1[-1, -1]
    cmis <- sample(1:(ncol(am1) -1), 1) ## if last row and colum, it'd be ok
    am1mis <- am1[-cmis, -cmis]

    am1mis <<- am1mis
    test_that("original OK, output",
              expect_true(inherits(OncoSimulR:::adjmat.to.restrictTable(am1, root = TRUE),
                               "restrictionTable")))
    test_that("original OK",
              expect_identical(OncoSimulR:::checkProperFullAdjMat(am1), NULL)
              )
    test_that("error from no root",
              expect_error(OncoSimulR:::checkProperFullAdjMat(am1nrnr))
              )
    test_that("error from non-square",
              expect_error(OncoSimulR:::checkProperFullAdjMat(am1nr))
              )
    test_that("error from missing column and row",
              expect_error(OncoSimulR:::checkProperFullAdjMat(am1mis))
              )
}

tmp <- replicate(numSimul, checkAdjMatOGraph())

## here set parents to 1
tmp <- replicate(numSimul, checkAdjMatOGraph(rangeNodes = 2:50,
                                             rangeParents = 1))

## no multilevel parents
tmp <- replicate(numSimul, checkAdjMatOGraph(rangeNodes = 2:50,
                                             rangeParents = 1:5,
                                             multilevelP = FALSE))





####### Lots of old code, with a former set of much more comprehensive
####### conversions that are never used, and were thus removed.

## test_that("adjmat to Poset, example 1", {
##     expect_identical(OncoSimulR:::adjMatToPoset(m1, dropRoot = FALSE),
##                      pm1)
## })

## test_that("adjmat to Poset, example 1, no root", {
##     expect_identical(OncoSimulR:::adjMatToPoset(m1, dropRoot = TRUE),
##                      pm1.nr)
## })

## test_that("adjmat to Poset, example 2", {
##     expect_identical(OncoSimulR:::adjMatToPoset(m1b, dropRoot = FALSE),
##                      pm1b)
## })

## test_that("adjmat to Poset, example 2, no root", {
##     expect_identical(OncoSimulR:::adjMatToPoset(m1b, dropRoot = TRUE),
##                      pm1b.nr)
## })

## adjmat -> poset -> adjMat



## For really exhaustive, set this to a big number
## numTests <- 5000

## createAndConvert <- function(rangeNodes = 4:30,
##                              rangeParents = 2:5,
##                              verbose = TRUE,
##                              names = 0) {
##     tp <- within(list( 
##         nodes = sample(rangeNodes, 1),
##         parents = sample(rangeParents, 1)
##     ),
##                  h <- sample(nodes, 1)
##                  )

##     if(verbose)
##         print(tp)

    
##     am1 <- OncoSimulR:::simOGraph(n = tp$nodes,
##                     h = tp$h,
##                     nparents = tp$parents,
##                     conjunction = TRUE,
##                     multilevelParent = TRUE,
##                     removeIndirect = TRUE,
##                     rootName = "0"
##                     )

##     am11 <<- am1
    
##     ## FIXME: for now, we only excercise names = 0.
##     ## Posets only deal with integer codes, and same for restrictTable.
    
##     if(names == 1) {
##         rn <- replicate(nrow(am1) - 1,
##                         paste(sample(letters, 10, replace = TRUE),
##                               collapse = ""))
##         rn <- c("Root", rn)
##     } else if(names == 2) {
##         rn <- replicate(nrow(am1) - 1,
##                         paste(sample(letters, 10, replace = TRUE),
##                               collapse = ""))
##         rn <- c("0", rn)
##     } else if(names == 3) {
##         rn <- paste0("G", seq.int(nrow(am1) - 1))
##         rn <- c("Root", rn)
##     } else if(names == 4) {
##         rn <- paste0("G", seq.int(nrow(am1) - 1))
##         rn <- c("0", rn)
##     }
##     if(names != 0)
##         rownames(am1) <- colnames(am1) <- rn

    
    
##     am1.nr <- am1[-1, -1]
##     gf1 <- as(am1, "graphNEL")
## ##    p1 <- OncoSimulR:::adjMatToPoset(am1) ## am1.To.p1

## ##     p1.nr <- OncoSimulR:::adjMatToPoset(am1, dropRoot = TRUE)
    
##     ## p1.am1 <- OncoSimulR:::posetToGraph(p1,
##     ##                                     names = 0:tp$nodes,
##     ##                                     addroot = TRUE,
##     ##                                     type = "adjmat")

## ##    p1.am1.p1 <- OncoSimulR:::adjMatToPoset(p1.am1)

    
##     ## p1.am1.nr <- OncoSimulR:::posetToGraph(p1,
##     ##                                        names = 1:tp$nodes,
##     ##                                        addroot = FALSE,
##     ##                                        type = "adjmat")
    

##     ## p1.nr.am1 <- OncoSimulR:::posetToGraph(p1.nr,
##     ##                                        names = 0:tp$nodes,
##     ##                                        addroot = TRUE,
##     ##                                        type = "adjmat")

##     ## p1.nr.am1.nr <- OncoSimulR:::posetToGraph(p1.nr,
##     ##                                           names = 1:tp$nodes,
##     ##                                           addroot = FALSE,
##     ##                                           type = "adjmat")


##     am1.gf1.am1 <- as(as(am1, "graphNEL"), "matrix")
##     gf1.am1.gf1 <- as(as(gf1, "matrix"), "graphNEL")

##     gf1.nr <-  as(am1[-1, -1], "graphNEL")
    
##     NR.am1.gf1.am1 <- as(as(am1.nr, "graphNEL"), "matrix")
##     NR.gf1.am1.gf1 <- as(as(gf1.nr, "matrix"), "graphNEL")

    
## ##    gf1.p1 <- OncoSimulR:::graphToPoset(gf1)
##     ## p1.gf1 <- OncoSimulR:::posetToGraph(p1,
##     ##                        names = 0:tp$nodes,
##     ##                        addroot = TRUE,
##     ##                        type = "graphNEL")
##     ## p1.gf1.p1 <- OncoSimulR:::graphToPoset(p1.gf1)

##     ## p1.gf1.nr <- OncoSimulR:::posetToGraph(p1,
##     ##                           names = 1:tp$nodes,
##     ##                           addroot = FALSE,
##     ##                           type = "graphNEL")

##     ## these cannot be done now, as adj matrices without 0 are disabled
##     ## gf1.nr.p1 <- OncoSimulR:::graphToPoset(gf1.nr)
##     ## gf1.nr.p1.B <- OncoSimulR:::graphToPoset(p1.gf1.nr)

##     ## p1.nr.gf1 <- OncoSimulR:::posetToGraph(p1.nr,
##     ##                           names = 0:tp$nodes,
##     ##                           addroot = TRUE,
##     ##                           type = "graphNEL")
##     ## p1.nr.gf1.nr <- OncoSimulR:::posetToGraph(p1.nr,
##     ##                           names = 1:tp$nodes,
##     ##                           addroot = FALSE,
##     ##                           type = "graphNEL")


##     am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1, root = TRUE)
##     am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)
##     ## p1.To.rt <- OncoSimulR:::poset.to.restrictTable(p1)
##     ## p1.To.rt.2 <- OncoSimulR:::poset.to.restrictTable(p1.nr)
   
   
##     ## p1.gf1.p1 <- graphToPoset(posetTograph(p1, keeproot = TRUE))
##     ## gf1.p1.gf1 <- posetToGraph(graphToPoset(gf1), keeproot = TRUE)

    

    
    
##     ## am1.To.p1.To.am1 <- OncoSimulR:::posetToGraph(am1.To.p1,
##     ##                                               keeproot = TRUE,
##     ##                                               type = "adjmat")

##     ## posets must have root
##     ## am1.To.p1.To.am1.No.Root <- OncoSimulR:::posetToGraph(am1.To.p1,
##     ##                                                       keeproot = FALSE,
##     ##                                                       type = "adjmat")

    
##     ## if(names == 0) {
##     ##     am1.To.rt <- OncoSimulR:::adjmat.to.restrictTable(am1)
##     ##     ## am1.To.rt.2 <- OncoSimulR:::adjmat.to.restrictTable(am1[-1, -1], root = FALSE)
##     ##     p1.To.rt <- OncoSimulR:::poset.to.restrictTable(p1)
##     ## } else {
##     ##     am1.To.rt <- NA
##     ##     ## am1.To.rt.2 <- NA
##     ##     p1.To.rt <- NA
##     ## }

##     return(list(
##         am1 = am1,
##         am1.nr = am1.nr,
##         gf1 = gf1,
## ##        p1 = p1,
## ##        p1.nr = p1.nr,
## ##        p1.am1 = p1.am1,
## ##        p1.am1.p1 = p1.am1.p1,
## ##        p1.am1.nr = p1.am1.nr,
## ##        p1.nr.am1 = p1.nr.am1,
## ##        p1.nr.am1.nr = p1.nr.am1.nr,
##         am1.gf1.am1 = am1.gf1.am1,
##         gf1.am1.gf1 = gf1.am1.gf1,
##         gf1.nr = gf1.nr,
##         NR.am1.gf1.am1 = NR.am1.gf1.am1,
##         NR.gf1.am1.gf1 = NR.gf1.am1.gf1,
##         ## gf1.p1 = gf1.p1,
## ##        p1.gf1 = p1.gf1,
##         ## p1.gf1.p1 = p1.gf1.p1,
## ##        p1.gf1.nr = p1.gf1.nr,
##         ## gf1.nr.p1 = gf1.nr.p1,
##         ## gf1.nr.p1.B = gf1.nr.p1.B,
## ##        p1.nr.gf1 = p1.nr.gf1,
## ##        p1.nr.gf1.nr = p1.nr.gf1.nr,
##         am1.To.rt   = am1.To.rt, 
##         am1.To.rt.2 = am1.To.rt.2 
##         ## p1.To.rt    = p1.To.rt, 
##         ## p1.To.rt.2  = p1.To.rt.2        
##     ))  
##     ## return(list(am1 = am1,
##     ##             gf1 = gf1,
##     ##             p1 = p1,
##     ##             p1.am1 = p1.am1,
##     ##             p1.am1.p1 = p1.am1.p1,
##     ##             am1.gf1.am1 = am1.gf1.am1,
##     ##             gf1.am1.gf1 = gf1.am1.gf1,
##     ##             p1.gf1.p1 = p1.gf1.p1,
##     ##             gf1.p1.gf1 = gf1.p1.gf1,
##     ##             am1.nr = am1.nr,
##     ##             gf1.nr = gf1.nr,
##     ##             NR.am1.gf1.am1 = NR.am1.gf1.am1,
##     ##             NR.gf1.am1.gf1 = NR.gf1.am1.gf1))
## }





## ## FIXME
## ## For now, any conversion to rT is from things with integer names
## ## between other things with arbitrary names??




## ## posetToGraph as graphnel output

## ## am -> graphNel  and graphNel -> am
## ## graphNel -> poset and poset -> graphNEL


## ##  am <-> poset
## ##  am  <-> graphNel
## ##  graphNel <-> poset

## ## Nice full circle: am -> graphNel -> poset -> am
## ## am -> poset -> graphNEL
## ## etc, starting at each point


## masterTestCall <- function(rangeNodes = 4:30,
##                            rangeParents = 2:5,
##                            verbose = TRUE,
##                            names = 0) {

##     out <- createAndConvert(rangeNodes = rangeNodes,
##                            rangeParents = rangeParents,
##                             verbose = verbose,
##                             names = names)


##     ## test_that("adjmat->poset->adjmat", {
##     ##     expect_identical(out$am1, out$p1.am1)
##     ## })
##     ## ## the next one is currently redundant
##     ## test_that("poset->adjmat->poset", {
##     ##     expect_identical(out$p1, out$p1.am1.p1)
##     ## })

##     ## test_that("adjmat->poset->adjmat, no root", {
##     ##     expect_identical(out$am1.nr, out$p1.am1.nr)
##     ## })


##     ## test_that("sizes return adj mats w/w.o nr", {
##     ##     expect_true( nrow(out$p1.am1) == (nrow(out$p1.am1.nr) + 1) ) 
##     ## })

##     ## test_that("different posets, p1, p1.nr", {
##     ##     expect_false(identical(out$p1 , out$p1.nr))
##     ## })

##     ## test_that("p1 larger than p1.nr", {
##     ##     expect_true( nrow(out$p1) > nrow(out$p1.nr) )
##     ## })
   
    
##     ## test_that("adjmat->poset.nr ->adjmat", {
##     ##     expect_identical(out$am1, out$p1.nr.am1)
##     ## })

##     ## test_that("adjmat->poset.nr ->adjmat, no root", {
##     ##     expect_identical(out$am1.nr, out$p1.nr.am1.nr)
##     ## })

##     ## next two are redundant
##     test_that("adjmat->graph->adjmat", {
##         ## equality, not identical, since storage.mode differs
##         expect_equal(out$am1, out$am1.gf1.am1)
##     })
##     test_that("graph->adjmat->graph", {
##         expect_equal(out$gf1, out$gf1.am1.gf1)
##     })
    
##     ## next two are redundant
##     test_that("adjmat->graph->adjmat, no root", {
##         expect_equal(out$am1.nr, out$NR.am1.gf1.am1)
##     })
##     test_that("graph->adjmat->graph, no root", {
##         expect_equal(out$gf1.nr, out$NR.gf1.am1.gf1)
##     })


##     ## test_that("poset = (graph->poset)", {
##     ##     expect_identical(out$p1, out$gf1.p1)
##     ## })
    
##     ## test_that("graph =   (poset -> graph) ", {
##     ##     expect_identical(out$gf1, out$p1.gf1)
##     ## })

##     ## ## redundant
##     ## test_that("poset -> graph->poset ", {
##     ##     expect_identical(out$p1, out$p1.gf1.p1)
##     ## })

##     ## test_that("poset -> graph.nr  = gf1.nr ", {
##     ##     expect_identical(out$gf1.nr, out$p1.gf1.nr)
##     ## })

##     ## No adjanceny matrices without 0 allowed
##     ## test_that("graph.nr -> poset  =  am.nr -> poset ", {
##     ##     expect_identical(out$gf1.nr.p1, out$p1.nr)
##     ## })

##     ## test_that("graph.nr -> poset =  poset-> graf.nr -> poset ", {
##     ##     expect_identical(out$gf1.nr.p1, out$gf1.nr.p1.B)
##     ## })

##     ## test_that("gf1 =  poset.nr ->graph", {
##     ##     expect_identical(out$gf1, out$p1.nr.gf1)
##     ## })

##     ## test_that("gf1.nr =  poset.nr ->graph nr", {
##     ##     expect_identical(out$gf1.nr, out$p1.nr.gf1.nr)
##     ## })

##     ## test_that(" from p1.nr to gf1 w and w/o root differ", {
##     ##     expect_false(identical(out$p1.nr.gf1, out$p1.nr.gf1.nr))
##     ## })

##     ## test_that("graph->poset->graph", {
##     ##     expect_identical(out$gf1, out$gf1.p1.gf1)
##     ## })
  
    
##     ## if(names == 0) {
##         ## test_that("restriction table: identical from adjmat and poset", {
##         ##     expect_identical(out$am1.To.rt, out$p1.To.rt)
##         ## })
##         test_that("restriction table: identical from adjmat w/wo root", {
##             expect_identical(out$am1.To.rt, out$am1.To.rt.2)
##         })
##         ## test_that("restriction table: identical from poset w/wo root", {
##         ##     expect_identical(out$p1.To.rt, out$p1.To.rt.2)
##         ## })
##   ##   }
##     ## cat("A full round of tests completed OK.\n")
##     return("OK")
## } 







## tmp <- replicate(numTests, masterTestCall() )


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

cat(paste("\n Ending poset-transformations tests", date(), "\n"))
