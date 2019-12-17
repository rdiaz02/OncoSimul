inittime <- Sys.time()
cat(paste("\n Starting LOD_POM at", date(), "\n"))
date()
test_that("Exercise LOD and POM code", {
    pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
    pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)
    pancr8 <- oncoSimulPop(6, pancr, model = "Exp", keepPhylog = TRUE,
                           detectionSize = 1e5,
                           mc.cores = 2)
    lop8 <- LOD(pancr8)
    OncoSimulR:::LOD_as_path(lop8)
    expect_true(inherits(POM(pancr1), "character"))
    require(igraph)
    ## expect_true(inherits(LOD(pancr1, strict = FALSE)$all_paths[[1]], "igraph.vs"))
    ## expect_true(inherits(LOD(pancr1), "igraph.vs"))
    expect_true(inherits(LOD(pancr1), "character"))
    expect_true(inherits(POM(pancr8), "list"))
    expect_true(inherits(LOD(pancr8), "list"))
    expect_true(inherits(diversityPOM(POM(pancr8)), "numeric"))
    expect_true(inherits(diversityLOD(LOD(pancr8)), "numeric"))
    expect_true(diversityPOM(POM(pancr8)) >= 0)
    expect_true(diversityLOD(LOD(pancr8)) >= 0)
    expect_error(diversityPOM(POM(pancr1)),
                 "Object must be a list", fixed = TRUE)
    expect_error(diversityLOD(LOD(pancr1)),
                 "Object must be a list", fixed = TRUE)
    pancr88 <- oncoSimulPop(8, pancr, model = "McFL",
                           keepPhylog = TRUE,
                           finalTime = 0.01,
                           max.num.tries = 1,
                           mc.cores = 2,
                           max.wall.time = 0.01,
                           detectionSize = 1e6)
    expect_warning(LOD(pancr88),
                   "Missing needed components.", fixed = TRUE)
    expect_warning(POM(pancr88),
                   "Missing needed components.", fixed = TRUE)
    pancr8 <- suppressWarnings(suppressMessages(oncoSimulPop(20,
                                                             pancr, model = "McFL",
                                            keepPhylog = TRUE,
                                            onlyCancer = FALSE,
                                            max.num.tries = 2,
                                            finalTime = 0.01,
                                            sampleEvery = 0.5,
                                            mu = 1e-8,
                                            mc.cores = 2,
                                            mutationPropGrowth = FALSE,
                                            initSize = 10)))
    lop8 <- suppressWarnings(LOD(pancr8))
    lop8b <- suppressWarnings(LOD(pancr8))
    OncoSimulR:::LOD_as_path(lop8[[1]])
    OncoSimulR:::LOD_as_path(lop8)
    gg <- allFitnessEffects(noIntGenes = rep(-.9, 100))
    pancr22 <- oncoSimulPop(6, gg,
                            model = "Exp",
                            keepPhylog = TRUE,
                            onlyCancer = FALSE,
                            max.num.tries = 2,
                            initSize = 1e3,
                            mu = 1e-2,
                            mc.cores = 2,
                            finalTime = 2.5)
    lp22 <- LOD(pancr22)
    ## There is like soooo remote chance this will fail
    ## and the previous exercises the code anyway.
    ## expect_true(any(unlist(lp22) %in% "WT_is_end"))
})
date()

test_that("Warnings when no descendants",  {
    ## cannot move from wt with this fitness landscape
    m1 <- cbind(A = c(0, 1), B = c(0, 1), Fitness = c(1, 1e-8))
    s1 <- oncoSimulIndiv(allFitnessEffects(genotFitness = m1),
                         mu = 1e-14, detectionSize = 1, initSize = 100,
                         keepPhylog = TRUE)
    expect_warning(LOD(s1),
                   "LOD structure has 0 rows:",
                   fixed = TRUE)
    s2 <- oncoSimulIndiv(allFitnessEffects(genotFitness = m1),
                         mu = 1e-14, detectionSize = 1, initSize = 100,
                         keepPhylog = FALSE)
    expect_warning(LOD(s2),
                   "LOD structure has 0 rows:",
                   fixed = TRUE)
})

## Done better below and make sure it runs with keepPhylog = FALSE
## date()
## test_that("LOD, strict, same as would be obtained from full structure", {
##     ## we are testing in an extremely paranoid way, against a
##     ## former version
##     n <- 10
##     for(i in 1:n) {
##         ng <- 6
##         rxx <- rfitness(ng)
##         rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
##         s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
##                              initSize = 1000, detectionSize = 1e6,
##                              keepPhylog = TRUE, mu = 1e-3)
##         lods <- LOD(s7)
##         loda <- OncoSimulR:::LOD.oncosimul2_pre_2.9.2(s7, strict = FALSE)
##         ## lods should be among the loda
##         if(!is.null(s7$pops.by.time)) {
##             expect_true(any(
##                 unlist(lapply(loda$all_paths,
##                               function(x) identical(names(x),
##                                                     names(lods))))))
##             if(length(loda$all_paths) == 1) {
##                 expect_true(identical(names(loda$lod_single),
##                                       names(lods)))
##             }
##             ## print(OncoSimulR:::LOD_as_path(lods))
##             ## print(OncoSimulR:::LOD_as_path(loda))
##         }
##     }
## })
## date()

## Done better below and make sure it runs with keepPhylog = FALSE
## date()
## test_that("LOD, strict, same as would be obtained from full structure, with initMutant", {
##     n <- 10
##     ## with initMutant
##     for(i in 1:n) {
##         rxx <- rfitness(6)
##         rxx[3, 7] <- 1.5        
##         s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
##                              initSize = 1000, detectionSize = 1e6,
##                              keepPhylog = TRUE, mu = 1e-3,
##                              initMutant = c("B"))
##         lods <- LOD(s7)
##         loda <- OncoSimulR:::LOD.oncosimul2_pre_2.9.2(s7, strict = FALSE)
##         ## lods should be among the loda
##         expect_true(any(
##             unlist(lapply(loda$all_paths,
##                           function(x) identical(names(x),
##                                                 names(lods))))))
##         if(length(loda$all_paths) == 1) {
##             expect_true(identical(names(loda$lod_single),
##                                   names(lods)))
##         }
##         ## print(loda)
##     }
## })
## date()

set.seed(NULL)
si <- runif(1, 1, 1e9)
print(si)
date()
test_that("LOD, strict, same as would be obtained from full structure, seed", {
    ## we are testing in an extremely paranoid way, against a
    ## former version
    n <- 10
    for(i in 1:n) {
        ng <- 8
        rxx <- rfitness(ng, c = 1, sd = 0.1,
                        reference = rep(1, ng))
        ## rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
        set.seed(si + i)
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 10, detectionSize = 1e5,
                             keepPhylog = FALSE, mu = 5e-3)
        set.seed(si + i)
        s7b <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 10, detectionSize = 1e5,
                             keepPhylog = TRUE, mu = 5e-3)
        lods <- LOD(s7)
        print(lods)
        loda <- OncoSimulR:::LOD.oncosimul2_pre_2.9.2(s7b, strict = FALSE)
        ## lods should be among the loda
        if(!is.null(s7$pops.by.time)) {
            expect_true(any(
                unlist(lapply(loda$all_paths,
                              function(x) identical(names(x),
                                                    lods)))))
            if(length(loda$all_paths) == 1) {
                expect_true(identical(names(loda$lod_single),
                                      lods))
            }
            ## print(OncoSimulR:::LOD_as_path(lods))
            ## print(OncoSimulR:::LOD_as_path(loda))
        }
    }
})
date()

set.seed(NULL)
si <- runif(1, 1, 1e9)
print(si)
date()
test_that("LOD, strict, same as would be obtained from full structure, with initMutant", {
    n <- 10
    ## with initMutant
    for(i in 1:n) {
        rxx <- rfitness(8, reference = rep(1, 8))
        rxx[4, ncol(rxx)] <- 1.5
        set.seed(2 * si + i)
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 100, detectionSize = 1e5,
                             keepPhylog = FALSE, mu = 5e-3,
                             initMutant = c("C"))
        set.seed(2 * si + i)
        s7b <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                              initSize = 100, detectionSize = 1e5,
                              keepPhylog = TRUE, mu = 5e-3,
                              initMutant = c("C"))
        lods <- LOD(s7)
        print(lods)
        attributes(lods) <- NULL
        loda <- OncoSimulR:::LOD.oncosimul2_pre_2.9.2(s7b, strict = FALSE)
        ## we need this because o.w. the old output it ain't an igraph
        ## object with names
        if(!any(grepl("_is_end", lods)) && !any(grepl("No_descendants", lods))) {
            ## lods should be among the loda
            if(!is.null(s7$pops.by.time)) {
                expect_true(any(
                    unlist(lapply(loda$all_paths,
                                  function(x) identical(names(x),
                                                        lods)))))
                if(length(loda$all_paths) == 1) {
                    expect_true(identical(names(loda$lod_single),
                                          lods))
                }
            }
        } else if (any(grepl("No_descendants", lods))) {
            expect_true(identical(loda$lod_single,
                                  lods)) 
        } else {
            if(!is.null(s7$pops.by.time)) {
                expect_true(any(
                    unlist(lapply(loda$all_paths,
                                  function(x) identical(x,
                                                        lods)))))
                if(length(loda$all_paths) == 1) {
                    expect_true(identical(loda$lod_single,
                                          lods))
                }
            }  
        }
    }
    
})
date()
set.seed(NULL)


date()
test_that("POM from C++ is the same as from the pops.by.time object", {
    ## Must make sure keepEvery = sampleEvery or granularity of
    ## C++ can be larger
    n <- 10
    for(i in 1:n) {
        ng <- 6
        rxx <- rfitness(ng)
        rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 1000, detectionSize = 1e6,
                             mu = 1e-3)
        pom <- OncoSimulR:::POM_pre_2.9.2(s7)
        if(!is.null(s7$pops.by.time) &&
           !any(apply(s7$pops.by.time[, -1], 1, function(x) length(which(x == max(x))) > 1)))
            expect_true(identical(s7$other$POM, pom))
    }
    ## with initMutant
    for(i in 1:n) {
        rxx <- rfitness(6)
        rxx[3, 7] <- 1.5        
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 1000, detectionSize = 1e6,
                             mu = 1e-3,
                             initMutant = c("B"))
        pom <- OncoSimulR:::POM_pre_2.9.2(s7)
        ## if(!is.null(s7$pops.by.time)) {
        if(!is.null(s7$pops.by.time) &&
           !any(apply(s7$pops.by.time[, -1, drop = FALSE], 1,
                      function(x) length(which(x == max(x))) > 1)))
            expect_true(identical(s7$other$POM, pom))
    }
    ## try to make extinction likely
    for(i in 1:n) {
        rxx <- rfitness(6)
        rxx[3, 7] <- 1e-8
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 10, detectionSize = 1e6,
                             mu = 1e-3,
                             max.num.tries = 3,
                             errorHitMaxTries = FALSE,
                             initMutant = c("B"))
        pom <- OncoSimulR:::POM_pre_2.9.2(s7)
        if(!is.null(s7$pops.by.time) &&
           !any(apply(s7$pops.by.time[, -1, drop = FALSE], 1,
                      function(x) length(which(x == max(x))) > 1))) {
            if(any(s7$other$POM == "_EXTINCTION_"))
                expect_true(identical(s7$other$POM[-length(s7$other$POM)], pom))
            else
                expect_true(identical(s7$other$POM, pom))
        }
    }
})
date()




## can be unpredictably slow. Not needed.
## date()
## test_that("POM from C++ is the same as from the pops.by.time object, McFL", {
##     ## Must make sure keepEvery = sampleEvery or granularity of
##     ## C++ can be larger
##     n <- 10
##     for(i in 1:n) {
##         ng <- 6
##         rxx <- rfitness(ng)
##         rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
##         s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
##                              initSize = 1000, detectionSize = 1e4,
##                              mu = 1e-3, model = "McFL")
##         pom <- OncoSimulR:::POM_pre_2.9.2(s7)
##         if(!is.null(s7$pops.by.time) &&
##            !any(apply(s7$pops.by.time[, -1], 1, function(x) length(which(x == max(x))) > 1)))
##             expect_true(identical(s7$other$POM, pom))
##     }
##     ## with initMutant
##     for(i in 1:n) {
##         rxx <- rfitness(6)
##         rxx[3, 7] <- 1.5        
##         s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
##                              initSize = 1000, detectionSize = 1e4,
##                              mu = 1e-3, model = "McFL",
##                              initMutant = c("B"))
##         pom <- OncoSimulR:::POM_pre_2.9.2(s7)
##         ## if(!is.null(s7$pops.by.time)) {
##         if(!is.null(s7$pops.by.time) &&
##            !any(apply(s7$pops.by.time[, -1], 1, function(x) length(which(x == max(x))) > 1)))
##             expect_true(identical(s7$other$POM, pom))
##     }
## })
## date()





## ## To see how the above works by returning LOD sensu stricto you can look
## ## at this code:



## pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
##                                                "TP53", "TP53", "MLL3"),
##                                            child = c("KRAS","SMAD4", "CDNK2A", 
##                                                "TP53", "MLL3",
##                                                rep("PXDN", 3), rep("TGFBR2", 2)),
##                                            s = 0.05,
##                                            sh = -0.3,
##                                            typeDep = "MN"))
     
     
## pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)

## ## All we need to get LOD sensu stricto (i.e., identical to Szendro)
## ## is keep pop size of receiving, or destination, genotype.
## ## Then, filter those where popSize > 0


## ## And use first (starting from bottom) of that path
## ## So find, for each child, the last event with popSize == 0,
## ## and keep only that row.

## ## Probably enough to run duplicated in reverse (on the df with popSize
## ## child = 0)

## ## The indices to keep: !rev(duplicated(rev(fg3[, 2])))

## fg1 <- pancr1$other$PhylogDF


## fg3 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
##                   child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
##                   time = 1:6,
##                   pop_size_child = c(0, 0, 0, 0, 0, 0),
##                   stringsAsFactors = FALSE)


## fg4 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
##                   child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
##                   time = 1:6,
##                   pop_size_child = c(0, 0, 0, 2, 0, 0),
##                   stringsAsFactors = FALSE)

## ## ## from phylogClone, key parts
## ## fpc <- function(df) {
## ##     tG <- unique(c(df[, 1], df[, 2]))
## ##     g <- igraph::graph.data.frame(df[, c(1, 2)])
## ##     nodesInP <- unique(unlist(igraph::neighborhood(g, order = 1e+09, 
## ##                                                    nodes = tG, mode = "in")))
## ##     allLabels <- unique(as.character(unlist(df[, c(1, 2)])))
## ##     nodesRm <- setdiff(allLabels, V(g)$name[nodesInP])
## ##     g <- igraph::delete.vertices(g, nodesRm)
## ##     tmp <- list(graph = g, df = df)
## ##     class(tmp) <- c(class(tmp), "phylogClone")
## ##     return(tmp)
## ## }

## ## ## Filter the PhylogDF so we obtain LOD, sensu stricto.
## ## filter_phylog_df_LOD_ <- function(x) {
## ##     x <- x[x$pop_size_child == 0, ]
## ##     keep <- !rev(duplicated(rev(x$child)))
## ##     return(x[keep, ])
## ## }

## require(igraph)
## all_simple_paths(OncoSimulR:::phcl_from_lod(OncoSimulR:::filter_phylog_df_LOD_with_n(fg3))$graph,
##                  from = "",
##                  to = "A, B, C",
##                  mode = "out")


## all_simple_paths(OncoSimulR:::phcl_from_lod(OncoSimulR:::filter_phylog_df_LOD_with_n(fg3))$graph,
##                  to = "",
##                  from = "A, B, C",
##                  mode = "in")

cat(paste("\n Ending LOD_POM at", date(), "\n"))

### Why we need to exclude some POMs in the testing

## with i = 2335   ## new one removes entries
##    because two have identical values
##    apply(s7$pops.by.time[, -1], 1, function(x) which(x == max(x)))
## with i = 10001421 ## different entries
##    same problem: a case of two pops with identical values
## with i = 15046  ## new one adds entries
##    same problem
##    OK if we account for that

## while(TRUE) {
##     i <- i + 1
##     set.seed(i)
##     ng <- 6
##     rxx <- rfitness(ng)
##     rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
##     s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
##                          initSize = 1000, detectionSize = 1e6,
##                          mu = 1e-3)
##     pom <- OncoSimulR:::POM_pre_2.9.2(s7)
##     if(!is.null(s7$pops.by.time))
##         expect_true(identical(s7$other$POM, pom))
## }

cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
