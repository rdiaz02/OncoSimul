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
    lop8 <- suppressWarnings(LOD(pancr8))
    OncoSimulR:::LOD_as_path(lop8)
    expect_true(inherits(POM(pancr1), "character"))
    require(igraph)
    expect_true(inherits(LOD(pancr1, strict = FALSE)$all_paths[[1]], "igraph.vs"))
    expect_true(is.na(LOD(pancr1, strict = TRUE)$all_paths))
    expect_true(inherits(LOD(pancr1)$lod_single, "igraph.vs"))
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
    pancr8 <- suppressWarnings(suppressMessages(oncoSimulPop(30, pancr, model = "McFL",
                                            keepPhylog = TRUE,
                                            onlyCancer = FALSE,
                                            finalTime = 0.01,
                                            sampleEvery = 0.5,
                                            mu = 1e-8,
                                            mc.cores = 2,
                                            mutationPropGrowth = FALSE,
                                            initSize = 10)))
    lop8 <- suppressWarnings(LOD(pancr8))
    lop8b <- suppressWarnings(LOD(pancr8, strict = TRUE))
    lop8c <- suppressWarnings(LOD(pancr8, strict = FALSE))
    OncoSimulR:::LOD_as_path(lop8[[1]])
    OncoSimulR:::LOD_as_path(lop8)
    ## what if this is commented out?
    ## expect_true(any(unlist(lapply(lop8,
    ##                               function(x) x$lod_single))
    ##                 %in% "No_descendants"))
    ## there are descendants but all go extinct
    ## pancr2 <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
    ##                                       "TP53", "TP53", "MLL3"),
    ##                                   child = c("KRAS","SMAD4", "CDNK2A", 
    ##                                       "TP53", "MLL3",
    ##                                       rep("PXDN", 3), rep("TGFBR2", 2)),
    ##                                   s = -0.9,
    ##                                   sh = -0.9,
    ##                                   typeDep = "MN"))
    gg <- allFitnessEffects(noIntGenes = rep(-.9, 100))
    pancr22 <- oncoSimulPop(10, gg,
                            model = "Exp",
                            keepPhylog = TRUE,
                            onlyCancer = FALSE,
                            initSize = 1e3,
                            mu = 1e-2,
                            mc.cores = 2,
                            finalTime = 2.5)
    lp22 <- LOD(pancr22)
    lp23 <- LOD(pancr22, strict = TRUE)
    lp24 <- LOD(pancr22, strict = FALSE)
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
    expect_warning(LOD(s1, strict = FALSE),
                   "There never was a descendant of WT",
                   fixed = TRUE)
    expect_warning(LOD(s1, strict = FALSE),
                   "PhylogDF has 0 rows:",
                   fixed = TRUE)
    ## expect_warning(LOD(s1, strict = TRUE),
    ##                "There never was a descendant of WT",
    ##                fixed = TRUE)
    expect_warning(LOD(s1, strict = TRUE),
                   "LOD structure has 0 rows:",
                   fixed = TRUE)
    s2 <- oncoSimulIndiv(allFitnessEffects(genotFitness = m1),
                         mu = 1e-14, detectionSize = 1, initSize = 100,
                         keepPhylog = FALSE)
    expect_warning(LOD(s2, strict = TRUE),
                   "LOD structure has 0 rows:",
                   fixed = TRUE)
})

date()
test_that("LOD, strict, same as would be obtained from full structure", {

    n <- 10
    for(i in 1:n) {
        ng <- 6
        rxx <- rfitness(ng)
        rxx[sample(2:(ng + 1)), ng + 1] <- 1.5 ## make sure we get going
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 1000, detectionSize = 1e6,
                             keepPhylog = TRUE, mu = 1e-3)
        lods <- LOD(s7, strict = TRUE)
        loda <- LOD(s7, strict = FALSE)
        ## lods should be among the loda
        if(!is.null(s7$pops.by.time)) {
            expect_true(any(
                unlist(lapply(loda$all_paths,
                              function(x) identical(names(x),
                                                    names(lods$lod_single))))))
            if(length(loda$all_paths) == 1) {
                expect_true(identical(names(loda$lod_single),
                                      names(lods$lod_single)))
            }
            ## print(OncoSimulR:::LOD_as_path(lods))
            ## print(OncoSimulR:::LOD_as_path(loda))
        }
    }
    ## with initMutant
    for(i in 1:n) {
        rxx <- rfitness(6)
        rxx[3, 7] <- 1.5        
        s7 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx),
                             initSize = 1000, detectionSize = 1e6,
                             keepPhylog = TRUE, mu = 1e-3,
                             initMutant = c("B"))
        lods <- LOD(s7, strict = TRUE)
        loda <- LOD(s7, strict = FALSE)
        ## lods should be among the loda
        expect_true(any(
            unlist(lapply(loda$all_paths,
                          function(x) identical(names(x),
                                                names(lods$lod_single))))))
        if(length(loda$all_paths) == 1) {
            expect_true(identical(names(loda$lod_single),
                                  names(lods$lod_single)))
        }
        ## print(loda)
    }
})
date()






## To see how the above works by returning LOD sensu stricto you can look
## at this code:



pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                               "TP53", "TP53", "MLL3"),
                                           child = c("KRAS","SMAD4", "CDNK2A", 
                                               "TP53", "MLL3",
                                               rep("PXDN", 3), rep("TGFBR2", 2)),
                                           s = 0.05,
                                           sh = -0.3,
                                           typeDep = "MN"))
     
     
pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)

## All we need to get LOD sensu stricto (i.e., identical to Szendro)
## is keep pop size of receiving, or destination, genotype.
## Then, filter those where popSize > 0


## And use first (starting from bottom) of that path
## So find, for each child, the last event with popSize == 0,
## and keep only that row.

## Probably enough to run duplicated in reverse (on the df with popSize
## child = 0)

## The indices to keep: !rev(duplicated(rev(fg3[, 2])))

fg1 <- pancr1$other$PhylogDF


fg3 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
                  child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
                  time = 1:6,
                  pop_size_child = c(0, 0, 0, 0, 0, 0),
                  stringsAsFactors = FALSE)


fg4 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
                  child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
                  time = 1:6,
                  pop_size_child = c(0, 0, 0, 2, 0, 0),
                  stringsAsFactors = FALSE)

## ## from phylogClone, key parts
## fpc <- function(df) {
##     tG <- unique(c(df[, 1], df[, 2]))
##     g <- igraph::graph.data.frame(df[, c(1, 2)])
##     nodesInP <- unique(unlist(igraph::neighborhood(g, order = 1e+09, 
##                                                    nodes = tG, mode = "in")))
##     allLabels <- unique(as.character(unlist(df[, c(1, 2)])))
##     nodesRm <- setdiff(allLabels, V(g)$name[nodesInP])
##     g <- igraph::delete.vertices(g, nodesRm)
##     tmp <- list(graph = g, df = df)
##     class(tmp) <- c(class(tmp), "phylogClone")
##     return(tmp)
## }

## ## Filter the PhylogDF so we obtain LOD, sensu stricto.
## filter_phylog_df_LOD_ <- function(x) {
##     x <- x[x$pop_size_child == 0, ]
##     keep <- !rev(duplicated(rev(x$child)))
##     return(x[keep, ])
## }

require(igraph)
all_simple_paths(OncoSimulR:::phcl_from_lod(OncoSimulR:::filter_phylog_df_LOD_with_n(fg3))$graph,
                 from = "",
                 to = "A, B, C",
                 mode = "out")


all_simple_paths(OncoSimulR:::phcl_from_lod(OncoSimulR:::filter_phylog_df_LOD_with_n(fg3))$graph,
                 to = "",
                 from = "A, B, C",
                 mode = "in")

