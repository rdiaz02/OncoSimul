inittime <- Sys.time()
cat(paste("\n Starting init-mutant tests", date(), "\n"))

RNGkind("Mersenne-Twister")
set.seed(NULL)

## Processing this file takes about 3 seconds on my laptop
test_that("initMutant crashes",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c(0.01, 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b, a") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b", "a") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("1", "2") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
          })


test_that("initMutant lexicog order",
{
    max.tries <- 4
    for(tries in 1:max.tries) {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c(0.01, 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      sampleEvery = 0.03,
                      keepEvery = 1,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_true( "d > m _ " %in% cn )
              expect_false( "m > d _" %in% cn )
              expect_false( "m > d > f _" %in% cn )
              ## this next one will be true almost always
              ## if there is pop growth. But very occasionally
              ## it might not. The above must ALWAYS be true,
              ## but this one we allow to repeat a couple of times
              T1 <- ( "d > m > f _ " %in% cn )
              if( T1 ) break;
              if( !T1 ) {
                  cat("\n pop in initMutant lexicog order\n ")
                  print(tmp)
        }
    }
})


test_that("initMutant lexicog order with noint",
{
    max.tries <- 4
    for(tries in 1:max.tries) {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      sampleEvery = 0.03,
                      keepEvery = 1,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m > z")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_true( "d > m _ z" %in% cn )
              expect_false( "m > d _" %in% cn )
              expect_false( "m > d > f _" %in% cn )
              
              T1 <- ( "d > m > f _ z" %in% cn )
              ## expect_true( "d, m_z" %in% cn )
              ## expect_true( "d, m, f_z" %in% cn )
              ## expect_false( "m, d_" %in% cn )
              ## expect_false( "m, d, f_" %in% cn )
              if( T1 ) break;
              if(! T1 ) {
                  cat("\n pop in initMutant lexicog order with no int\n ")
                  print(tmp)
              }
    }
    ## cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})


test_that("initMutant non lexicog order", {
    max.tries <- 4
    for(tries in 1:max.tries) {
        o3 <- allFitnessEffects(orderEffects = c(
                                    "M > D > F" = 0.99,
                                    "D > M > F" = 0.2,
                                    "D > M"     = 0.1,
                                    "M > D"     = 0.9),
                                noIntGenes = c("u" = 0.01, "z" = 0.01),
                                geneToModule =
                                    c("Root" = "Root",
                                      "M" = "m",
                                      "F" = "f",
                                      "D" = "d") )
        tmp <- oncoSimulIndiv(o3, model = "McFL",
                              mu = 5e-5, finalTime = 500,
                              detectionDrivers = 3,
                              sampleEvery = 0.03,
                              keepEvery = 1,
                              onlyCancer = FALSE,
                              initSize = 1000,
                              keepPhylog = TRUE
                            , initMutant = c("m > d")
                              )
        cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                             returnGraph = TRUE),
                                             sparse = FALSE))
        expect_false( "d > m _" %in% cn )
        expect_false( "d > m > f _" %in% cn )
        expect_true( "m > d _ " %in% cn )
        ## this next one will be true almost always
        ## if there is pop growth. But very occasionally
        ## it might not. The above must ALWAYS be true,
        ## but this one we allow to repeat a couple of times
        T1 <- ( "m > d > f _ " %in% cn )
        if( T1 ) break;
        if( !T1 ) {
            cat("\n pop in initMutant non lexicog order\n ")
            print(tmp)
        }
    }
    ## cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})
  

test_that("initMutant non lexicog order",
{
    max.tries <- 4
    for(tries in 1:max.tries) {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      sampleEvery = 0.03,
                      keepEvery = 1,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("m > u > d")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_false( "d > m _" %in% cn )
              expect_false( "d > m > f _" %in% cn )
              expect_true( "m > d _ u" %in% cn )
              ## this next one will be true almost always
              ## if there is pop growth. But very occasionally
              ## it might not. The above must ALWAYS be true,
              ## but this one we allow to repeat a couple of times
              
              T1 <- ( "m > d > f _ u" %in% cn )
              if( T1 ) break;
              if(! T1 ) {
                  cat("\n pop in initMutant non lexicog order\n ")
                  print(tmp)
              }
    }
    ## cat(paste("\n done tries", tries, "\n"))
    expect_true(T1)
})

## FIXME: we could use stronger test: we will never see M > D
test_that("initMutant with oncoSimulSample", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = 0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d"),
                        drvNames = c("m", "f", "d"))
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 2,
                        onlyCancer = TRUE,
                        initSize = 500,
                        sampleEvery = 0.03,
                        initMutant = c("z > d"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("d", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("d", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("d",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant with oncoSimulSample, 2", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9,
                            "M > A"     = 0.25,
                            "A > H"     = 0.2,
                            "A > G"     = 0.3),
                        noIntGenes = c("u" = 0.1, 
                                       "v" = 0.2,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "A" = "a",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d",
                              "H" = "h",
                              "G" = "g"),
                        drvNames = c("m", "f", "d", "a", "h", "g"))
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        sampleEvery = 0.03,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > a"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("a", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("a", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("a",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant with oncoSimulPop", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.2,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d"),
                        drvNames = c("m", "f", "d"))
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        keepPhylog = TRUE,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        initSize = 500,
                        initMutant = c("d > m > y"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "y"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "d > m _ y"))
})



test_that("initMutant with oncoSimulPop, 2", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                            noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d"),
                        drvNames = c("m", "f", "d"))
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 70,
                        detectionDrivers = 4, ## yes, reach end
                        onlyCancer = FALSE,
                        keepPhylog = TRUE,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        initSize = 100,
                        initMutant = c("m > v > d"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "v"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "m > d _ v"))
##                   as.character(x$other$PhylogDF[1, 1])) == "m, d_v"))
})




test_that("initMutant with oncoSimulPop, McFL", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.2,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d"),
                        drvNames = c("m", "f", "d"))
    ospI <- oncoSimulPop(4, 
                        o3init, model = "McFL",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        keepPhylog = TRUE,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        detectionProb = NA,
                        initSize = 500,
                        initMutant = c("d > m > y"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "y"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "d > m _ y"))
})



test_that("initMutant with oncoSimulPop, Bozic", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.2,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d"),
                        drvNames = c("m", "f", "d"))
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Bozic",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        keepPhylog = TRUE,
                        sampleEvery = 0.03,
                        keepEvery = 1,
                        initSize = 500,
                        initMutant = c("d > m > y"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "y"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "d > m _ y"))
})



test_that("initMutant with oncoSimulSample, 2, McFL", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9,
                            "M > A"     = 0.25,
                            "A > H"     = 0.2,
                            "A > G"     = 0.3),
                        noIntGenes = c("u" = 0.1, 
                                       "v" = 0.2,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "A" = "a",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d",
                              "H" = "h",
                              "G" = "g"),
                        drvNames = c("m", "f", "d", "a", "h", "g")
                        )
    ossI <- oncoSimulSample(4, 
                        o3init, model = "McFL",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        sampleEvery = 0.03,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > a"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("a", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("a", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("a",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant with oncoSimulSample, 2, Bozic", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9,
                            "M > A"     = 0.25,
                            "A > H"     = 0.2,
                            "A > G"     = 0.3),
                        noIntGenes = c("u" = 0.1, 
                                       "v" = 0.2,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "A" = "a",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d",
                              "H" = "h",
                              "G" = "g"),
                        drvNames = c("m", "f", "d", "a", "h", "g"))
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Bozic",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        sampleEvery = 0.03,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > a"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("a", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("a", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("a",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant crashes if > number of genes", {
    o1 <- allFitnessEffects(
        noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))
    expect_error(oncoSimulIndiv(o1, initMutant = "b, a, c, d"),
                 "For driver or initMutant you have passed genes not in the fitness table",
                 fixed = TRUE)
})

test_that("initMutant crashes if not in fitness table even if fewer", {
     o1 <- allFitnessEffects(
         noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))
     expect_error(oncoSimulIndiv(o1, initMutant = "a, d"),
                  "For driver or initMutant you have passed genes not in the fitness table",
                  fixed = TRUE)
})


test_that("initMutant works if == number of genes", {
    o1 <- allFitnessEffects(
        noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))
    expect_silent(ooox <- oncoSimulIndiv(o1, initMutant = "b, a, c"))
                 ## "No mutable positions. Mutation set to dummyMutationRate",
                 ## fixed = TRUE)
    
    ## This used to crash
    set.seed(5)
    o2 <- allFitnessEffects(genotFitness = rfitness(2))
    oncoSimulIndiv(o2, initMutant = "B, A")

    ## Test a few
    for(i in 1:5){
        set.seed(i)
        o2i <- allFitnessEffects(genotFitness = rfitness(2))
        suppressWarnings(o2io <- oncoSimulIndiv(o2i, initMutant = "B, A",
                               onlyCancer = FALSE))
        if(!is.null(o2io$other$ExceptionMessage)) {
            expect_false(
                grepl("Trying to obtain a mutation when nonmutated.size is 0",
                      o2io$other$ExceptionMessage)
            )
        }
    }
    set.seed(NULL)
})


test_that("initMutant with freq-dep-fitness"  , {
    r <- data.frame(rfitness(2))
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

    r[, "Fitness"] <- c("f_ - f_1 - f_2 - f_1_2", 
                        "max(100*f_1, 10)", 
                        "max(100*f_2, 10)", 
                        "max((200*(f_1 + f_2) + 50*f_1_2), 1)")
    suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                             frequencyDependentFitness = TRUE))

    expect_s3_class(
        os1 <- oncoSimulIndiv(afe, 
                              model = "McFL",
                              initMutant = "A",
                              onlyCancer = FALSE, 
                              finalTime = 50, 
                              verbosity = 0, 
                              mu = 1e-6,
                              initSize = 500, 
                              keepPhylog = FALSE,
                              seed = NULL, 
                              errorHitMaxTries = TRUE, 
                              errorHitWallTime = TRUE),
        "oncosimul2")
    
    expect_true(is.null(os1$other$ExceptionMessage))
    
    expect_s3_class(os2 <- oncoSimulIndiv(afe, 
                          model = "McFL",
                          initMutant = "B",
                          onlyCancer = FALSE, 
                          finalTime = 50, 
                          verbosity = 0, 
                          mu = 1e-6,
                          initSize = 500, 
                          keepPhylog = FALSE,
                          seed = NULL, 
                          errorHitMaxTries = TRUE, 
                          errorHitWallTime = TRUE),
                    "oncosimul2")

    expect_s3_class(os3 <- oncoSimulIndiv(afe, 
                          model = "McFL",
                          initMutant = "A, B",
                          onlyCancer = FALSE, 
                          finalTime = 50, 
                          verbosity = 0, 
                          mu = 1e-6,
                          initSize = 500, 
                          keepPhylog = FALSE,
                          seed = NULL, 
                          errorHitMaxTries = TRUE, 
                          errorHitWallTime = TRUE),
                    "oncosimul2")
    
    expect_true(is.null(os3$other$ExceptionMessage))
    
})


test_that("WT initMutant simulation equiv. to no init mutant", {

    set.seed(1)
    r <- data.frame(rfitness(2))
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

    r[, "Fitness"] <- c("f_ - f_1 - f_2 - f_1_2", 
                        "max(100*f_1, 10)", 
                        "max(100*f_2, 10)", 
                        "max((200*(f_1 + f_2) + 50*f_1_2), 1)")
    suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                             frequencyDependentFitness = TRUE))
    null <- capture.output({
    set.seed(1)
    of1 <- oncoSimulIndiv(afe, 
                          model = "McFL",
                          initMutant = "WT",
                          onlyCancer = FALSE, 
                          finalTime = 50, 
                          verbosity = 0, 
                          mu = 1e-6,
                          initSize = 500, 
                          keepPhylog = FALSE,
                          seed = NULL, 
                          errorHitMaxTries = TRUE, 
                          errorHitWallTime = TRUE)

    set.seed(1)
    of2 <- oncoSimulIndiv(afe, 
                          model = "McFL",
                          onlyCancer = FALSE, 
                          finalTime = 50, 
                          verbosity = 0, 
                          mu = 1e-6,
                          initSize = 500, 
                          keepPhylog = FALSE,
                          seed = NULL, 
                          errorHitMaxTries = TRUE, 
                          errorHitWallTime = TRUE)
    })
    expect_true(of1$InitMutant == "WT")
    expect_true(of2$InitMutant == "")
    expect_identical(of1[!(names(of1) == "InitMutant")],
                     of2[!(names(of2) == "InitMutant")])

    null <- capture.output({
    set.seed(2)
    o2 <- allFitnessEffects(genotFitness = rfitness(2))
    set.seed(2)
    os1 <- oncoSimulIndiv(o2, initMutant = "WT",
                   model = "McFLD",
                   onlyCancer = FALSE, 
                   finalTime = 50, 
                   verbosity = 0, 
                   mu = 1e-6,
                   initSize = 500, 
                   keepPhylog = FALSE,
                   seed = NULL, 
                   errorHitMaxTries = TRUE, 
                   errorHitWallTime = TRUE)
    set.seed(2)
    os2 <- oncoSimulIndiv(o2, 
                   model = "McFLD",
                   onlyCancer = FALSE, 
                   finalTime = 50, 
                   verbosity = 0, 
                   mu = 1e-6,
                   initSize = 500, 
                   keepPhylog = FALSE,
                   seed = NULL, 
                   errorHitMaxTries = TRUE, 
                   errorHitWallTime = TRUE)
    })
    expect_true(os1$InitMutant == "WT")
    expect_true(os2$InitMutant == "")
    expect_identical(os1[!(names(of1) == "InitMutant")],
                     os2[!(names(of2) == "InitMutant")])
    
    
})

test_that("initMutant does not accept integer vectors", {
    ## Yes, evalGenotype does accept them
    set.seed(2)
    o2 <- allFitnessEffects(genotFitness = rfitness(2))
    set.seed(2)
    expect_error(oncoSimulIndiv(o2, initMutant = 1,
                   model = "McFLD",
                   onlyCancer = FALSE, 
                   finalTime = 50, 
                   verbosity = 0, 
                   mu = 1e-6,
                   initSize = 500, 
                   keepPhylog = FALSE,
                   seed = NULL, 
                   errorHitMaxTries = TRUE, 
                   errorHitWallTime = TRUE),
                 "initMutant must be a (list of) character string(s)",
                 fixed = TRUE)
})


## zz4:
test_that("initMutant: multiple pops, basic", {
    o1 <- allFitnessEffects(
        noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))
    
    expect_silent(out <- oncoSimulIndiv(o1, initMutant = c("c, b", "b"),
                                        initSize = c(300, 20),
                                        onlyCancer = FALSE,
                                        seed = NULL))

    null <- capture.output({
    outx <- oncoSimulIndiv(o1, initMutant = c("c, b", "WT"),
                   initSize = c(300, 20),
                   onlyCancer = FALSE,
                   seed = NULL)

    outx <- oncoSimulIndiv(o1, initMutant = c("c, b", "WT", "a"),
                   initSize = c(300, 20, 10),
                   onlyCancer = FALSE,
                   seed = NULL)
    })
    
    ## Pass all genotypes
    null <- capture.output(oncoSimulIndiv(o1, initMutant = c("WT", "a", "b", "c",
                                      "a, b", "a, c", "b, c",
                                      "a, b, c"),
                   initSize = c(300, 20, 10, 6, 9, 5, 4, 6),
                   onlyCancer = FALSE,
                   seed = NULL))

    set.seed(1)
    r2 <- rfitness(2)
    r2[4, 3] <- 0
    o2 <- allFitnessEffects(genotFitness = r2)
    ## oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
    ##                initSize = c(300, 20),
    ##                onlyCancer = FALSE)
    expect_warning( oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
                   initSize = c(300, 20),
                   onlyCancer = FALSE),
                   "Init Mutant with birth == 0.0", fixed = TRUE)
    set.seed(1)
    r2 <- rfitness(2)
    o2 <- allFitnessEffects(genotFitness = r2)
    expect_silent(out <- oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
                   initSize = c(300, 20),
                   onlyCancer = FALSE,
                   seed = NULL))
})


test_that("lengths initSize and initMutant must match", {
    set.seed(1)
    r2 <- rfitness(2)
    o2 <- allFitnessEffects(genotFitness = r2)
    expect_error(oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
                          initSize = c(20),
                          onlyCancer = FALSE),
                 "Lengths of initSize and initMutant differ",
                 fixed = TRUE)
    })


test_that("multiple init mutants: pops of 0 size is error", {
    set.seed(1)
    r2 <- rfitness(2)
    o2 <- allFitnessEffects(genotFitness = r2)
    expect_error(oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
                          initSize = c(0, 20),
                          onlyCancer = FALSE),
                 "At least one initSize <= 0",
                 fixed = TRUE)
    ## Now caught in R
    ## expect_true(out$other$ExceptionMessage
    ##              "Unrecoverable exception: ti: popSize <= 0",
    ##              fixed = TRUE)
})


## with fdf and without fdf
test_that("multiple init mutants: different species", {

    num_reps <- 10
    for(i in 1:num_reps) {
        r2 <- rfitness(6)
        ## Make sure these always viable for interesting stuff
        r2[2, 7] <- 1 + runif(1) # A
        r2[4, 7] <- 1 + runif(1) # C
        r2[8, 7] <- 1 + runif(1) # A, B
        ## Two species, first with pos 1 and 2
        ## all with both mutated is 0
        r2[ (r2[, 1] == 0) & (r2[, 2] == 1),  7] <- 0
        r2[ (r2[, 3] == 0) & apply(r2[, 3:6] == 1, 1, any),  7] <- 0
        
        r2[ (r2[, 1] == 1) & apply(r2[, 3:6] == 1, 1, any),  7] <- 0
        r2[ (r2[, 3] == 1) & apply(r2[, 1:2] == 1, 1, any),  7] <- 0
        r2 <- r2[r2[, 7] > 0, ]
        o2 <- allFitnessEffects(genotFitness = r2)
        ag <- evalAllGenotypes(o2)
        ## Set mutation rate of the "species indicator" close to 1e-10
        out1 <- oncoSimulIndiv(o2, initMutant = c("A", "C"),
                               initSize = c(1e3, 1e4),
                               finalTime = 300,
                               mu = c(A = 1e-10, C = 1e-10,
                                      B = 1e-5, D = 1e-5, E = 1e-5, F = 1e-5),
                               onlyCancer = FALSE)
        not_possible <- c("", ag$Genotype[ag$Fitness == 0])
        expect_false(any(not_possible %in%  out1$GenotypesLabels))
    }


    num_reps <- 5
    for(i in 1:num_reps) {
        r2 <- rfitness(6)
        ## Make sure these always viable for interesting stuff
        r2[2, 7] <- 1 + runif(1) # A
        r2[4, 7] <- 1 + runif(1) # C
        r2[8, 7] <- 1 + runif(1) # A, B
        ## Two species, first with pos 1 and 2
        ## all with both mutated is 0
        r2[ (r2[, 1] == 0) & (r2[, 2] == 1),  7] <- 0
        r2[ (r2[, 3] == 0) & apply(r2[, 3:6] == 1, 1, any),  7] <- 0
        
        r2[ (r2[, 1] == 1) & apply(r2[, 3:6] == 1, 1, any),  7] <- 0
        r2[ (r2[, 3] == 1) & apply(r2[, 1:2] == 1, 1, any),  7] <- 0
        r2 <- r2[r2[, 7] > 0, ]
        o2 <- allFitnessEffects(genotFitness = r2)
        ag <- evalAllGenotypes(o2)
        ## Set mutation rate of the "species indicator" close to 1e-10
        ## But this ain't really enough: what about A, D? A is mutating
        ## to A, D with rate given by D, not by C.
        out1 <- oncoSimulIndiv(o2, initMutant = c("A", "C"),
                   model = "McFLD",
                   initSize = c(1e3, 1e4),
                   finalTime = 500,
                   mu = c(A = 1e-10, C = 1e-10,
                          B = 1e-5, D = 1e-5, E = 1e-5, F = 1e-5),
                   onlyCancer = FALSE)
        not_possible <- c("", ag$Genotype[ag$Fitness == 0])
        expect_false(any(not_possible %in%  out1$GenotypesLabels))
    }
})




test_that("multiple init mutants: cannot have descendants of absent parents", {
    num_reps <- 10
    for(i in 1:num_reps) {
        r2 <- rfitness(6)
        ## Make sure these always viable 
        r2[8, 7] <- 1 + runif(1) # A, B
        r2[19, 7] <- 1 + runif(1) # C, F
        o2 <- allFitnessEffects(genotFitness = r2)
        ag <- evalAllGenotypes(o2)
        ## Set mutation rate of the "species indicator" close to 1e-10
        out1 <- oncoSimulIndiv(o2, initMutant = c("A, B", "C, F"),
                               initSize = c(100, 200),
                               onlyCancer = FALSE)
        not_possible <- unique(c("", LETTERS[1:6],
                          paste0("A, ", LETTERS[3:6]),
                          paste0("B, ", LETTERS[3:6]),
                          paste0("C, ", LETTERS[4:5]),
                          paste0("D, ", LETTERS[5:6]),
                          "E, F",
                          ag$Genotype[ag$Fitness == 0]))
        expect_false(any(not_possible %in%  out1$GenotypesLabels))
    }
})


test_that("multiple init mutants: different species, FDF", {
    gffd0 <- data.frame(
        Genotype = c(
                     "A", "A, B",
                     "C", "C, D", "C, E"),
        Fitness = c(
                    "1.3",
                    "1.4",
                    "1.4",
                    "1.1 + 0.7*((f_1 + f_A_B) > 0.3)",
                    "1.2 + sqrt(f_A + f_C + f_C_D)"),        
    stringsAsFactors = FALSE)
    suppressWarnings(afd0 <- allFitnessEffects(genotFitness = gffd0,
                             frequencyDependentFitness = TRUE))

    suppressWarnings(eag0 <- evalAllGenotypes(afd0, spPopSizes = 1:5))


    
    gffd <- data.frame(
        Genotype = c("WT",
                     "A", "A, B",
                     "C", "C, D", "C, E"),
        Fitness = c("0",
                    "1.3",
                    "1.4",
                    "1.4",
                    "1.1 + 0.7*((f_1 + f_1_2) > 0.3)",
                    "1.2 + sqrt(f_1 + f_3 + f_3_4)"),        
    stringsAsFactors = FALSE)
    suppressWarnings(afd <- allFitnessEffects(genotFitness = gffd, 
                             frequencyDependentFitness = TRUE))

    suppressWarnings(eag1 <- evalAllGenotypes(afd, spPopSizes = 0:5))

    ## No wildtype
    gffd2 <- data.frame(
        Genotype = c("A", "A, B",
                     "C", "C, D", "C, E"),
        Fitness = c("1.3",
                    "1.4",
                    "1.4",
                    "1.1 + 0.7*((f_1 + f_1_2) > 0.3)",
                    "1.2 + sqrt(f_1 + f_3 + f_3_4)"),        
    stringsAsFactors = FALSE)
    suppressWarnings(afd2 <- allFitnessEffects(genotFitness = gffd2, 
                             frequencyDependentFitness = TRUE))

    suppressWarnings(eag2 <- evalAllGenotypes(afd2, spPopSizes = 1:5))
    expect_identical(eag1, eag2)
    expect_identical(eag0, eag2)

    set.seed(1)
    os1 <- oncoSimulIndiv(afd, initMutant = "A", seed = NULL,
                          finalTime = 20, initSize = 1e6,                         
                          onlyCancer = FALSE, model = "McFLD")
    set.seed(1)
    os2 <- oncoSimulIndiv(afd2, initMutant = "A", seed = NULL,
                          finalTime = 20, initSize = 1e6,
                          onlyCancer = FALSE, model = "McFLD")
    set.seed(1)
    os0 <- oncoSimulIndiv(afd0, initMutant = "A", seed = NULL,
                          finalTime = 20, initSize = 1e6,
                          onlyCancer = FALSE, model = "McFLD")

    expect_equivalent(os1, os2)
    expect_equivalent(os1, os0)
})



test_that("multiple init mutants: different species, FDF, check fitness", {
    mspecF <- data.frame(
        Genotype = c("A",
                     "A, a1", "A, a2", "A, a1, a2",
                     "B",
                     "B, b1", "B, b2", "B, b3",
                     "B, b1, b2", "B, b1, b3", "B, b1, b2, b3"),
        Fitness = c("1 + f_A_a1",
                    "1 + f_A_a2",
                    "1 + f_A_a1_a2",
                    "1 + f_B",
                    "1 + f_B_b1",
                    "1 + f_B_b2",
                    "1 + f_B_b3",
                    "1 + f_B_b1_b2",
                    "1 + f_B_b1_b3",
                    "1 + f_B_b1_b2_b3",
                    "1 + f_A")
    )
    suppressWarnings(fmspecF <- allFitnessEffects(genotFitness = mspecF,
                                 frequencyDependentFitness = TRUE))
    ## Remeber, spPopSizes correspond to the genotypes
    ## shown in
    fmspecF$full_FDF_spec
    ## in exactly that order

    suppressWarnings(afmspecF <- evalAllGenotypes(fmspecF, spPopSizes = 1:11))

    ## Show only viable ones
    afmspecF[afmspecF$Fitness >= 1, ]

    ## Expected values
    exv <- 1 + c(3, 5, 4, 8, 6, 7, 9, 2, 10, 11, 1)/sum(1:11)
    stopifnot(isTRUE(all.equal(exv, afmspecF[afmspecF$Fitness >= 1, ]$Fitness)))


})


test_that("multiple init mutants: different species, FDF, exprtk crash if not in fitness table", {
    ## Crash, as f_2 is not present
    gffd <- data.frame(
        Genotype = c("WT",
                     "A", "A, B",
                     "C", "C, D", "C, E"),
        Fitness = c("0",
                    "1.3",
                    "1.4",
                    "1.4",
                    "1.1 + 0.7*((f_1 + f_2) > 0.3)",
                    "1.2 + sqrt(f_1 + f_3 + f_2)"),        
    stringsAsFactors = FALSE)
    suppressWarnings(afd <- allFitnessEffects(genotFitness = gffd, 
                             frequencyDependentFitness = TRUE))
							 
    suppressWarnings(expect_error(evalAllGenotypes(afd, spPopSizes = rep(10, 6))))
    ### FIXME: catch this exact string"Undefined symbol: 'f_2'", fixed = TRUE)
})

set.seed(NULL)
cat(paste("\n Ending init-mutant tests", date(), "\n"))


cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)



 








