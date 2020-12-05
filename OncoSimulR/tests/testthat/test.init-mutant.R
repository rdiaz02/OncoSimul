inittime <- Sys.time()
cat(paste("\n Starting init-mutant tests", date(), "\n"))

## RNGkind("Mersenne-Twister")

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
    cat(paste("\n done tries", tries, "\n"))
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
    cat(paste("\n done tries", tries, "\n"))
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
    cat(paste("\n done tries", tries, "\n"))
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
    for(i in 1:10) {
        set.seed(i)
        o2i <- allFitnessEffects(genotFitness = rfitness(2))
        o2io <- oncoSimulIndiv(o2i, initMutant = "B, A")
        if(!is.null(o2io$other$ExceptionMessage)) {
            expect_false(
                grepl("Trying to obtain a mutation when nonmutated.size is 0",
                      o2io$other$ExceptionMessage)
            )
        }
    }
})


## ## zz4:
## test_that("initMutant: multiple pops", {

##     o1 <- allFitnessEffects(
##         noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))

##     ## zz: test with the above fitness spec
##     o2 <- allFitnessEffects(genotFitness = rfitness(2))

##     oncoSimulIndiv(o2, initMutant = c("B, A", "A"),
##                    initSize = c(300, 20),
##                    onlyCancer = FALSE)

    
## })



test_that("initMutant with freq-dep-fitness"  , {
    r <- data.frame(rfitness(2))
    r[, "Fitness"] <- c("f_ - f_1 - f_2 - f_1_2", 
                        "max(100*f_1, 10)", 
                        "max(100*f_2, 10)", 
                        "max((200*(f_1 + f_2) + 50*f_1_2), 1)")
    afe <- allFitnessEffects(genotFitness = r, 
                             frequencyDependentFitness = TRUE)

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


cat(paste("\n Ending init-mutant tests", date(), "\n"))




cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)



