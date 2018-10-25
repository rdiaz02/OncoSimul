### Benchmark runs. This produces "benchmark_1.RData"


rm(list = ls())
set.seed(NULL)

library(OncoSimulR)

######################################################################
######################################################################
######################################################################
######################################################################
system_summary <- function() {
    return(list(versioninfo = version,
                memimfo = system("free", intern = TRUE),
                cpuinfo = system("cat /proc/cpuinfo | grep 'model name'", intern = TRUE),
                packageinfo = paste("OncoSimulR, ", packageVersion("OncoSimulR")),
##                nodeinfo = Sys.info()$nodename,
##                nodelinuxinfo = paste(Sys.info()$sysname, Sys.info()$release),
                dateinfo = date()))
}

stats_simuls <- function(sim) {
    ## sim is an oncoSimulPop output
    trf <- function(x) {
        tt <- try(c(NumClones = x$NumClones,
                    NumIter = x$NumIter,
                    FinalTime = x$FinalTime,
                    TotalPopSize = x$TotalPopSize,
                    Attempts = x$other$attemptsUsed))
        if(!inherits(tt, "try-error")) {
            return(tt)
        } else {
            return(c(NumClones = NA,
                    NumIter = NA,
                    FinalTime = NA,
                    TotalPopSize = NA,
                    Attempts = NA))
        }
    }
    tmp <- try(do.call("rbind", lapply(sim, trf)))
    unlist(lapply(as.data.frame(tmp), summary))
}

all_sim_stats <- function(...) {
    ## Returns a single data frame with all benchmark info
    ## in ... pass the names of the objects
    names <- as.character(as.list(match.call())[-c(1)])
    m1 <- do.call("rbind", lapply(list(...), stats_simuls))
    ## the function will break next if something failed in previous
    rownames(m1) <- names
    Numindiv <- unlist(lapply(list(...), length))
    time <- unlist(lapply(paste0("t_", names), get))
    size <- unlist(lapply(list(...), object.size))
    df <- data.frame(m1, size, time, Numindiv)
    df$time_per_simul <- df$time/df$Numindiv
    df$size_mb_per_simul <- df$size/(df$Numindiv * 1024^2)
    attributes(df)$system_summary <- system_summary()
    return(df)
}


all_sim_stats_single <- function(sim, tsim, name) {
    ## Like previous, but for a single simulation. We avoid the dangerous
    ## dynGet that we would nedd to obtain the "t_"

    ## Returns a single data
    ## frame with all benchmark info in ... pass the names of the objects
    ## names <- as.character(as.list(match.call())[-c(1)])
    names <- name
    m1 <- data.frame(as.list(stats_simuls(sim)))
    ## the function will break next if something failed in previous
    rownames(m1) <- names
    Numindiv <- length(sim)
    time <- tsim
    size <- as.numeric(object.size(sim))
    df <- data.frame(m1, size, time, Numindiv)
    df$time_per_simul <- df$time/df$Numindiv
    df$size_mb_per_simul <- df$size/(df$Numindiv * 1024^2)
    df$name <- names
    ## attributes(df)$system_summary <- system_summary()
    return(df)
}


######################################################################
######################################################################
######################################################################
######################################################################
####################                             #####################
####################   Fitness specifications    #####################
####################                             #####################
######################################################################
######################################################################
######################################################################
######################################################################


pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"),
    drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53", 
	             "MLL3", "TGFBR2", "PXDN"))


## Random fitness landscape with 6 genes 
rfl6 <- rfitness(6, min_accessible_genotypes = 50)
attributes(rfl6)$accessible_genotypes ## How many actually accessible
rf6 <- allFitnessEffects(genotFitness = rfl6)


## Random fitness landscape with 12 genes
rfl12 <- rfitness(12, min_accessible_genotypes = 200)
attributes(rfl12)$accessible_genotypes ## How many actually accessible
rf12 <- allFitnessEffects(genotFitness = rfl12)




## Independent genes; positive fitness from exponential distribution
## mean around 0.1, and negative from exponential with mean around 0.02.
## Half positive, half negative
ng <- 200
re_200 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10), -rexp(ng/2, 50)))

ng <- 500
re_500 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10), -rexp(ng/2, 50)))

ng <- 2000
re_2000 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10), -rexp(ng/2, 50)))

ng <- 4000
re_4000 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10), -rexp(ng/2, 50)))



cases <- expand.grid(Model = c("Exp", "McFL"),
                     fitness = c("pancr", "rf6", "rf12", 
                                 "re_200", "re_500", "re_2000", "re_4000"),
                     stringsAsFactors = FALSE)



sim_bench <- function(Model, fitness, Nindiv, ...) {
    cat("\n\n\n")
    cat("***  Doing ", Model, " ", fitness)
    
    t_tmp <- system.time({
        if(Model == "Exp") {
            tmp <- oncoSimulPop(Nindiv,
                                get(fitness),
                                detectionProb = NA, 
                                detectionSize = 1e5,
                                initSize = 500,
                                detectionDrivers = NA,
                                keepPhylog = TRUE,
                                model = "Exp",
                                errorHitWallTime = FALSE,
                                errorHitMaxTries = FALSE,
                                finalTime = 25000,
                                onlyCancer = TRUE,
                                mc.cores = 1,
                                sampleEvery = 0.5,
                                ...)
        } else {
            initSize <- 1000
            tmp <- oncoSimulPop(Nindiv,
                                get(fitness),
                                detectionProb = c(
                                    PDBaseline = 1.4 * initSize,
                                    n2 = 2 * initSize,
                                    p2 = 0.1,
                                    checkSizePEvery = 4),
                                initSize = initSize,
                                detectionSize = NA,
                                detectionDrivers = NA,
                                keepPhylog = TRUE,
                                model = "McFL",
                                errorHitWallTime = FALSE,
                                errorHitMaxTries = FALSE,
                                finalTime = 25000,
                                max.wall.time = 10,
                                onlyCancer = TRUE,
                                mc.cores = 1,
                                ...)
        }
    })["elapsed"]

    cat("\n\n\n t_tmp = ", t_tmp, "\n")
    print(object.size(tmp)/(Nindiv * 1024^2))
    cat("\n\n")
    print(summary(unlist(lapply(tmp, "[[", "NumClones"))))
    print(summary(unlist(lapply(tmp, "[[", "NumIter"))))
    print(summary(unlist(lapply(tmp, "[[", "FinalTime"))))
    print(summary(unlist(lapply(tmp, "[[", "TotalPopSize"))))
    name <- paste(Model, fitness, sep = "_")
    df <- all_sim_stats_single(tmp, t_tmp, name)
    df$Model <- Model
    df$fitness <- fitness
    return(df)
}



Nindiv <- 100

benchmark_3 <- dplyr::bind_rows(Map(sim_bench,
                             cases[, 1], cases[, 2], Nindiv,
                             keepEvery = 1))

attributes(benchmark_3)$system_summary <- system_summary()

save(file = "../../data/benchmark_3.RData", benchmark_3)


gc()
