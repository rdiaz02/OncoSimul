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
               s = 0.05,
               sh = -0.9,
               typeDep = "MN"),
    drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53", 
	             "MLL3", "TGFBR2", "PXDN"))


Nindiv <- 100


## keepEvery = 1
t_exp1 <- system.time(
    exp1 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = "default", 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = 1,
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]


t_mc1 <- system.time(
    mc1 <- oncoSimulPop(Nindiv, pancr, 
                           detectionProb = "default", 
                           detectionSize = NA,
                           detectionDrivers = NA,
                           finalTime = NA,
                           keepPhylog = TRUE, keepEvery = 1,                                  
                           model = "McFL", 
                           mc.cores = 1))["elapsed"]

## keepEvery = NA
t_exp2 <- system.time(
    exp2 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = "default", 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]


t_mc2 <- system.time(
    mc2 <- oncoSimulPop(Nindiv, pancr, 
                           detectionProb = "default", 
                           detectionSize = NA,
                           detectionDrivers = NA,
                           finalTime = NA,
                           keepPhylog = TRUE, keepEvery = NA,
                           model = "McFL", 
                           mc.cores = 1))["elapsed"]



### exp3 to exp6
t_exp3 <- system.time(
    exp3 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e4,
                                              p2 = 0.1, n2 = 5e5,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = 1, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]

t_exp4 <- system.time(
    exp4 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e4,
                                              p2 = 0.1, n2 = 5e5,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]



t_exp5 <- system.time(
    exp5 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e5,
                                              p2 = 0.1, n2 = 5e7,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = 1, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]

t_exp6 <- system.time(
    exp6 <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e5,
                                              p2 = 0.1, n2 = 5e7,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]


####################################################
###### onlyCancer = FALSE


## keepEvery = 1
t_exp1_noc <- system.time(
    exp1_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = "default", 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = 1,
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]
t_mc1_noc <- system.time(
    mc1_noc <- oncoSimulPop(Nindiv, pancr, 
                           detectionProb = "default", 
                           detectionSize = NA,
                           detectionDrivers = NA,
                           finalTime = NA, onlyCancer = FALSE,
                           keepPhylog = TRUE, keepEvery = 1,                                  
                           model = "McFL", 
                           mc.cores = 1))["elapsed"]
## keepEvery = NA
t_exp2_noc <- system.time(
    exp2_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = "default", 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]
t_mc2_noc <- system.time(
    mc2_noc <- oncoSimulPop(Nindiv, pancr, 
                           detectionProb = "default", 
                           detectionSize = NA,
                           detectionDrivers = NA,
                           finalTime = NA, onlyCancer = FALSE,
                           keepPhylog = TRUE, keepEvery = NA,
                           model = "McFL", 
                           mc.cores = 1))["elapsed"]
### exp3_noc to exp6_noc
t_exp3_noc <- system.time(
    exp3_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e4,
                                              p2 = 0.1, n2 = 5e5,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = 1, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]
t_exp4_noc <- system.time(
    exp4_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e4,
                                              p2 = 0.1, n2 = 5e5,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]
t_exp5_noc <- system.time(
    exp5_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e5,
                                              p2 = 0.1, n2 = 5e7,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = 1, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]
t_exp6_noc <- system.time(
    exp6_noc <- oncoSimulPop(Nindiv, pancr, 
                            detectionProb = c(PDBaseline = 5e5,
                                              p2 = 0.1, n2 = 5e7,
                                              checkSizePEvery = 20), 
                            detectionSize = NA,
                            detectionDrivers = NA,
                            finalTime = NA, onlyCancer = FALSE,
                            keepPhylog = TRUE, keepEvery = NA, 
                            model = "Exp", 
                            mc.cores = 1))["elapsed"]

benchmark_1_0.05 <- all_sim_stats(exp1, mc1, exp2, mc2, exp3, exp4, exp5, exp6,
                         exp1_noc, mc1_noc, exp2_noc, mc2_noc, exp3_noc, exp4_noc, exp5_noc, exp6_noc
                         )

## Add the info about key settings
benchmark_1_0.05$keepEvery <- c(1, 1, NA, NA, 1, NA, 1, NA,
                           1, 1, NA, NA, 1, NA, 1, NA)
benchmark_1_0.05$PDBaseline <- c(rep(1.2 * 500, 4), 5e4, 5e4, 5e5, 5e5,
                            rep(1.2 * 500, 4), 5e4, 5e4, 5e5, 5e5)
benchmark_1_0.05$n2 <- c(rep(2 * 500, 4), 5e5, 5e5, 5e7, 5e7,
                    rep(2 * 500, 4), 5e5, 5e5, 5e7, 5e7)
benchmark_1_0.05$onlyCancer <- c(rep(TRUE, 8), rep(FALSE, 8))
save(file = "../../data/benchmark_1_0.05.RData", benchmark_1_0.05)


gc()

## ######################################################################
## ######################################################################
## ###############                                         ##############
## ###############      If you do it by hand ...           ##############
## ###############                                         ##############
## ######################################################################
## ######################################################################



## cat("\n\n\n t_exp1 = ", t_exp1, "\n")
## object.size(exp1)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp1, "[[", "NumClones")))
## summary(unlist(lapply(exp1, "[[", "NumIter")))
## summary(unlist(lapply(exp1, "[[", "FinalTime")))
## summary(unlist(lapply(exp1, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_mc1 = ", t_mc1, "\n")
## object.size(mc1)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(mc1, "[[", "NumClones")))
## summary(unlist(lapply(mc1, "[[", "NumIter")))
## summary(unlist(lapply(mc1, "[[", "FinalTime")))
## summary(unlist(lapply(mc1, "[[", "TotalPopSize")))

## cat("\n\n")
## cat("\n\n\n t_exp2 = ", t_exp2, "\n")
## object.size(exp2)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp2, "[[", "NumClones")))
## summary(unlist(lapply(exp2, "[[", "NumIter")))
## summary(unlist(lapply(exp2, "[[", "FinalTime")))
## summary(unlist(lapply(exp2, "[[", "TotalPopSize")))

## cat("\n\n")
## cat("\n\n\n t_mc2 = ", t_mc2, "\n")
## object.size(mc2)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(mc2, "[[", "NumClones")))
## summary(unlist(lapply(mc2, "[[", "NumIter")))
## summary(unlist(lapply(mc2, "[[", "FinalTime")))
## summary(unlist(lapply(mc2, "[[", "TotalPopSize")))



## cat("\n\n")
## cat("\n\n\n t_exp3 = ", t_exp3, "\n")
## object.size(exp3)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp3, "[[", "NumClones")))
## summary(unlist(lapply(exp3, "[[", "NumIter")))
## summary(unlist(lapply(exp3, "[[", "FinalTime")))
## summary(unlist(lapply(exp3, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp4 = ", t_exp4, "\n")
## object.size(exp4)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp4, "[[", "NumClones")))
## summary(unlist(lapply(exp4, "[[", "NumIter")))
## summary(unlist(lapply(exp4, "[[", "FinalTime")))
## summary(unlist(lapply(exp4, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp5 = ", t_exp5, "\n")
## object.size(exp5)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp5, "[[", "NumClones")))
## summary(unlist(lapply(exp5, "[[", "NumIter")))
## summary(unlist(lapply(exp5, "[[", "FinalTime")))
## summary(unlist(lapply(exp5, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp6 = ", t_exp6, "\n")
## object.size(exp6)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp6, "[[", "NumClones")))
## summary(unlist(lapply(exp6, "[[", "NumIter")))
## summary(unlist(lapply(exp6, "[[", "FinalTime")))
## summary(unlist(lapply(exp6, "[[", "TotalPopSize")))


## ## Median runs until cancer

## lapply(list(exp1, mc1, exp2, mc2, exp3, exp4, exp5, exp6),
##        function(y) median(unlist(lapply(y, function(x) x$other$attemptsUsed))))






## cat("\n\n\n t_exp1_noc = ", t_exp1_noc, "\n")
## object.size(exp1_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp1_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp1_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp1_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp1_noc, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_mc1_noc = ", t_mc1_noc, "\n")
## object.size(mc1_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(mc1_noc, "[[", "NumClones")))
## summary(unlist(lapply(mc1_noc, "[[", "NumIter")))
## summary(unlist(lapply(mc1_noc, "[[", "FinalTime")))
## summary(unlist(lapply(mc1_noc, "[[", "TotalPopSize")))

## cat("\n\n")
## cat("\n\n\n t_exp2_noc = ", t_exp2_noc, "\n")
## object.size(exp2_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp2_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp2_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp2_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp2_noc, "[[", "TotalPopSize")))

## cat("\n\n")
## cat("\n\n\n t_mc2_noc = ", t_mc2_noc, "\n")
## object.size(mc2_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(mc2_noc, "[[", "NumClones")))
## summary(unlist(lapply(mc2_noc, "[[", "NumIter")))
## summary(unlist(lapply(mc2_noc, "[[", "FinalTime")))
## summary(unlist(lapply(mc2_noc, "[[", "TotalPopSize")))



## cat("\n\n")
## cat("\n\n\n t_exp3_noc = ", t_exp3_noc, "\n")
## object.size(exp3_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp3_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp3_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp3_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp3_noc, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp4_noc = ", t_exp4_noc, "\n")
## object.size(exp4_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp4_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp4_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp4_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp4_noc, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp5_noc = ", t_exp5_noc, "\n")
## object.size(exp5_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp5_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp5_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp5_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp5_noc, "[[", "TotalPopSize")))


## cat("\n\n")
## cat("\n\n\n t_exp6_noc = ", t_exp6_noc, "\n")
## object.size(exp6_noc)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp6_noc, "[[", "NumClones")))
## summary(unlist(lapply(exp6_noc, "[[", "NumIter")))
## summary(unlist(lapply(exp6_noc, "[[", "FinalTime")))
## summary(unlist(lapply(exp6_noc, "[[", "TotalPopSize")))


## ## Median runs until cancer

## lapply(list(exp1_noc, mc1_noc, exp2_noc, mc2_noc, exp3_noc, exp4_noc, exp5_noc, exp6_noc),
##        function(y) median(unlist(lapply(y, function(x) x$other$attemptsUsed))))





## ## bench1 <- data.frame(
## ##     time = c(t_exp1, t_mc1, t_exp2, t_mc2, t_exp3, t_exp4, t_exp5, t_exp6),
## ##     size = unlist(lapply(list(exp1, mc1, exp2, mc2, exp3, exp4, exp5, exp6),
## ##                          object.size))/(Nindiv * 1024^2),
    
## ##     lapply(list(exp1, mc1, exp2, mc2, exp3, exp4, exp5, exp6),
## ##            function(y) median(unlist(lapply(y, function(x) x$other$attemptsUsed))))
## ## )




