library(OncoSimulR)
timecomp <- vector()
timev <- vector()
for (ng in 4:15){  
    for (contime in 1:20){
        r2 <- rfitness(ng)
        cat("\n", "Doing iter ", contime, "/20 for ", ng, " genes", sep="")
        source("MeasuringTimeOption_1.R")
        time1 <- system.time(s2 <- epistasis(r2))
        source("MeasuringTimeOption_2.R")
        time2 <- system.time(s2 <- epistasis(r2))
        timecomp[contime] <- 100-time2[3]*100/time1[3]                                                    
    }
    timev[ng-3] <- mean(timecomp)
}
genes <- 4:15

## smoothingSpline = smooth.spline(genes, timev, spar=0.23)
plot(genes,timev, xlab = "genes", ylab = "% faster cpp vs R")
## lines(smoothingSpline)
