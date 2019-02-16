## This is a short version of test.fitness-preds-long.R But it encompasses
## the most complex cases. It is expected that this will fail occasionally
## (based on p-values); thus, we loop twice if needed.


## Here we compare against expected population size for different fitness
## specifications 

## FIXME: We might do equivalence testing, here and in its long
## version. But since we have huge sample sizes, this is not bad. And we'd
## need to be more careful about distributions (what transformation, if
## any, etc).

cat(paste("\n Starting fitness preds at", date(), "\n"))
## RNGkind("Mersenne-Twister") ## but this is irrelevant now.


## rm(list = ls())


expe <- function(no, s, ft) {
    no * exp(((1 + s) - 1) * ft)
}

mce <- function(s, K) {
    ## yes, at equilibrium
    return( K * (exp(1 + s) - 1))
}


expected <- function(no, s, ft, model) {
    if(model == "Exp")
        return(expe(no = no, s = s, ft = ft))
    if(model == "McFL") ## this K has to be the same as in OandE, of course.        
        return(mce(s = s, K = no/(exp(1) - 1) ))
}

OandE <- function(fe, s, ft,  model, initMutant, no,
                  reps, mu, verbose = FALSE,
                  sampleEvery = sEvery) {
    ## sampleEevery can be large, as we stop on ft.
    ## But not too large, to avoid numerical issues
    ## with large s
    E <- expected(no = no, s = s, ft = ft, model = model)
    O <- oncoSimulPop(reps,
                      fe,
                      mu = mu,
                      initSize = no,
                      ## K = no,
                      model = model,
                      detectionProb = NA,
                      detectionDrivers = 99,
                      finalTime = ft,
                      detectionSize = 1e12,
                      sampleEvery = sampleEvery,
                      keepEvery = ft,
                      initMutant = initMutant,
                      mutationPropGrowth = FALSE,
                      onlyCancer = FALSE,
                      mc.cores = 2)
    if(verbose) {
        print(E)
        print(summary(O)[, c(1:3, 8:9)])
    }
    return(c(E,
             unlist(lapply(O, function(x) x$TotalPopSize))))
}


## A comment about the McFL model and what we do in OandE

##  We set K as given by default, so initSize/(exp(1) - 1).  This does not
##  have a major effect iff you run the thing for long enough. Why?
##  Because if you set K to, say, initSize, then, especially with the s =
##  0, you need to run it for long for it to reach equilibrium. In
##  addition, I am testing here in the most general, usual, circumstances,
##  given that especial ones (other values of K) work if you use long
##  enough finalTimes.


verboseOandE <- FALSE
sEvery <- 0.05

## we could now use try_again
date()
test_that("Observed vs expected, case III", {
    
    cat("\n Observed vs expected, case III\n")
    genmodule <- function(l, num) {
        paste(paste0(l, 1:num), collapse = ", ")
    }
    max.tries <- 5 ## yes, like 1 over 10000 or 20000 we get up to 4. We do lots of p-value tests.
    for(tries in 1:max.tries) {
        ## Like II, but we combine several effects.
        reps <- 25
        mu <- 1e-7
        nig <- 20 ## there are lots of genes with modules with many genes
        out <- NULL
        outNO <- NULL
        so <- 0.2
        sni <- 0.1
        niG <- c(sni, rep(0, nig - 1))
        names(niG) <- paste0("nint", 1:nig)
        feo <- allFitnessEffects(orderEffects = c("A > B" = so,
                                                  "B > A" = -5),
                                 noIntGenes = niG,
                                 geneToModule = c(A = genmodule("a", 19),
                                                  B = genmodule("b", 15)))
        so1 <- so + sni + so * sni
        out <- rbind(out, OandE(feo, so1, 10, "Exp", "a1 > b2 > nint1", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 30, "McFL", "a2 > nint1 > b8", 1.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 15, "Exp", "a9 > b3 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 40, "McFL", "nint3 > nint1 > a3 > b2", 1.4e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 20, "Exp", "a4 > nint9 > nint1 > b1", 1e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 40, "McFL", "nint1 > nint8 > a7 > b2", 3.1e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 12, "Exp", "a2 > b13 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 50, "McFL", "a19 > b10 > nint6 > nint1", 3.1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 12, "Exp", "a2 > b13", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 50, "McFL", "a19 > b10 > nint6", 1.9e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 12, "Exp", "a2 > nint4 > b13", 5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 50, "McFL", "a19 > b10 > nint6 > nint8", 7.1e4, reps, mu, verboseOandE))
        so <- 0.17
        niG <- c(sni, rep(0, nig - 1))
        names(niG) <- paste0("nint", 1:nig)
        so1 <- so + sni + so * sni
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, so),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 geneToModule = c(
                                     "Root" = "Root",
                                     A = genmodule("a", 10),
                                     B = genmodule("b", 5),
                                     C = genmodule("c", 7)),
                                 noIntGenes = niG)
        out <- rbind(out, OandE(feo, so1, 29, "Exp", "a2, c3, nint1", 1.8e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 40, "McFL", "a5, nint1, c1", 1.6e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 10, "Exp", "nint1, a2, c3", 1.7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 36, "McFL", "nint6, nint1, a5, c1", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 32, "Exp", "nint1 > a2 > c3", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 46, "McFL", "nint6 > a5 > c1 > nint1", 4.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 18, "Exp", "a2 > nint2 > nint1 > c3", 3.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 46, "McFL", "a5 > nint9 > c1 > nint1", 1.97e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 20, "Exp", "c1, a7, c2, a3, a9, a8, nint1, nint7", 1.76e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 40, "McFL", "a5, nint1, c6, a2, nint3, c5", 3.3e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 40, "Exp", "a2, c6", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 40, "McFL", "b4, c2", 1e3, reps, mu, verboseOandE))
        s <- rep(0, 4)
        s[1] <- .05
        s[2] <- .07
        s[3] <- .08
        s[4] <- .09
        st <- prod(1 + s) - 1
        niG <- c(s[4], rep(0, nig - 1))
        names(niG) <- paste0("nint", 1:nig)
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, s[1]),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 orderEffects = c("F > E" = s[2],
                                                  "E > F" = -5),
                                 epistasis = c("G : H" = s[3]),
                                 geneToModule = c(
                                     "Root" = "Root",
                                     A = genmodule("a", 8),
                                     B = genmodule("b", 9),
                                     C = genmodule("c", 4),
                                     E = genmodule("e", 5),
                                     F = genmodule("f", 6),
                                     G = genmodule("g", 7),
                                     H = genmodule("h", 4)
                                 ),
                                 noIntGenes = niG)
        out <- rbind(out, OandE(feo, st, 15, "Exp",
                                "a2 > f3 > e2 > c1 > g4 > h2 > nint1", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 15, "Exp",
                                "a2 > a4 > f3 > e2 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 15, "Exp",
                                "a6 > f1 > e3 > g4 > c1 > nint1 > h2 > nint3", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 20, "Exp",
                                "f3 > e2 > c1 > g4 > h2 > nint1 > a4", 4.8e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 30, "Exp",
                                "h2 > f3 > e2 > c1 > g4 > a3 > nint1 > nint4", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 28, "Exp",
                                "a7 > f3 > e2 > c1 > g4 > a6 > h2 > nint1 > nint3", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 45, "McFL",
                                "a2 > f3 > e2 > c1 > g4 > h2 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 35, "McFL",
                                "a2 > a4 > f3 > e2 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4.7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 35, "McFL",
                                "a6 > f1 > e3 > g4 > c1 > nint1 > h2 > nint3", 5.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 30, "McFL",
                                "f3 > e2 > c1 > g4 > h2 > nint1 > a4", 4.21e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 40, "McFL",
                                "h2 > f3 > e2 > c1 > g4 > a3 > nint1 > nint4", 5.1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 38, "McFL",
                                "a7 > f3 > e2 > c1 > g4 > a6 > h2 > nint1 > nint3", 1.16e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 15, "Exp",
                                    "f3 > e2 > c1 > g4 > h2 > nint1", 2.5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 15, "Exp",
                                    "a2 > a4 > e3 > f2 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4.3e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 15, "Exp",
                                    "a6 > f1 > e3 > nint1 > g5 > h2 > nint3", 5.2e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 20, "Exp",
                                    "f3 > e2 > c1 > h2 > nint1 > a4", 4.1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 30, "Exp",
                                    "h2 > e2 > c1 > g4 > a3 > nint1 > nint4", 5e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 28, "Exp",
                                    "a7 > f3 > c1 > g4 > a6 > h2 > nint1 > nint3", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 45, "McFL",
                                    "a2 > e2 > c1 > g4 > h2 > nint1", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 35, "McFL",
                                    "a2 > a4 > e3 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 35, "McFL",
                                    "a6 > f1 > e3 > c1 > nint1 > h2 > nint3", 5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 30, "McFL",
                                    "f3 > e2 > c1 > g4 > nint1 > a4", 4e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 40, "McFL",
                                    "h2 > e2 > c1 > g4 > a3 > nint1 > nint4", 5e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 38, "McFL",
                                    "a7 > f3 > c1 > g4 > a6 > h2 > nint1 > nint3", 1e3, reps, mu, verboseOandE))
        out <- data.frame(out)
        colnames(out) <- c("Expected", paste0("Observed_", 1:reps))
        d1 <- data.frame(Expected = out[, 1], Observed = rowMeans(out[, -1]))
        p.fail <- 0.01
        lm1 <- lm(log(Observed) ~ log(Expected), data = d1)
        
        ## For NO, do a t.test by row.
        no.t <- apply(outNO, 1,
                      function(x) t.test(log(x[2:(reps + 1)] + 1), mu = log(x[1] + 1))$p.value
                      )
        ## And repeat by row for the ones expected OK
        yes.t <- apply(out, 1,
                       function(x) t.test(log(x[2:(reps + 1)] + 1), mu = log(x[1] + 1))$p.value
                       )
        p.value.threshold <- 1e-6

        T.not <- (all(no.t < p.value.threshold))
        T.yest <- (min(p.adjust(yes.t, method = "BH")) > p.fail)
        T.lm <- (car::linearHypothesis(lm1, diag(2), c(0, 1))[["Pr(>F)"]][2] >
                 p.fail)
        ## so a difference of 0.1 would be large enough, and this would
        ## fail even if large sd in estimates
        T.lm.diff <- all(abs(coefficients(lm1) - c(0, 1) ) < 0.1)
        if( T.not && T.yest && T.lm && T.lm.diff ) break;
        if(! (T.not && T.yest && T.lm && T.lm.diff) ) {
            cat("\n T.not is \n"); print(T.not)
            print(no.t)
            cat("\n T.yest \n"); print(T.yest)
            print(yes.t)
            cat("\n larger values of T.yest \n")
            print(which(p.adjust(yes.t, method = "BH") <= p.fail))
            cat("\n lm and plot\n")
            print(summary(lm1)); print(car::linearHypothesis(lm1, diag(2), c(0, 1)))
            plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        }
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true((T.not && T.yest && T.lm && T.lm.diff) )
})
date()

cat(paste("\n Ending fitness preds long at", date(), "\n"))
