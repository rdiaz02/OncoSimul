## Here we compare against expected population size for different fitness
## specifications . We also verify that minor differences are detected
## (e.g., leaving out a mutation of relatively small effect).

cat(paste("\n Starting fitness preds long at", date(), "\n"))
cat(paste("\n            a runif ", runif(1), "\n"))
## RNGkind("Mersenne-Twister") ## but this is irrelevant now.


rm(list = ls())


expe <- function(no, s, ft) {
    no * exp(((1 + s) - 1) * ft)
}

mce <- function(s, K) {
    ## yes, at equilibrium, when birth = death
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
                      detectionDrivers = 99,
                      finalTime = ft,
                      detectionSize = 1e12,
                      sampleEvery = sampleEvery,
                      keepEvery = ft,
                      initMutant = initMutant,
                      mutationPropGrowth = FALSE,
                      onlyCancer = FALSE, detectionProb = NA,
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


## Yes, the linear model is not really correct, as variance is different
## for Exp and McFL. And with Exp, variance is actually larger the samller
## the initial population size. And we could do a lack of fit test, but we
## already test each case with a t.test And I have different levels of
## variation: between simulations within params, and between params. So
## just take average of all simuls. Sure, could be done better.

verboseOandE <- FALSE
sEvery <- 0.05

date()
test_that("Observed vs expected, case I", {
    max.tries <- 5
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n Observed vs expected, case I\n")
        ## Create a set of scenarios where we know what to expect
        ## We write a small set of helper functions.
        reps <- 120
        mu <- 1e-7
        nig <- 50
        out <- NULL
        outNO <- NULL
        so <- 0.2
        feo <- allFitnessEffects(orderEffects = c("a > b" = so,
                                                  "b > a" = -5),
                                 noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a > b", 1.4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 2* 40, "McFL", "a > b", 1.4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a > b", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 2* 40, "McFL", "a > b", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a > b", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 2 * 40, "McFL", "a > b", 1.56e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "Exp", "a > b", 1.34e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 2 * 30, "McFL", "a > b", 1.27e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 5.29e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 2 * 40, "McFL", "b > a", 1.32e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 2 * 40, "McFL", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 5e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 2 * 40, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 30, "Exp", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 2 * 30, "McFL", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 400, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 500, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        so <- 0.15
        niG <- c(so, rep(0, nig))
        names(niG) <- replicate(nig + 1, paste(sample(letters, 12, replace = TRUE),
                                               collapse = ""))
        names(niG)[1] <- "ThisisA"
        feo <- allFitnessEffects(noIntGenes = niG)
        out <- rbind(out, OandE(feo, so, 40, "Exp", "ThisisA", 1.36e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 80, "McFL", "ThisisA", 1.14e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", names(nig)[2], 2.1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 80, "McFL", names(nig)[2], 2.41e3, reps, mu, verboseOandE))
        so <- 0.32
        niG <- c(so, rep(0, nig))
        names(niG) <- replicate(nig + 1, paste(sample(letters, 12, replace = TRUE),
                                               collapse = ""))
        names(niG)[1] <- "ThisisA"
        feo <- allFitnessEffects(noIntGenes = niG)
        out <- rbind(out, OandE(feo, so, 40, "Exp", "ThisisA", 1.46e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 80, "McFL", "ThisisA", 1.271e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", names(niG)[2], 3.1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 80, "McFL", names(niG)[2], 2.1e3, reps, mu, verboseOandE))
        so <- 0.17
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, so),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "A, C", 1.153e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 80, "McFL", "A, C", 1.674e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "B, C", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 80, "McFL", "B, C", 1e3, reps, mu, verboseOandE))
        so <- 0.09
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, so),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "A, C", 2.15e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 80, "McFL", "A, C", 3.14e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "B, C", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 80, "McFL", "B, C", 1e3, reps, mu, verboseOandE))
        ft <- 60
        sa02 <- 0.1
        fe02 <- allFitnessEffects(epistasis = c("A" = sa02),
                                  noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 2 * ft, "McFL", "A", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 5.4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 2 * ft, "McFL", "A", 5.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 1.6e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 2 * ft, "McFL", "A", 1.3e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 30, "Exp", "A", 1.75e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 130, "McFL", "A", 1.09e4, reps, mu, verboseOandE))
        out <- data.frame(out)
        colnames(out) <- c("Expected", paste0("Observed_", 1:reps))
        ## d1 <- tidyr::gather(out, key = "Expected", value = "Observed", 2:(reps + 1))
        d1 <- data.frame(Expected = out[, 1], Observed = rowMeans(out[, -1]))
        p.fail <- 0.05
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
        ## T.not <- (sum(no.t < p.value.threshold) > (nrow(outNO) * 0.9))
        T.not <- (all(no.t < p.value.threshold))
        T.yest <- (min(p.adjust(yes.t, method = "BH")) > p.fail)
        T.lm <- (car::linearHypothesis(lm1, diag(2), c(0, 1))[["Pr(>F)"]][2] >
                 p.fail)
        T.lm.diff <- all(abs(coefficients(lm1) - c(0, 1) ) < 0.1)
        if(! (T.not && T.yest && T.lm && T.lm.diff) ) {
            cat("\n T.not is \n"); print(T.not)
            print(no.t)
            cat("\n T.yest \n"); print(T.yest)
            print(yes.t)
            cat("\n larger values of T.yest I \n")
            print(which(p.adjust(yes.t, method = "BH") <= p.fail))
            cat("\n lm and plot\n")
            print(summary(lm1)); print(car::linearHypothesis(lm1, diag(2), c(0, 1)))
            plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        }
        if ((T.not && T.yest && T.lm && T.lm.diff) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true((T.not && T.yest && T.lm && T.lm.diff) )

})
date()

## This one will occasionally fail. I think there is a very slight
## tendency to understimate. We do not see that in the individual comparisons,
## but on the intercept of the lm: instead of 0, it seems to be sometimes around
## -0.003 or -0.004. This is of course undetectable by eye, etc, and irrelevant.

date()
test_that("Observed vs expected, case II", {
    max.tries <- 6
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n Observed vs expected, case II\n")
        genmodule <- function(l, num) {
            paste(paste0(l, 1:num), collapse = ", ")
        }
        ## Create a set of scenarios where we know what to expect
        ## Here we add modules and we add two additional epistasis scenarios
        ## (includig two genes, no epistasis). We also add order where it does
        ## not matter, and no interaction genes in the init mutants even when
        ## they have s = 0.
        reps <- 120
        mu <- 1e-7
        nig <- 20 ## there are lots of genes with modules with many genes
        out <- NULL
        outNO <- NULL
        so <- 0.2
        niG <- rep(0, nig)
        names(niG) <- paste0("nint", 1:nig)
        feo <- allFitnessEffects(orderEffects = c("A > B" = so,
                                                  "B > A" = -5),
                                 noIntGenes = niG,
                                 geneToModule = c(A = genmodule("a", 19),
                                                  B = genmodule("b", 15)))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a1 > b2 > nint9", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 190, "McFL", "a2 > nint6 > b8", 2.61e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a9 > b3", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 190, "McFL", "nint3 > a3 > b2", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 20, "Exp", "a4 > nint9 > nint2 > b1", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 190, "McFL", "nint3 > nint8 > a7 > b2", 1.5e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "Exp", "a2 > b13 > nint7", 1.67e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 199, "McFL", "a19 > b10 > nint6", 1.67e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b13 > a2", 5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 199, "McFL", "b1 > a2", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b14 > a8", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 199, "McFL", "b8 > a6", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b7 > a2", 5e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 199, "McFL", "b12 > a9", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 30, "Exp", "b15 > a4", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 170, "McFL", "b5 > a5", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 200, "McFL", "b2 > a1", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 300, "McFL", "b1 > a2", 1e2, reps, mu, verboseOandE))
        so <- 0.32
        niG <- c(so, rep(0, nig))
        names(niG) <- replicate(nig + 1, paste(sample(letters, 12, replace = TRUE),
                                               collapse = ""))
        names(niG)[1] <- "ThisisA"
        feo <- allFitnessEffects(noIntGenes = niG)
        ## add some fun 
        out <- rbind(out, OandE(feo, so, 30, "Exp", paste("ThisisA >", names(niG)[2]),
                                5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 130, "McFL", paste("ThisisA >", names(niG)[5]),
                                1.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 33, "Exp",
                                paste(names(niG)[3], "> ThisisA >", names(niG)[2]),
                                4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 133, "McFL",
                                paste(names(niG)[5], "> ThisisA >", names(niG)[5]),
                                1.6e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp",
                                    paste(names(niG)[2], ">", names(niG)[4]),
                                    5.32e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 140, "McFL",
                                    paste(names(niG)[2], ">", names(niG)[5]),
                                    1.33e3, reps, mu, verboseOandE))
        so <- 0.17
        niG <- rep(0, nig)
        names(niG) <- paste0("nint", 1:nig)
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
        out <- rbind(out, OandE(feo, so, 29, "Exp", "a2, c3", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 199, "McFL", "a5, c1", 1.2e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "Exp", "nint1, a2, c3", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 186, "McFL", "nint6, a5, c1", 1.4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 42, "Exp", "nint1 > a2 > c3", 3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 196, "McFL", "nint6 > a5 > c1", 1.8e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 43, "Exp", "a2 > nint2 > c3", 2e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 246, "McFL", "a5 > nint9 > c1", 1.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "c1, a7, c2, a3, a9, a8, nint6", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 340, "McFL", "a5, c6, a2, nint3, c5", 1.21e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b2, c6", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 240, "McFL", "b4, c2", 2e3, reps, mu, verboseOandE))
        so <- 0.29 ## 0.09
        niG <- rep(0, nig)
        names(niG) <- paste0("nint", 1:nig)
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, so),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 geneToModule = c(
                                     "Root" = "Root",
                                     A = genmodule("a", 8),
                                     B = genmodule("b", 9),
                                     C = genmodule("c", 4)),
                                 noIntGenes = niG)
        out <- rbind(out, OandE(feo, so, 45, "Exp", "nint6, a6, a8, c1, c2", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", "a1, c3, nint4, c3", 3.1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b2, b3, c3, nint3", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "McFL", "b1, b4, nint2, c4", 1e3, reps, mu, verboseOandE))
        ft <- 60
        sa02 <- 0.1
        fe02 <- allFitnessEffects(epistasis = c("A" = sa02),
                                  noIntGenes = niG,
                                  geneToModule = c(A = genmodule("a", 10)))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "a2, nint3", 4.4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 5 * ft, "McFL", "nint1, a7", 3.21e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "a1, a4, nint6", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 5 * ft, "McFL", "a6, nint4, a3", 5.22e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "a4, a1, a3", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 5 * ft, "McFL", "a2, a8", 1.17e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 30, "Exp", "a3, a6", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 6 * 30, "McFL", "a7, a5, a4, a2, a3", 4.3e4, reps, mu, verboseOandE))
        ft <- 15
        sa02 <- 0.1
        sb02 <- 0.2
        sc <- sa02 + sb02 + sa02 * sb02 ## we pass the single one of the double mutant, when no epist
        fe02 <- allFitnessEffects(epistasis = c("A" = sa02,
                                                "B" = sb02),
                                  noIntGenes = niG,
                                  geneToModule = c(A = genmodule("a", 10),
                                                   B = genmodule("b", 5)))
        out <- rbind(out, OandE(fe02, sc, ft, "Exp", "nint3, nint6, a2, a3, b2", 3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 10 * ft, "McFL", "a7, b5", 1.8e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, ft, "Exp", "a1 > b2 > b3 > a4", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 10 * ft, "McFL", "a6 > b2 > b4 > nint2 > a3", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, ft, "Exp", "a4, a1, a3, b4", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 10 * ft, "McFL", "a2 > a8 > b5 > b1 > b3", 1.8e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 17, "Exp", "a3 > a6 > nint5 >  b3", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 10 * 19, "McFL", "a7 >  a8 > a9 > b5 > b4 > a2 > a3", 1e4, reps, mu, verboseOandE))
        sd <- sa02 * sb02
        outNO <- rbind(outNO, OandE(fe02, sd, ft, "Exp", "a4, a1, a3, b4", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, ft, "McFL", "a2, a8, b5, b1, b3", 1.9e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, 30, "Exp", "a3, a6, b3", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, 30, "McFL", "a7, a8, a9, b5, b4, a2, a3", 1e4, reps, mu, verboseOandE))
        ft <- 60
        sa02 <- 0.3
        sb02 <- 0.6
        sc <- sa02 * sb02 ## we pass the single one of the double mutant, when no epist
        fe02 <- allFitnessEffects(epistasis = c("A : B" = sa02 * sb02),
                                  noIntGenes = niG,
                                  geneToModule = c(A = genmodule("a", 10),
                                                   B = genmodule("b", 5)))
        out <- rbind(out, OandE(fe02, sc, ft, "Exp", "a2 > a3 > b2 > nint5", 3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 9 * ft, "McFL", "a7 > nint3 > b5", 3.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 8 * ft, "McFL", "a7 > nint3 > b5", 5e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 17, "Exp", "a1, b2, b3, a4, nint2", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 10 * 49, "McFL", "a6, b2, b4, nint4, a3", 5.2e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 34, "Exp", "a4 > a1 > a3 > b4 > nint6", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 2 * 60, "McFL", "a4 > a1 > a3 > b4 > nint6", 4e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 9 * ft, "McFL", "a2 > a8 > nint3 > b5 > b1 > b3", 1.45e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 30, "Exp", "a3, a6, b3", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sc, 9 * 30, "McFL", "a7, a8, a9, b5, b4, a2, a3", 1.26e4, reps, mu, verboseOandE))
        sd <- sa02 + sb02 + sa02 * sb02
        outNO <- rbind(outNO, OandE(fe02, sd, ft, "Exp", "a4, a1, a3, b4, nint5", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, 8 * ft, "McFL", "a2 > nint1 > a8 > b5 > b1 > b3", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, 30, "Exp", "a3 > a6 > b3", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(fe02, sd, 90, "McFL", "a7, a8, a9, b5, b4, a2, a3", 1e4, reps, mu, verboseOandE))
        out <- data.frame(out)
        colnames(out) <- c("Expected", paste0("Observed_", 1:reps))
        d1 <- data.frame(Expected = out[, 1], Observed = rowMeans(out[, -1]))
        p.fail <- 0.05
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
        T.lm.diff <- all(abs(coefficients(lm1) - c(0, 1) ) < 0.1)
        if(! (T.not && T.yest && T.lm && T.lm.diff) ) {
            cat("\n T.not is \n"); print(T.not)
            print(no.t)
            cat("\n T.yest \n"); print(T.yest)
            print(yes.t)
            cat("\n larger values of T.yest II \n")
            print(which(p.adjust(yes.t, method = "BH") <= p.fail))
            cat("\n lm and plot\n")
            print(summary(lm1)); print(car::linearHypothesis(lm1, diag(2), c(0, 1)))
            plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        }
        plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        if ((T.not && T.yest && T.lm && T.lm.diff) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true((T.not && T.yest && T.lm && T.lm.diff) )
})
date()


date()
test_that("Observed vs expected, case III", {
    max.tries <- 5
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n Observed vs expected, case III\n")
        genmodule <- function(l, num) {
            paste(paste0(l, 1:num), collapse = ", ")
        }
        ## Like II, but we combine several effects.
        reps <- 120
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
        out <- rbind(out, OandE(feo, so1, 130, "McFL", "a2 > nint1 > b8", 1.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 15, "Exp", "a9 > b3 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 140, "McFL", "nint3 > nint1 > a3 > b2", 1.4e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 20, "Exp", "a4 > nint9 > nint1 > b1", 1e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 140, "McFL", "nint1 > nint8 > a7 > b2", 3.1e2, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 12, "Exp", "a2 > b13 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 250, "McFL", "a19 > b10 > nint6 > nint1", 3.1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 12, "Exp", "a2 > b13", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 250, "McFL", "a19 > b10 > nint6", 1.9e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 12, "Exp", "a2 > nint4 > b13", 5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so1, 250, "McFL", "a19 > b10 > nint6 > nint8", 7.1e4, reps, mu, verboseOandE))
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
        out <- rbind(out, OandE(feo, so1, 140, "McFL", "a5, nint1, c1", 1.6e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 10, "Exp", "nint1, a2, c3", 1.7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 236, "McFL", "nint6, nint1, a5, c1", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 32, "Exp", "nint1 > a2 > c3", 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 246, "McFL", "nint6 > a5 > c1 > nint1", 4.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 18, "Exp", "a2 > nint2 > nint1 > c3", 3.1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 146, "McFL", "a5 > nint9 > c1 > nint1", 1.97e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 20, "Exp", "c1, a7, c2, a3, a9, a8, nint1, nint7", 1.76e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so1, 140, "McFL", "a5, nint1, c6, a2, nint3, c5", 3.3e3, reps, mu, verboseOandE))
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
        out <- rbind(out, OandE(feo, st, 145, "McFL",
                                "a2 > f3 > e2 > c1 > g4 > h2 > nint1", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 135, "McFL",
                                "a2 > a4 > f3 > e2 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4.7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 235, "McFL",
                                "a6 > f1 > e3 > g4 > c1 > nint1 > h2 > nint3", 5.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 230, "McFL",
                                "f3 > e2 > c1 > g4 > h2 > nint1 > a4", 4.21e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 240, "McFL",
                                "h2 > f3 > e2 > c1 > g4 > a3 > nint1 > nint4", 5.1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, st, 238, "McFL",
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
        outNO <- rbind(outNO, OandE(feo, st, 145, "McFL",
                                    "a2 > e2 > c1 > g4 > h2 > nint1", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 235, "McFL",
                                    "a2 > a4 > e3 > f1 > e3 > nint6 > c4 > c1 > g4 > g1 > h2 > nint1", 4e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 235, "McFL",
                                    "a6 > f1 > e3 > c1 > nint1 > h2 > nint3", 5e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 230, "McFL",
                                    "f3 > e2 > c1 > g4 > nint1 > a4", 4e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 240, "McFL",
                                    "h2 > e2 > c1 > g4 > a3 > nint1 > nint4", 5e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, st, 238, "McFL",
                                    "a7 > f3 > c1 > g4 > a6 > h2 > nint1 > nint3", 1e3, reps, mu, verboseOandE))
        out <- data.frame(out)
        colnames(out) <- c("Expected", paste0("Observed_", 1:reps))
        d1 <- data.frame(Expected = out[, 1], Observed = rowMeans(out[, -1]))
        p.fail <- 0.05
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
        T.lm.diff <- all(abs(coefficients(lm1) - c(0, 1) ) < 0.1)
        if(! (T.not && T.yest && T.lm && T.lm.diff) ) {
            cat("\n T.not is \n"); print(T.not)
            print(no.t)
            cat("\n T.yest \n"); print(T.yest)
            print(yes.t)
            cat("\n larger values of T.yest III \n")
            print(which(p.adjust(yes.t, method = "BH") <= p.fail))
            cat("\n lm and plot\n")
            print(summary(lm1)); print(car::linearHypothesis(lm1, diag(2), c(0, 1)))
            plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        }
        if ((T.not && T.yest && T.lm && T.lm.diff) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true((T.not && T.yest && T.lm && T.lm.diff) )

})
date()



date()
test_that("Init mutant no effects if fitness is 0", {
    max.tries <- 5
    for(tries in 1:max.tries) {
        T1 <- T2 <- T3 <- T4 <- T5 <- T6 <- T7 <- T8 <- TRUE
        cat("\n Init mutant no effect if fitness is 0\n")
        ## This can give false positives often.
        ## s should never be relevant concern, so increase reps but
        ## keep fast by decreasing sampling
        sEvery <- 1
        reps <- 200 
        mu <- 1e-7
        nig <- 50
        out <- NULL
        outNO <- NULL
        so <- 0
        feo <- allFitnessEffects(orderEffects = c("a > b" = so,
                                                  "b > a" = -5),
                                 noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(feo, so, 20, "Exp", "a > b", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", "a > b", 4e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a > b", 1.27e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", "a > b", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", "a > b", 1e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", "a > b", 1e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "Exp", "a > b", 1.3e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 130, "McFL", "a > b", 1.3e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 20, "Exp", NULL, 7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", NULL, 7e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", NULL, 1.5e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", NULL, 1.5e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", NULL, 1.6e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 140, "McFL", NULL, 1.6e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "Exp", NULL, 1.8e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 30, "McFL", NULL, 1.8e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 140, "McFL", "b > a", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 140, "McFL", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "b > a", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 30, "Exp", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 30, "McFL", "b > a", 1e4, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 200, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 300, "McFL", "b > a", 1e2, reps, mu, verboseOandE))
        so <- 0
        niG <- c(so, rep(0, nig))
        names(niG) <- replicate(nig + 1, paste(sample(letters, 12, replace = TRUE),
                                               collapse = ""))
        names(niG)[1] <- "ThisisA"
        feo <- allFitnessEffects(noIntGenes = niG)
        out <- rbind(out, OandE(feo, so, 40, "Exp", "ThisisA", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", "ThisisA", 4.6e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", NULL, 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", NULL, 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "Exp", names(nig)[2], 4e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", names(nig)[2], 4e4, reps, mu, verboseOandE))
        so <- 0
        feo <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                            child  = c("A", "B", "C"),
                                            s = c(0, -1, so),
                                            sh = rep(-1, 3),
                                            typeDep = "MN"
                                            ),
                                 noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(feo, so, 20, "Exp", "A, C", 3e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", "A, C", 3e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 10, "Exp", NULL, 3e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(feo, so, 40, "McFL", NULL, 3e5, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "Exp", "B, C", 1e3, reps, mu, verboseOandE))
        outNO <- rbind(outNO, OandE(feo, so, 40, "McFL", "B, C", 1e3, reps, mu, verboseOandE))
        ft <- 10
        sa02 <- 0
        fe02 <- allFitnessEffects(epistasis = c("A" = sa02),
                                  noIntGenes = rep(0, nig))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 2e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", "A", 2e4, reps, mu, verboseOandE)) ## 28
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", "A", 5e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", "A", 1e5, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", "A", 4.8e4, reps, mu, verboseOandE)) ## 32
        out <- rbind(out, OandE(fe02, sa02, 30, "Exp", "A", 1e4, reps, mu, verboseOandE)) ## 33
        out <- rbind(out, OandE(fe02, sa02, 30, "McFL", "A", 1e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", NULL, 2e3, reps, mu, verboseOandE))  
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", NULL, 1e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", NULL, 5.3e3, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", NULL, 5.3e3, reps, mu, verboseOandE)) ## 38
        out <- rbind(out, OandE(fe02, sa02, ft, "Exp", NULL, 5.7e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, ft, "McFL", NULL, 5.7e4, reps, mu, verboseOandE)) ## 40
        out <- rbind(out, OandE(fe02, sa02, 30, "Exp", NULL, 1.4e4, reps, mu, verboseOandE))
        out <- rbind(out, OandE(fe02, sa02, 30, "McFL", NULL, 1.4e4, reps, mu, verboseOandE))
        out <- data.frame(out)
        colnames(out) <- c("Expected", paste0("Observed_", 1:reps))
        d1 <- data.frame(Expected = out[, 1], Observed = rowMeans(out[, -1]))
        p.fail <- 0.05
        lm1 <- lm(log(Observed) ~ log(Expected), data = d1)
        ## summary(lm1)
        ## plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        ## abline(lm1); abline(a = 0, b = 1, col = "red")
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
        T.lm.diff <- all(abs(coefficients(lm1) - c(0, 1) ) < 0.1)
        if(! (T.not && T.yest && T.lm && T.lm.diff) ) {
            cat("\n T.not is \n"); print(T.not)
            print(no.t)
            cat("\n T.yest \n"); print(T.yest)
            print(yes.t)
            ytbug <<- yes.t
            outbug <<- out
            cat("\n larger values of T.yest 0 \n")
            print(which(p.adjust(yes.t, method = "BH") <= p.fail))
            cat("\n lm and plot\n")
            print(summary(lm1)); print(car::linearHypothesis(lm1, diag(2), c(0, 1)))
            plot(log(Observed) ~ log(Expected), data = d1); abline(lm1); abline(a = 0, b = 1, col = "red")
        }
        if ((T.not && T.yest && T.lm && T.lm.diff) ) break;
    }
    cat(paste("\n done tries", tries, "\n"))
    expect_true((T.not && T.yest && T.lm && T.lm.diff) )

})
date()

cat(paste("\n            a final runif ", runif(1), "\n"))
cat(paste("\n Ending fitness preds long at", date()))


## Leave III in regular tests, move rest to long
