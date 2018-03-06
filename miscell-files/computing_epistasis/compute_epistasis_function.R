rm(list=ls(all=TRUE))
library(gtools)
library(OncoSimulR)
library(compiler)
options(digits=8)
enableJIT(3)
p <- 0.0000001   ## p is the accuracy of the comparation in CalcEpis function

## RDU: FIXME
##      - do not use gtools anymore?
##      - the p: use all.equal?
##      - tests
##         - use files  example-magellan-bug-sign-epistasis.R
##                      magellan_rsign_epistasis_examples.R

genes <- 9
frfitness <- rfitness(genes)

## while(TRUE){
##     iter <- iter+1
##     cat("\n Doing iter ", iter, "\n")

                                        #for (genes in 3:3){   
## for (qwe in 1:100){
## gnames <- LETTERS[1:5]
## (loci_bg <- combinations(length(gnames), length(gnames) - 2, v = gnames))

## names of background loci -> matrix with all genotypes for background
matrix_genotypes_background <- function(z) {
    y <- OncoSimulR:::generate_matrix_genotypes(length(z))
    colnames(y) <- z
    return(y)
}

## matrix of all 4 genotype combinations of the epistatic loci
the_epist_loci <- function(genes) {
    tmp <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
    colnames(tmp) <- genes
    return(tmp)
}

## From
## https://stat.ethz.ch/pipermail/r-help/2007-February/125739.html
bin_to_int <- function(x) sum(2^(which(rev(as.logical(x))) - 1))



## vector of background, genenames, fitness landscape, index into fitness landscape -> vector of fitness of all four genotypes (and its epistasis)
epist_for_background <- function(x, gnames, fitnessl, intind) {
    u2 <- matrix(rep(x, 4), nrow = 4, byrow = TRUE)
    colnames(u2) <- names(x)
    u3 <- cbind(u2, the_epist_loci(sort(setdiff(gnames, names(x)))))
    ogn <- order(colnames(u3))
    u3 <- u3[, ogn]

    binintu3 <- apply(u3, 1, bin_to_int)
    fpos <- match(binintu3, intind)

    the_four_fitness <- fitnessl[fpos, "Fitness"]
    ## do calculations epistasis and return them
    return(the_four_fitness)
}




compute_epistasis <- function(fitnessl) {
    ## Make sure fitness landscape has columns ordered by gene name
    ocn <- order(colnames(fitnessl)[-ncol(fitnessl)])
    r2 <- cbind(fitnessl[, ocn], Fitness = fitnessl[, "Fitness"])
    gnames <- colnames(r2)[-ncol(r2)]

    (loci_bg <- combinations(length(gnames), length(gnames) - 2, v = gnames))

    tg <- length(gnames)
    output <- matrix(nrow = (tg * (tg -1) * 0.5 * (2^(tg-2))),
                     ncol = 4 )
    k <- 0

    
    squeleton <- OncoSimulR:::generate_matrix_genotypes(tg-2)

    ## create the index for genotypes in fitnesslandscape
    intind <- apply(fitnessl[, -ncol(fitnessl)], 1,
                    bin_to_int)
    stopifnot(!any(duplicated(intind)))
    
    for(i in 1:nrow(loci_bg)) {
        u <- squeleton
        colnames(u) <- loci_bg[i, ]
        for(j in 1:nrow(u)) {
            k <- k + 1
            output[k, ] <- (epist_for_background(u[j, ], gnames, r2, intind))
        }
    }
    return(output)
}

CalcEpis <- function(calc){
    epistasia_de_signo <- 0
    epistasia_de_signo_reciproco <- 0
    magnitud_epistasis <- 0

    A <- calc[1]
    B <- calc[2]
    C <- calc[3]
    D <- calc[4]

    Q1 <- abs(B-A + D-C)
    Q2 <- abs(B-A) + abs(D-C)
    Q3 <- abs(C-A + D-B)
    Q4 <- abs(C-A) + abs(D-B)
    Q5 <- B - A
    Q6 <- D - C
    
    if (Q5 != Q6 && Q1 + p > Q2 && Q1 - p < Q2 && Q3 + p > Q4 && Q3 - p < Q4){magnitud_epistasis <- magnitud_epistasis + 1}
    
    if (Q1 < Q2 - p && Q3 + p > Q4 || Q3 < Q4 - p && Q1 + p > Q2){epistasia_de_signo <- epistasia_de_signo + 1}
    
    if (Q1 < Q2 - p && Q3 < Q4 - p){
        epistasia_de_signo_reciproco <- epistasia_de_signo_reciproco + 1}

    epis <- c(magnitud_epistasis,epistasia_de_signo,epistasia_de_signo_reciproco)
    return(epis)}

## r1 <- rfitness(5)

ComputeEpistasis <- function(frfitness){
    start.time <- Sys.time()    
    calc <- compute_epistasis(frfitness)
    Pepistasis <- apply((apply(calc, 1, CalcEpis)), 1, sum) / nrow(calc)
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    
    start.time2 <- Sys.time()
    epis_magellan <- OncoSimulR:::Magellan_stats(frfitness, use_log=FALSE, verbose=TRUE) 
    q <- c(epis_magellan[[11]], epis_magellan[[12]], epis_magellan[[13]])    
    end.time2 <- Sys.time()
    time.taken2 <- end.time2 - start.time2

    output <- rbind(c(Pepistasis,time.taken),c(q,time.taken2))
    colnames(output) <- c("Magnitud","Sign","Reciprocal Sign","Time taken")
    rownames(output) <- c("Our code", "Magellan code")
    return(output)
}    


Types_of_epistasis <- ComputeEpistasis(frfitness)











