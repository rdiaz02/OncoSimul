rm(list=setdiff(ls(), c("timev","ng","r2", "time1", "contime", "timecomp")))
library(gtools)
library(OncoSimulR)
library(compiler)
options(digits=8)
enableJIT(3)
library(testthat)
library(Rcpp)
option=1

## FIXME:
## allow to work if "sparse fitness landscape" (no entries with fitness <= 0)?
## Make it faster? recode in C++ the inner function?


## RDU: don't do this. Use machine accuracy
## p <- 0.0000001   ## p is the accuracy of the comparation in CalcEpis function

## RDU: FIXME
##      - do not use gtools anymore?
##      - the p: use a tolerance (and I use the one from all.equal)
##      - tests
##         - use files  example-magellan-bug-sign-epistasis.R
##                      magellan_rsign_epistasis_examples.R

## genes <- 9
## frfitness <- rfitness(genes)

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
    tmp <- matrix(c(0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L),
                  ncol = 2, byrow = TRUE)
    colnames(tmp) <- genes
    return(tmp)
}

## From
## https://stat.ethz.ch/pipermail/r-help/2007-February/125739.html
bin_to_int <- function(x) sum(2^(which(rev(as.logical(x))) - 1))




## RDU: but this does not compute epistasis. Misleading name.
## vector of background, genenames, fitness landscape, index into fitness landscape -> vector of fitness of all four genotypes (and its epistasis)
## epist_for_background <- function(x, gnames, fitnessl, intind) {
four_fitn_for_bg <- function(x, gnames, fitnessl, intind) {
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


########################################################################
########################################################################
####
####               These functions are no longer used directly
####              (modified versions of them are used)
####              
########################################################################
########################################################################


## this does not compute epistasis either. Misleading name
## compute_epistasis <- function(fitnessl) {
fitness_of_four_genotypes <- function(fitnessl) {
    browser()
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
            ## output[k, ] <- (epist_for_background(u[j, ], gnames, r2, intind))
            output[k, ] <- (four_fitn_for_bg(u[j, ], gnames, r2, intind))
        }
    }
    return(output)
}


## CalcEpis <- function(calc,
##                      tolerance = sqrt(.Machine$double.eps)){
##     sign_epist <- 0
##     r_sign_epist <- 0
##     magn_epist <- 0

##     A <- calc[1]
##     B <- calc[2]
##     C <- calc[3]
##     D <- calc[4]

##     Q1 <- abs(B-A + D-C)
##     Q2 <- abs(B-A) + abs(D-C)
##     Q3 <- abs(C-A + D-B)
##     Q4 <- abs(C-A) + abs(D-B)
##     Q5 <- B - A
##     Q6 <- D - C
    
##     if ( (Q5 != Q6) && (Q1 + tolerance > Q2) &&
##          (Q1 - tolerance < Q2 ) &&
##          (Q3 + tolerance > Q4 ) &&
##          (Q3 - tolerance < Q4) ){
##         magn_epist <- magn_epist + 1
##     }

##     ## RDU: what is the order of evaluation here if you do not use parentheses?
##     ## never relay on possibly obscure syntactic rules about operator precedence
##     if ( ((Q1 < Q2 - tolerance) &&
##           (Q3 + tolerance > Q4) ) ||
##          ((Q3 < Q4 - tolerance) &&
##           (Q1 + tolerance > Q2) )){
##         sign_epist <- sign_epist + 1
##     }
    
##     if ( (Q1 < Q2 - tolerance) && (Q3 < Q4 - tolerance)) {
##         r_sign_epist <- r_sign_epist + 1
##     }
    
##     epis <- c(magnitude = magn_epist,
##               sign = sign_epist,
##               r_sign = r_sign_epist)
##     return(epis)
## }

## r1 <- rfitness(5)

## RDU: this mixes computation of epistasis with comparison with Magellan
## ComputeEpistasis <- function(frfitness){ 
##     ## RDU: do not time this way. Use system.time
##     start.time <- Sys.time()    
##     ## calc <- compute_epistasis(frfitness)
##     calc <- fitness_of_four_genotypes(frfitness)
##     Pepistasis <- apply((apply(calc, 1, CalcEpis)), 1, sum) / nrow(calc)
##     end.time <- Sys.time()
##     time.taken <- end.time - start.time
    
##     start.time2 <- Sys.time()
##     epis_magellan <- OncoSimulR:::Magellan_stats(frfitness, use_log=FALSE, verbose=TRUE) 
##     q <- c(epis_magellan[[11]], epis_magellan[[12]], epis_magellan[[13]])    
##     end.time2 <- Sys.time()
##     time.taken2 <- end.time2 - start.time2

##     output <- rbind(c(Pepistasis,time.taken),c(q,time.taken2))
##     colnames(output) <- c("Magnitud","Sign","Reciprocal Sign","Time taken")
##     rownames(output) <- c("Our code", "Magellan code")
##     return(output)
## }    


## Types_of_epistasis <- ComputeEpistasis(frfitness)


######################################################################
######################################################################
######################################################################
######################################################################

## Maybe rewrite in C++?

## Based on CalcEpis, with some checks and no addition
## but direct assignment
## vector of fitness -> epistasis (0/1) for sign, rsign, magnitude
##     vector: background, single mutant in one locus,
##                         single mutant in the other, double mutant
##     The _q_ is for quartet
if (option == 1){
    single_q_epis <- function(calc,
                              tolerance = sqrt(.Machine$double.eps)){
        stopifnot(is.vector(calc))
        sign_epist <- 0L
        r_sign_epist <- 0L
        magn_epist <- 0L

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
        
        if ( (Q5 != Q6) && (Q1 + tolerance > Q2) &&
             (Q1 - tolerance < Q2 ) &&
             (Q3 + tolerance > Q4 ) &&
             (Q3 - tolerance < Q4) ){
            ## magn_epist <- magn_epist + 1
            magn_epist <- 1L
        }

        ## RDU: what is the order of evaluation here if you do not use parentheses?
        ## never relay on possibly obscure syntactic rules about operator precedence
        if ( ((Q1 < Q2 - tolerance) &&
              (Q3 + tolerance > Q4) ) ||
             ((Q3 < Q4 - tolerance) &&
              (Q1 + tolerance > Q2) )){
            ## sign_epist <- sign_epist + 1
            sign_epist <- 1L
        }
        
        if ( (Q1 < Q2 - tolerance) && (Q3 < Q4 - tolerance)) {
            ## r_sign_epist <- r_sign_epist + 1
            r_sign_epist <- 1L
        }
        
        epis <- c(magnitude = magn_epist,
                  sign = sign_epist,
                  r_sign = r_sign_epist)
        return(epis)
    }
}
if (option == 2){
   cppFunction('std::vector<double> single_q_epis(double x1, double x2, double x3, double x4, double tolerance){
#include <vector>
#include <iostream>
using namespace std;
double sign_epist = 0;
double r_sign_epist = 0;
double magn_epist = 0;
double Q1 = abs(x2-x1 + x4-x3);
double Q2 = abs(x2-x1) + abs(x4-x3);
double Q3 = abs(x3-x1 + x4-x2);
double Q4 = abs(x3-x1) + abs(x4-x2);
double Q5 = x2 - x1;
double Q6 = x4 - x3;
    if ((Q5 != Q6) &&
(Q1 + tolerance > Q2) &&
 (Q1 - tolerance < Q2 ) &&
 (Q3 + tolerance > Q4) &&
 (Q3 - tolerance < Q4) )
{magn_epist = magn_epist + 1;}

    if (((Q1 < Q2 - tolerance) &&
          (Q3 + tolerance > Q4) ) ||
         ((Q3 < Q4 - tolerance) &&
          (Q1 + tolerance > Q2) )){
        sign_epist = sign_epist + 1;
    }
    
    if ( (Q1 < Q2 - tolerance) && (Q3 < Q4 - tolerance)) {
        r_sign_epist = r_sign_epist + 1;
    }

std::vector<double> epis{magn_epist, sign_epist, r_sign_epist};
return epis;
}')
}



## Based on compute_epistasis

## fitness landscape, tolerance -> epistasis (number and fraction)
epistasis <- function(fitnessl, tolerance = sqrt(.Machine$double.eps)) {
    ## Make sure fitness landscape has columns ordered by gene name
    ocn <- order(colnames(fitnessl)[-ncol(fitnessl)])
    r2 <- cbind(fitnessl[, ocn], Fitness = fitnessl[, "Fitness"])
    gnames <- colnames(r2)[-ncol(r2)]

    (loci_bg <- combinations(length(gnames), length(gnames) - 2, v = gnames))

    tg <- length(gnames)

    ## Decrease memory usage by returning the epistasis directly. Note,
    ## though, that if we were to want other measures, keeping the quartet
    ## output might be better. Uncomment and comment as needed.
    
    ## If we were to return the quartet
    ## output <- matrix(nrow = (tg * (tg -1) * 0.5 * (2^(tg-2))),
    ##                  ncol = 4 )
    ## Return the 3-vector epistasis
    output <- matrix(nrow = (tg * (tg -1) * 0.5 * (2^(tg-2))),
                     ncol = 3)

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
            ## output[k, ] <- four_fitn_for_bg(u[j, ], gnames, r2, intind)
if (option == 1){
            output[k, ] <- single_q_epis(four_fitn_for_bg(u[j, ], gnames,
                                                          r2, intind),
                                         tolerance)}
            if (option == 2){
                calc <- four_fitn_for_bg(u[j, ], gnames,r2, intind)
                #print(calc)
            output[k, ] <- single_q_epis(calc[1],calc[2],calc[3],calc[4],sqrt(.Machine$double.eps)) 
                }
        }
    }
    ## epi <- t(apply(output, 1, single_q_epis))
    ## epis <- colSums(epi)
    ## Names for compatibility with MAGELLAN
    colnames(output) <- c("epist_magn", "epist_sign", "epist_rsign")
    epis <- colSums(output)
    epis <- c(epis, epis/nrow(output))
    names(epis) <- c(paste0("num_", colnames(output)),
                     colnames(output))
    return(epis)
}





## ## an alternative

## ## vector of background, genenames, fitness landscape, index into fitness landscape -> epistasis
## epist_for_bg <- function(x, gnames, fitnessl, intind) {
##     u2 <- matrix(rep(x, 4), nrow = 4, byrow = TRUE)
##     colnames(u2) <- names(x)
##     u3 <- cbind(u2, the_epist_loci(sort(setdiff(gnames, names(x)))))
##     ogn <- order(colnames(u3))
##     u3 <- u3[, ogn]

##     binintu3 <- apply(u3, 1, bin_to_int)
##     fpos <- match(binintu3, intind)

##     the_four_fitness <- fitnessl[fpos, "Fitness"]
##     ## do calculations epistasis and return them
##     return(CalcEpis(the_four_fitness))
## }


