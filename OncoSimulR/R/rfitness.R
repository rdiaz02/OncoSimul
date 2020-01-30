## Copyright 2013, 2014, 2015, 2016, 2017 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


create_eq_ref <- function(g) {
    ## "random" gives more prob. to genotypes with
    ## number of mutated genes close to g/2.
    ## This gives equal prob to having the reference
    ## be of any of the possible number of mutated genes.
    nm <- sample(g, 1)
    ref <- c(rep(1, nm), rep(0, g - nm))
    sample(ref)
}



get_magellan_binaries <-
    function(names_binaries = c("fl_statistics", "fl_generate", "fl_genchains"))
{
    rarch <- Sys.getenv('R_ARCH')
    nn_names_binaries <- names_binaries
    if(.Platform$OS.type == "windows")
        names_binaries <- paste0(names_binaries, ".exe")
    if(nzchar(rarch)) {
        rarch <- sub("^/", "", rarch)
        magellan_binaries <-  system.file(package = "OncoSimulR", "exec",
                                          rarch, names_binaries)
    } else {
        magellan_binaries <-  system.file(package = "OncoSimulR", "exec",
                                          names_binaries)
    }
    names(magellan_binaries) <- nn_names_binaries
    return(magellan_binaries)
}

## pkg.env <- new.env()
## ## The next will not install with error
## ## ERROR: hard-coded installation path:
## please report to the package maintainer and use ‘--no-staged-install’
## pkg.env <- c(pkg.env, get_magellan_binaries())



fl_statistics_binary <- function() get_magellan_binaries("fl_statistics")
fl_generate_binary <- function() get_magellan_binaries("fl_generate")


rfitness <- function(g, c= 0.5,
                     sd = 1,
                     mu = 1,
                     reference = "random", ## "random", "max", or the vector,
                                     ## e.g., rep(1, g). If random, a
                                     ## random genotype is chosen as
                     ## reference. If "max" this is rep(1, g)
                     scale = NULL, ## a two-element vector: min and max
                     wt_is_1 = c("subtract", "divide", "force", "no"),
                     ## wt_is_1 = TRUE, ## wt has fitness 1
                     log = FALSE, ## log only makes sense if all values >
                                 ## 0. scale with min > 0, and/or set
                                 ## wt_is_1 = divide
                     min_accessible_genotypes = NULL,
                     accessible_th = 0,
                     truncate_at_0 = TRUE,
                     K = 1,
                     r = TRUE,
                     i = 0, # Ising, cost incompatibility
                     I = -1, # Ising, sd for "i"
                     circular = FALSE, # Ising, circular arrangement
                     e = 0, # Eggbox, +/- e
                     E = -1, # Eggbox, noise on "e"
                     H = -1, # HoC stdev
                     s = 0.1, # mean multiplivative
                     S = -1, # SD multiplicative
                     d = 0, # disminishing/increasing for multiplicative
                     o = 0, # mean optimum
                     O = -1, # sd optimum
                     p = 0, # mean production for non 0 allele (optimum)
                     P = -1, # sd for p
                    
                     model = c("RMF", "Additive", 
                               "NK", "Ising", "Eggbox", "Full")) {
    ## Like Franke et al., 2011 and others of Krug. Very similar to Greene
    ## and Crona, 2014. And this allows moving from HoC to purely additive
    ## changing c and sd.

    ## FIXME future: do this with order too?
    ##    - do not generate a matrix of genotypes but a matrix of ordered genot.
    wt_is_1 = match.arg(wt_is_1)
    model = match.arg(model)
    if(is_null_na(g)) stop("Number of genes argument (g) is null or NA")
    m <- generate_matrix_genotypes(g)
    done <- FALSE
    ## attempts <- 0 ## for debugging/tracking purposes
    while(!done) {
        if(model == "RMF") {
            ## attempts <- attempts + 1
            f_r <- rnorm(nrow(m), mean = mu, sd = sd)
            if(inherits(reference, "character") && length(reference) == 1) {
                if(reference == "random") {
                    referenceI <- m[sample(nrow(m), 1), ]
                } else if(reference == "max") {
                    referenceI <- rep(1, g)
                } else if(reference == "random2") {
                    referenceI <- create_eq_ref(g)
                }
            } else {
                referenceI <- reference
            }
            d_reference <- apply(m, 1, function(x) sum(abs(x - referenceI)))
            f_det <- -c * d_reference
            ## f_det <- rowSums(m) * slope/nrow(m) ## this is Greene and Krona
            fi <- f_r + f_det
        } else if (model == "Additive") {
      ## get fitness effect for mutations in each gene
      mutants <-rep(1,g)
            ## FIXME: Why not just?
            ## f_single_mut <- rnorm(g, mean = mu, sd = sd)
            f_single_mut <- sapply(mutants, FUN = function(x) 
                                            rnorm(x, mean = mu, sd = sd))
      ## find which gene is mutated 
      m2 <- m == 1
      ## Sum the fitness effect of that mutation to generate a vector fi with
      ## the fitness for each mutant condition
      fi <- apply(m2, MARGIN = 1, FUN = function (x) sum(x*f_single_mut))
      ## remove unnecessary variables
      rm (f_single_mut, m2)
        } else if(model == "NK") {
            if(K >= g) stop("It makes no sense to have K >= g")
            argsnk <- paste0("-K ", K,
                             ifelse(r, " -r ", " "),
                             g, " 2")
            fl1 <- system2(fl_generate_binary(), args = argsnk, stdout = TRUE)[-1]
        } else if (model == "Ising") {
      argsIsing <- paste0("-i ", i, " -I ", I ,
                          ifelse(circular, " -c ", " "),
                          g, " 2")
      fl1 <- system2(fl_generate_binary(), args = argsIsing, stdout = TRUE)[-1]
    } else if (model == "Eggbox") {
      argsEgg <- paste0("-e ", e, " -E ", E," ", g, " 2")
      fl1 <- system2(fl_generate_binary(), args = argsEgg, stdout = TRUE)[-1]
    } else if (model == "Full") {
      if(K >= g) stop("It makes no sense to have K >= g")
      argsFull <- paste0("-K ", K, ifelse(r, " -r ", " "),
                         "-i ", i, " -I ", I , ifelse(circular, " -c ", " "),
                         "-e ", e, " -E ", E, " ",
                         "-H ", H, " ",
                         "-s ", s, " -S ", S, " -d ", d, " ",
                         "-o ", o, " -O ", O, " -p ", p, " -P ", P, " ",
                         g, " 2")
      fl1 <- system2(fl_generate_binary(), args = argsFull, stdout = TRUE)[-1]
    }
    if (model == "Eggbox" || model == "Ising" || model == "Full" || model == "NK") {
            fl1 <- matrix(
                as.numeric(unlist(strsplit(paste(fl1, collapse = " "), " "))),
                ncol = g + 1, byrow = TRUE)
            m1 <- fl1[, 1:g]
            fi <- fl1[, g + 1]

            ## For scaling, etc, all that matters, if anything, is the wildtype

            ## We could order by doing this
            ## But I am not 100% sure this will always be the same as
            ## generate_matrix_genotypes
            ## oo <- do.call(order,
            ##               c(list(muts),
            ##                 as.list(data.frame(m1[, rev(1:ncol(m1))]))))
            ## m2 <- m1[oo, ]
            ## Or create an id and order by it?

            ## When we get to 20 genes, this is slow, about 18 seconds each
            ## apply. Matching is fast (< 0.5 seconds)
            gtstring <- apply(m, 1, function(x) paste0(x, collapse = ""))
            gtstring2 <- apply(m1, 1, function(x) paste0(x, collapse = ""))
            oo <- match(gtstring, gtstring2)
            fi <- fi[oo]
            ## make sure no left overs
            rm(gtstring, gtstring2, oo, fl1, m1)

            ## Had we not ordered, do this!!!
            ## Which one is WT?
            ## muts <- rowSums(m1)
            ## w_wt <- which(muts == 0)
            ## if(w_wt != 1) {
            ##     f_a <- fi[1]
            ##     fi[1] <- fi[w_wt]
            ##     fi[w_wt] <- f_a
            ##     rm(f_a)
            ## }
            ## m[] <- m1
            ## rm(m1)
            ## rm(fl1)
        }

        if(!is.null(scale)) {
            fi <- (fi - min(fi))/(max(fi) - min(fi))
            fi <- scale[1] + fi * (scale[2] - scale[1])
        }
        if(wt_is_1 == "divide") {
            ## We need to shift to avoid ratios involving negative numbers and
            ## we need to avoid having any fitness at 0, specially the wt.  If
            ## any negative value, add the min, and shift by the magnitude of
            ## the min to avoid any at 0.

            ## If you use scale and wt_is_1, this will move the scale. It is
            ## not possible to obtain a linear transformation that will keep
            ## the min and max of the scale, and wt at 1.
            min_fi <- min(fi)
            if(min_fi < 0)
                fi <- fi + 2 * abs(min(fi))
            fi <- fi/fi[1]
        } else if (wt_is_1 == "subtract") {
            fi <- fi - fi[1] + 1.0
        } else if(wt_is_1 == "force") {
            fi[1] <- 1.0
            if(!is.null(scale)) {
                if( (1 > scale[2]) || (1 < scale[1]))
                    warning("Using wt_is_1 = force and scale, but scale does ",
                            "not include 1")
            }
        }
        if(log) {
            ## If you want logs, you certainly do not want
            ## the log of a negative number
            fi[fi < 0] <- 0
            if(wt_is_1 == "no") {
                fi <- log(fi)
            } else {
                ## by decree, fitness of wt is 1. So shift everything
                fi <- log(fi) + 1
            }
            ## former expression, but it was more confusing
            ## fi <- log(fi/fi[1]) + 1
            
        }
        if(truncate_at_0) {
            ## yes, truncate but add noise to prevent identical
            fi[fi <= 0] <- runif(sum(fi <= 0),
                                 min = 1e-10,
                                 max = 1e-9)
        }
        m <- cbind(m, Fitness = fi)
        if(!is_null_na(min_accessible_genotypes)) {
            ## num_accessible_genotypes <- count_accessible_g(m, accessible_th)
            ## Recall accessibleGenotypes includes the wt, if accessible.
            num_accessible_genotypes <- length(wrap_accessibleGenotypes(m, accessible_th)) - 1
            ## message("\n     num_accessible_genotypes = ", num_accessible_genotypes, "\n")
            if(num_accessible_genotypes >= min_accessible_genotypes) {
                done <- TRUE
                attributes(m) <- c(attributes(m),
                                   accessible_genotypes = num_accessible_genotypes,
                                   accessible_th = accessible_th)
            } else {
                ## Cannot start again with a fitness column
                m <- m[, -ncol(m), drop = FALSE]
            }
        } else {
            done <- TRUE
        }
    }
    ## message("\n number of attempts = ", attempts, "\n")
    class(m) <- c(class(m), "genotype_fitness_matrix")
    return(m)
}










## rfitness <- function(g, c= 0.5,
##                      sd = 1,
##                      mu = 1,
##                      reference = "random", ## "random", "max", or the vector,
##                                      ## e.g., rep(1, g). If random, a
##                                      ## random genotype is chosen as
##                                      ## reference. If "max" this is rep(1, g)
##                      scale = NULL, ## a two-element vector: min and max
##                      wt_is_1 = c("subtract", "divide", "force", "no"),
##                      ## wt_is_1 = TRUE, ## wt has fitness 1
##                      log = FALSE, ## log only makes sense if all values >
##                                  ## 0. scale with min > 0, and/or set
##                                  ## wt_is_1 = divide
##                      min_accessible_genotypes = NULL,
##                      accessible_th = 0,
##                      truncate_at_0 = TRUE) {
##     ## Like Franke et al., 2011 and others of Krug. Very similar to Greene
##     ## and Crona, 2014. And this allows moving from HoC to purely additive
##     ## changing c and sd.

##     ## FIXME future: do this with order too?
##     ##    - do not generate a matrix of genotypes but a matrix of ordered genot.
##     wt_is_1 = match.arg(wt_is_1)
##     if(is_null_na(g)) stop("Number of genes argument (g) is null or NA")
##     m <- generate_matrix_genotypes(g)
##     done <- FALSE
##     ## attempts <- 0 ## for debugging/tracking purposes
##     while(!done) {
##         ## attempts <- attempts + 1
##         f_r <- rnorm(nrow(m), mean = mu, sd = sd)
##         if(inherits(reference, "character") && length(reference) == 1) {
##             if(reference == "random") {
##                 referenceI <- m[sample(nrow(m), 1), ]
##             } else if(reference == "max") {
##                 referenceI <- rep(1, g)
##             } else if(reference == "random2") {
##                 referenceI <- create_eq_ref(g)
##             }
##         } else {
##             referenceI <- reference
##             }
##         d_reference <- apply(m, 1, function(x) sum(abs(x - referenceI)))
##         f_det <- -c * d_reference
##         ## f_det <- rowSums(m) * slope/nrow(m) ## this is Greene and Krona
##         fi <- f_r + f_det

##         if(!is.null(scale)) {
##             fi <- (fi - min(fi))/(max(fi) - min(fi))
##             fi <- scale[1] + fi * (scale[2] - scale[1])
##         }
##         if(wt_is_1 == "divide") {
##             ## We need to shift to avoid ratios involving negative numbers and
##             ## we need to avoid having any fitness at 0, specially the wt.  If
##             ## any negative value, add the min, and shift by the magnitude of
##             ## the min to avoid any at 0.

##             ## If you use scale and wt_is_1, this will move the scale. It is
##             ## not possible to obtain a linear transformation that will keep
##             ## the min and max of the scale, and wt at 1.
##             min_fi <- min(fi)
##             if(min_fi < 0)
##                 fi <- fi + 2 * abs(min(fi))
##             fi <- fi/fi[1]
##         } else if (wt_is_1 == "subtract") {
##             fi <- fi - fi[1] + 1.0
##         } else if(wt_is_1 == "force") {
##             fi[1] <- 1.0
##             if(!is.null(scale)) {
##                 if( (1 > scale[2]) || (1 < scale[1]))
##                     warning("Using wt_is_1 = force and scale, but scale does ",
##                             "not include 1")
##             }
##         }
##         if(truncate_at_0) {
##             fi[fi <= 0] <- 1e-9
##         }
##         if(log) {
##             fi <- log(fi/fi[1]) + 1
##         }
##         m <- cbind(m, Fitness = fi)
##         if(!is_null_na(min_accessible_genotypes)) {
##             ## num_accessible_genotypes <- count_accessible_g(m, accessible_th)
##             ## Recall accessibleGenotypes includes the wt, if accessible.
##             num_accessible_genotypes <- length(wrap_accessibleGenotypes(m, accessible_th)) - 1
##             ## message("\n     num_accessible_genotypes = ", num_accessible_genotypes, "\n")
##             if(num_accessible_genotypes >= min_accessible_genotypes) {
##                 done <- TRUE
##                 attributes(m) <- c(attributes(m),
##                                    accessible_genotypes = num_accessible_genotypes,
##                                    accessible_th = accessible_th)
##             } else {
##                 ## Cannot start again with a fitness column
##                 m <- m[, -ncol(m), drop = FALSE]
##             }
##         } else {
##             done <- TRUE
##         }
##     }
##     ## message("\n number of attempts = ", attempts, "\n")
##     class(m) <- c(class(m), "genotype_fitness_matrix")
##     return(m)
## }

