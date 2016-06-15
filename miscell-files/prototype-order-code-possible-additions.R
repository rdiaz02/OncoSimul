## Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

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


## For now, I stop the work on the ordered part.

## This needs to be thought better, and there needs to be some demand.
## Things missing now:
##    1. Plotting fitness landscapes with ordered genotypes
##    2. Reading a fitness spec with order
##    3. Generating random fitness when order effects

## A general issue has to do with interactions. 

## Some of the code below is a nice idea, but can't really work because
## suppose A, B are ordered, and C is not.

## We get A, B, AB, BA, genotypes, each with and without C.

## But for coefficients, in OncoSimul, we can only specify A, B, BA, AB,
## and C, CA, CB, CAB. But CBA is missing from the previous, and cannot be
## specified in OncoSimulR. We would need to expand the model, and
## consider C as part of the order, and then make many terms identical
## (e.g., BAC same as BCA same as CBA).

## The problem is that order creates like new events (BA different from
## AB) and thus we can conceivably say that BAC has a different
## coefficient from ABC.


## What code needs to be written:

## 1. Plotting:

##    - the adjacency matrix dance.
##           - find genotypes with 1 more mutation (candidates)

##           - prune these by considering only those with same order up to
##             the new mutation. There is code below to see if an order
##             effect included. But we ask for more. So probably see if a
##             gsub(the_genotype, the_candidate) only returns " > a_gene"

##           - what do we do with epistasis and no interactions? Think this.


## 2. See comments above and the prototype code below.

##   An additional issue is that if we allow order, the nubmer of
##   genotypes blows up. Users might want to say "A > B _ C, D" which
##   means order affects A and B, but not C and D. If they pass that, then
##   we must neverthelss expand to a full thing with all four ordered, but this would be like the user said:

##   I am saying A > B and a C and a D has a given fitness. So
##   A > C > D > B  and C > A > B > D, etc, do have the same fitness.

## Can be accommodated the following way in the code below:


##    Take all genotypes given.
##    For each genotype, generate all permutations of the ordered and unordered parts.
##    Exclude those that have the wrong order for the originally ordered part.
##    Return the rest as genotypes with identical fitness.
##         How do I exclude/find out if order OK? See the order_term_present function below.
##     This becomes, say, the g2 list.
##    Finally, generate a huge list of ALL possible genotypes.: g3.
##    Assign to g3 all a fitness of 1.
##    Replace the fitness of g3 with those of g2, as appropriate.
##    Solve the system.





## 3. Should be easy. There is lots of code to generate all
## combinations. The problem is what do do next and if we really want all
## combinations.



## FIXME: is the way of getting the epist from fitness fully sound? I
## think so, but ... double check. (With tests)






## Where is the next function useful? To get the smallest set of genotypes
## for plotting landscapes, for instance. Or to obtain the list of
## genotypes.

generateAllGenotypes_minimal <- function(fitnessEffects, max = 256) {
    ## This generates the strictly needed set of genotypes: it uses order
    ## when it should, and does not when it shouldn't. And adapts usage to
    ## the different parts of the fitnessEffects specification.

    ## We really go for the minimal tables. What blows things up are the
    ## order effects, so we only grab them for genes in order effects, and
    ## we generate all permutations. And then do the Cartesian product
    ## with any other combinations of genes.

    ## Next one is different in generateAllGenotypes because there
    ## starting from 0 means adding 1 entry so not very relevant. Not
    ## multiplying.

    ## We are assuming that there are no interactions between order
    ## effects and epistasis. See below.
    
    tot_o <- function(n) {sum(sapply(c(0, seq.int(n)),
                                     function(x) choose(n, x) * factorial(x)))}
    tot_no <- function(n) {2^n}

    order <- order_no_int <- no_order <- FALSE

   
    if(!is.null(fitnessEffects$orderEffects)) {
        order <- TRUE
        ## Only modules involved in order effects
        all_o_modules <- unique(unlist(lapply(fitnessEffects$long.orderEffects,
                                              function(x) x$NumID)))
        ## ... but what we need are the genes
        all_o_genes <- fitnessEffects$geneModule$GeneNumID[
            which(fitnessEffects$geneModule$ModuleNumID %in% all_o_modules )]

        ## n_o <- nrow(fitnessEffects$geneModule) -1
        n_o <- length(all_o_genes)
        ## There could be genes in epist, etc, not in order
        all_no_genes <- setdiff(fitnessEffects$geneModule$GeneNumID[-1], all_o_genes)
        n_no <- length(all_no_genes)
        
        if(!is.null(fitnessEffects$noIntGenes)) {
            order_no_int <- TRUE
            n_no <- n_no + nrow(fitnessEffects$long.geneNoInt)
        }
        tnn <- tot_o(n_o) * tot_no(n_no)
    } else {
        no_order <- TRUE
        n_no <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
        tnn <- tot_no(n_no)
    }

    if(tnn > max) {
        m <- paste("There are", tnn, "genotypes.")
        m <- paste(m, "This is larger than max.")
        m <- paste(m, "Adjust max and rerun if you want")
        stop(m)
    }

    vido <- vidno <- vector("integer")
    if(order) {
        ## vido <- fitnessEffects$geneModule$GeneNumID[-1] ## 0 not in there
        vido <- all_o_genes
        vidno <- all_no_genes
        if(order_no_int)
            vidno <- c(vidno, fitnessEffects$long.geneNoInt$GeneNumID)
    } else if(no_order) {
        vidno <- allNamedGenes(fitnessEffects)$GeneNumID
    }

    ## Yes, I keep thinking there's got to be a simpler way :-)
    ## Anyway, we obtain a list of the vectors
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    
    ## Function to call with ordered part
    f1_o <- function(n) {
        lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vido))}
    ## Ditto in non-ordered part
    f1_no <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vidno))}

    names <- c(fitnessEffects$geneModule$Gene[-1],
               fitnessEffects$long.geneNoInt$Gene)

    f_names_o <- function(ggnums) {
        unlist(lapply(lapply(ggnums, function(x) names[x]),
                      function(z)
                          paste(z,
                                collapse = " > " )))
    }
    f_names_no <- function(ggnums) {
        unlist(lapply(lapply(ggnums, function(x) names[x]),
                      function(z)
                          paste(z,
                                collapse = ", " )))
    }
    
    my.na.omit <- function(x) {
        ## I na.omit a lot but don't need the attributes
        z <- na.omit(x)
        attributes(z)$na.action <- NULL
        z
    }

    if(order) {
        genotNums_o <- f1_o(n_o)
        genotNums_o <- list.of.vectors(genotNums_o)
        genotNames_o <- f_names_o(genotNums_o)
        if(n_no > 0) {
        ## if(order_no_int) {
            genotNums_no <- f1_no(n_no)
            genotNums_no <- list.of.vectors(genotNums_no)
            genotNames_no <- f_names_no(genotNums_no)
            ## No ordered mutations or no non-ordered mutations: The WT is
            ## added by fiat. Not the ones where only the ordered XOR the
            ## non-ordered part is not present. So we add NAs, and then
            ## remove them where needed. We also remove the very first
            ## entry, which is the WT proper.
            genotNums_no <- c(NA, genotNums_no)
            genotNums_o <- c(NA, genotNums_o)
            genotNames_no <- c(NA, genotNames_no)
            genotNames_o <- c(NA, genotNames_o)
            ## Combine them
            no_part <- rep(genotNums_no, rep(length(genotNums_o),
                                           length(genotNums_no)))
            o_part <- rep(genotNums_o, length(genotNums_no))
            genotNums_o <- lapply(mapply(c, o_part, no_part), my.na.omit)[-1]
            genotNames_o <- apply(
                expand.grid(genotNames_o, genotNames_no, stringsAsFactors = FALSE), 1,
                function(z) gsub("NA", "", paste(z, collapse = " _ ")))[-1]
            rm(genotNames_no, genotNums_no) ## paranoid removal
        }
    } else if(no_order) {
        genotNums_no <- f1_no(n_no)
        genotNums_no <- list.of.vectors(genotNums_no)
        genotNames_no <- f_names_no(genotNums_no)
    }

    ## This ain't the best, as we repeatedly read and convert
    ## fitnessEffects.  If this were slow, prepare C++ function that takes
    ## vectors of genotypes and uses same fitnessEffects.

    return(list(genotNums = if(order) {genotNums_o} else{genotNums_no},
                genotNames = if(order) {genotNames_o} else{genotNames_no}))

    ## ## With mutator, the ids of genes need not go from 1:n
    ## vid <- allNamedGenes(fitnessEffects)$GeneNumID
    ## if(order) {
    ##     f1 <- function(n) {
    ##         lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vid))
    ##     }
    ## } else {
    ##     f1 <- function(n) {
    ##         lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vid))}
        
    ## }
   
    ## genotNums <- f1(nn)
    ## list.of.vectors <- function(y) {
    ##     ## there's got to be a simpler way
    ##     lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
    ##            function(m) m[[1]])
    ## }
    ## genotNums <- list.of.vectors(genotNums)
    ## names <- c(fitnessEffects$geneModule$Gene[-1],
    ##            fitnessEffects$long.geneNoInt$Gene)
    
    ## genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]),
    ##                             function(z)
    ##                                 paste(z,
    ##                                       collapse = if(order){" > "} else {", "} )))
    ## This ain't the best, as we repeatedly read and convert
    ## fitnessEffects.  If this were slow, prepare C++ function that takes
    ## vectors of genotypes and uses same fitnessEffects.
    ## return(list(genotNums = genotNums,
    ##             genotNames = genotNames
    ##             ))
}



## generateAllGenotypes_minimal <- function(fitnessEffects, max = 256) {
##     ## This generates the strictly needed set of genotypes: it uses order
##     ## when it should, and does not when it shouldn't. And adapts usage to
##     ## the different parts of the fitnessEffects specification.

##     ## We use a simplifying assumption that might give somewhat larger
##     ## tables: if order is ever used, then we obtain all ordered genotypes
##     ## for everything not in geneNoInt. And geneNoInt as unordered. And
##     ## then we do the Cartesian product.

##     ## Why shouldn't we try to simplify any further? Because users can
##     ## pass the same genes in order and epistasis as well as partially
##     ## overlapping lists.

##     ## Next one is different in generateAllGenotypes because there
##     ## starting from 0 means adding 1 entry so not very relevant. Not
##     ## multiplying.
##     tot_o <- function(n) {sum(sapply(c(0, seq.int(n)),
##                                      function(x) choose(n, x) * factorial(x)))}
##     tot_no <- function(n) {2^n}

##     order <- order_no_int <- no_order <- FALSE
    
##     if(!is.null(fitnessEffects$orderEffects)) {
##         order <- TRUE
##         n_o <- nrow(fitnessEffects$geneModule) -1

##         if(!is.null(fitnessEffects$noIntGenes)) {
##             order_no_int <- TRUE
##             n_no <- nrow(fitnessEffects$long.geneNoInt)
##         }
##         tnn <- tot_o(n_o) * tot_no(n_no)
##     } else {
##         no_order <- TRUE
##         n_no <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
##         tnn <- tot_no(n_no)
##     }

##     if(tnn > max) {
##         m <- paste("There are", tnn, "genotypes.")
##         m <- paste(m, "This is larger than max.")
##         m <- paste(m, "Adjust max and rerun if you want")
##         stop(m)
##     }

##     vido <- vidno <- vector("integer")
##     if(order) {
##         vido <- fitnessEffects$geneModule$GeneNumID[-1] ## 0 not in there
##         if(order_no_int)
##             vidno <- fitnessEffects$long.geneNoInt$GeneNumID
##     } else if(no_order) {
##         vidno <- allNamedGenes(fitnessEffects)$GeneNumID
##     }

##     ## Yes, I keep thinking there's got to be a simpler way :-)
##     ## Anyway, we obtain a list of the vectors
##     list.of.vectors <- function(y) {
##         ## there's got to be a simpler way
##         lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
##                function(m) m[[1]])
##     }
    
##     ## Function to call with ordered part
##     f1_o <- function(n) {
##         lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vido))}
##     ## Ditto in non-ordered part
##     f1_no <- function(n) {
##         lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vidno))}

##     names <- c(fitnessEffects$geneModule$Gene[-1],
##                fitnessEffects$long.geneNoInt$Gene)

##     f_names_o <- function(ggnums) {
##         unlist(lapply(lapply(ggnums, function(x) names[x]),
##                       function(z)
##                           paste(z,
##                                 collapse = " > " )))
##     }
##     f_names_no <- function(ggnums) {
##         unlist(lapply(lapply(ggnums, function(x) names[x]),
##                       function(z)
##                           paste(z,
##                                 collapse = ", " )))
##     }
    
##     my.na.omit <- function(x) {
##         ## I na.omit a lot but don't need the attributes
##         z <- na.omit(x)
##         attributes(z)$na.action <- NULL
##         z
##     }

##     if(order) {
##         genotNums_o <- f1_o(n_o)
##         genotNums_o <- list.of.vectors(genotNums_o)
##         genotNames_o <- f_names_o(genotNums_o)
##         if(order_no_int) {
##             genotNums_no <- f1_no(n_no)
##             genotNums_no <- list.of.vectors(genotNums_no)
##             genotNames_no <- f_names_no(genotNums_no)
##             ## No ordered mutations or no non-ordered mutations: The WT is
##             ## added by fiat. Not the ones where only the ordered XOR the
##             ## non-ordered part is not present. So we add NAs, and then
##             ## remove them where needed. We also remove the very first
##             ## entry, which is the WT proper.
##             genotNums_no <- c(NA, genotNums_no)
##             genotNums_o <- c(NA, genotNums_o)
##             genotNames_no <- c(NA, genotNames_no)
##             genotNames_o <- c(NA, genotNames_o)
##             ## Combine them
##             no_part <- rep(genotNums_no, rep(length(genotNums_o),
##                                            length(genotNums_no)))
##             o_part <- rep(genotNums_o, length(genotNums_no))
##             genotNums_o <- lapply(mapply(c, o_part, no_part), my.na.omit)[-1]
##             genotNames_o <- apply(
##                 expand.grid(genotNames_o, genotNames_no, stringsAsFactors = FALSE), 1,
##                 function(z) gsub("NA", "", paste(z, collapse = " _ ")))[-1]
##             rm(genotNames_no, genotNums_no) ## paranoid removal
##         }
##     } else if(no_order) {
##         genotNums_no <- f1_no(n_no)
##         genotNums_no <- list.of.vectors(genotNums_no)
##         genotNames_no <- f_names_no(genotNums_no)
##     }

##     ## This ain't the best, as we repeatedly read and convert
##     ## fitnessEffects.  If this were slow, prepare C++ function that takes
##     ## vectors of genotypes and uses same fitnessEffects.

##     return(list(genotNums = if(order) {genotNums_o} else{genotNums_no},
##                 genotNames = if(order) {genotNames_o} else{genotNames_no}))

##     ## ## With mutator, the ids of genes need not go from 1:n
##     ## vid <- allNamedGenes(fitnessEffects)$GeneNumID
##     ## if(order) {
##     ##     f1 <- function(n) {
##     ##         lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vid))
##     ##     }
##     ## } else {
##     ##     f1 <- function(n) {
##     ##         lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vid))}
        
##     ## }
   
##     ## genotNums <- f1(nn)
##     ## list.of.vectors <- function(y) {
##     ##     ## there's got to be a simpler way
##     ##     lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
##     ##            function(m) m[[1]])
##     ## }
##     ## genotNums <- list.of.vectors(genotNums)
##     ## names <- c(fitnessEffects$geneModule$Gene[-1],
##     ##            fitnessEffects$long.geneNoInt$Gene)
    
##     ## genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]),
##     ##                             function(z)
##     ##                                 paste(z,
##     ##                                       collapse = if(order){" > "} else {", "} )))
##     ## This ain't the best, as we repeatedly read and convert
##     ## fitnessEffects.  If this were slow, prepare C++ function that takes
##     ## vectors of genotypes and uses same fitnessEffects.
##     ## return(list(genotNums = genotNums,
##     ##             genotNames = genotNames
##     ##             ))
## }




ordered_part <- function(x) {

}


get_the_two_parts <- function(x) {
    tmp <- strsplit(x$Genotype, "_", fixed = TRUE)
    return(list(ordered = lapply(tmp, function(z) z[1]),
                unordered = lapply(tmp, function(z) z[2])
                ))
}

order_term_present <- function(O, G) {
    ## Similar logic to my match_order_effects C++ function
    ## Check if the order term in this genotype
    if(length(G) < length(O)) return(0L)
    mm <- match(O, G)
    if(any(is.na(mm))) return(0L)
    if(is.unsorted(mm)) {
        return(0L)
    } else {
        return(1L)
    }
}

unordered_term_present <- function(T, G) {
    ## Cannot tell if epistasis or something else
    if(length(G) < length(T)) return(0L)
    mm <- match(T, G)
    if(any(is.na(mm))) {
        return(0L)
    } else {
        return(1L)
    }
}



## FIXME check what happens with some fitness = 0

## We assume a term is either in the order or in the no order part. Not in
## both. If in both, there are many mappings possible

## If anything is ordered, everything is ordered. There is no simple way
## around that. See below.



## The general idea here is to use a system of equations and solve:
## log fitness = x*a
## the coefficients are for log(1 + si)
## we solve, so we then obtaini the si as exp(log(1 + si)) -1


genot_fitness_to_order <- function(x) {
    anywt <- which(x[, 1] == "WT")
    if(length(anywt) > 1) stop("More than 1 WT")
    if(length(anywt) == 1) {
        fwt <- x[anywt, 2]
        x <- x[-anywt, ]
    } else {
        fwt <- 1
    }

    if(any(x$Fitness) < 0)
        stop("Fitness cannot be less than 0")
    ## FIXME if there is a WT, use it to set fitness by making it a divisor

    two_parts <- get_the_two_parts(x)
    
    splitted_genots_ordered <- lapply(two_parts$ordered,
                                      function(z) OncoSimulR:::nice.vector.eo(z, ">"))
    splitted_genots_unordered <- lapply(two_parts$unordered,
                                      function(z) OncoSimulR:::nice.vector.eo(z, ","))

    ## And put all of them together. Period.
    all_o_genes <- na.omit(unique(unlist(splitted_genots_ordered)))
    all_no_genes <- na.omit(unique(unlist(splitted_genots_unordered)))
    attributes(all_o_genes)$na.action <- NULL
    attributes(all_no_genes)$na.action <- NULL

    f1_o <- function(n) {
        lapply(seq.int(n), function(x) permutations(n = n, r = x, v = all_o_genes))}
    f1_no <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x, v = all_no_genes))}
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    ordered_terms <- list.of.vectors(f1_o(length(all_o_genes)))
    unordered_terms <- list.of.vectors(f1_no(length(all_no_genes)))

    ## Create the set of all possible genotypes
    ## First is for checking against coefficients
    ordered_terms <- c(NA, ordered_terms)
    unordered_terms <- c(NA, unordered_terms)
    no_part <- rep(unordered_terms, rep(length(ordered_terms),
                                     length(unordered_terms)))
    o_part <- rep(ordered_terms, length(unordered_terms))
    all.genotypes <- mapply(c, o_part, no_part)

    ## FIXME we need more terms for "unordered". The interaction of any
    ## unordered with any ordered, as a single thing. Otherwise, I have
    ## fewer interactions than I should.

    ## Make a note of this: if A, B are ordered, and C is unordered, in
    ## addition to A, B, AB, BA, C, CA, CB I need CAB and CBA (note, CAB
    ## and CBA)

    
    if(nrow(x) > length(ordered_terms) * length(unordered_terms))
        stop("You must have repeated entries or a problem of specification."
             " There are more entries than (sum of permutations) + (sum of combinations).")

    ncoeff <- length(ordered_terms) + length(unordered_terms)  ## remove the constant
    coeffs_ordered <- matrix(0L, nrow = ncoeff, ncol = length(ordered_terms))
    coeffs_unordered <- matrix(0L, nrow = ncoeff, ncol = length(unordered_terms))

    for(i in seq.int(ncoeff)) {
        coeffs_ordered[i, ] <- unlist(lapply(ordered_terms,
                                             function(x)
                                                 order_term_present(x,
                                                                    o_part[[i]])))
        coeffs_unordered[i, ] <- unlist(lapply(unordered_terms,
                                               function(x)
                                                   unordered_term_present(x,
                                                                          no_part[[i]])))
    }

  
    ## This is for later taking the diff w.r.t. given ones
    pasted_ordered_terms <- unlist(lapply(ordered_terms, function(x) paste(x, collapse = " > " ) ))
    pasted_unordered_terms <- unlist(lapply(unordered_terms, function(x) paste(x, collapse = " , " ) ))
    all.genotypes.string <- apply(
        expand.grid(pasted_ordered_terms, pasted_unordered_terms, stringsAsFactors = FALSE), 1,
        function(z) gsub("NA", "", paste(z, collapse = " _ ")))

    ## We need to do it again, because of the possible formatting differences
    gg_o <- unlist(lapply(splitted_genots_ordered, function(x) paste(x, collapse = " > " ) ))
    gg_no <- unlist(lapply(splitted_genots_unordered, function(x) paste(x, collapse = " , " ) ))
    these.genotypes.string <- paste(gsub("NA", "", gg_o), gsub("NA", "", gg_no), sep = " _ ")

    coeffs <- cbind(coeffs_ordered, coeffs_unordered)
    colnames(coeffs) <- c(pasted_ordered_terms, pasted_unordered_terms)

    log_fitness <- rep(0, ncoeff)
    log_fitness[match( these.genotypes.string, all.genotypes.string)] <- log(x$Fitness)
        
    
    solve(coeffs, log_fitness)
    
    
    
    
    num_genes <- unlist(lapply(splitted_genots, length))
    total_genes <- length(unique(unlist(splitted_genots)))
    total_genots <- factorial(total_genes)
    

    ## Here is where we do the heavy lifting:
    ## order rows by number of genes
    ## get coefficients of lower order
    ## do the math
    ## repeat

    if(nrow(x) < factorial(total_genots))
        warning("Number of genotypes < factorial(number of genes). ",
                "This is unlikely to do what you want, ",
                "and you'll have to do the math to understand the fitness ",
                "of the specified and the missing genotypes",
                " (we cannot guess your intentions).")
        
    x <- x[order(num_genes),  ]
    if(!isTRUE(all.equal(fwt, 1))) {
        message("Fitness of wildtype != 1. ",
                "Dividing all fitnesses by fitness of wildtype.")
        x[, 2] <- x[, 2]/fwt
    }
    s <- x[, 2] - 1
    names(s) <- x[, 1]
    return(s)
}


## The code below is a nice idea, but can't really work because suppose A,
## B are ordered, and C is not.

## We get A, B, AB, BA, genotypes, each with and without C.

## But for coefficients, in OncoSimul, we can only specify A, B, BA, AB,
## and C, CA, CB, CAB. But CBA is missing from the previous, and cannot be
## specified in OncoSimulR. We would need to expand the model, and
## consider C as part of the order, and then make many terms identical
## (e.g., BAC same as BCA same as CBA).

## The problem is that order creates like new events (BA different from
## AB) and thus we can conceivably say that BAC has a different
## coefficient from ABC.



genot_fitness_to_order <- function(x) {
    anywt <- which(x[, 1] == "WT")
    if(length(anywt) > 1) stop("More than 1 WT")
    if(length(anywt) == 1) {
        fwt <- x[anywt, 2]
        x <- x[-anywt, ]
    } else {
        fwt <- 1
    }

    if(any(x$Fitness) < 0)
        stop("Fitness cannot be less than 0")
    ## FIXME if there is a WT, use it to set fitness by making it a divisor

    two_parts <- get_the_two_parts(x)
    
    splitted_genots_ordered <- lapply(two_parts$ordered,
                                      function(z) OncoSimulR:::nice.vector.eo(z, ">"))
    splitted_genots_unordered <- lapply(two_parts$unordered,
                                      function(z) OncoSimulR:::nice.vector.eo(z, ","))
   
    all_o_genes <- na.omit(unique(unlist(splitted_genots_ordered)))
    all_no_genes <- na.omit(unique(unlist(splitted_genots_unordered)))
    attributes(all_o_genes)$na.action <- NULL
    attributes(all_no_genes)$na.action <- NULL

    f1_o <- function(n) {
        lapply(seq.int(n), function(x) permutations(n = n, r = x, v = all_o_genes))}
    f1_no <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x, v = all_no_genes))}
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    ordered_terms <- list.of.vectors(f1_o(length(all_o_genes)))
    unordered_terms <- list.of.vectors(f1_no(length(all_no_genes)))

    ## Create the set of all possible genotypes
    ## First is for checking against coefficients
    ordered_terms <- c(NA, ordered_terms)
    unordered_terms <- c(NA, unordered_terms)
    no_part <- rep(unordered_terms, rep(length(ordered_terms),
                                     length(unordered_terms)))
    o_part <- rep(ordered_terms, length(unordered_terms))
    all.genotypes <- mapply(c, o_part, no_part)

    ## FIXME we need more terms for "unordered". The interaction of any
    ## unordered with any ordered, as a single thing. Otherwise, I have
    ## fewer interactions than I should.

    ## Make a note of this: if A, B are ordered, and C is unordered, in
    ## addition to A, B, AB, BA, C, CA, CB I need CAB and CBA (note, CAB
    ## and CBA)

    
    if(nrow(x) > length(ordered_terms) * length(unordered_terms))
        stop("You must have repeated entries or a problem of specification."
             " There are more entries than (sum of permutations) + (sum of combinations).")

    ncoeff <- length(ordered_terms) + length(unordered_terms)  ## remove the constant
    coeffs_ordered <- matrix(0L, nrow = ncoeff, ncol = length(ordered_terms))
    coeffs_unordered <- matrix(0L, nrow = ncoeff, ncol = length(unordered_terms))

    for(i in seq.int(ncoeff)) {
        coeffs_ordered[i, ] <- unlist(lapply(ordered_terms,
                                             function(x)
                                                 order_term_present(x,
                                                                    o_part[[i]])))
        coeffs_unordered[i, ] <- unlist(lapply(unordered_terms,
                                               function(x)
                                                   unordered_term_present(x,
                                                                          no_part[[i]])))
    }

  
    ## This is for later taking the diff w.r.t. given ones
    pasted_ordered_terms <- unlist(lapply(ordered_terms, function(x) paste(x, collapse = " > " ) ))
    pasted_unordered_terms <- unlist(lapply(unordered_terms, function(x) paste(x, collapse = " , " ) ))
    all.genotypes.string <- apply(
        expand.grid(pasted_ordered_terms, pasted_unordered_terms, stringsAsFactors = FALSE), 1,
        function(z) gsub("NA", "", paste(z, collapse = " _ ")))

    ## We need to do it again, because of the possible formatting differences
    gg_o <- unlist(lapply(splitted_genots_ordered, function(x) paste(x, collapse = " > " ) ))
    gg_no <- unlist(lapply(splitted_genots_unordered, function(x) paste(x, collapse = " , " ) ))
    these.genotypes.string <- paste(gsub("NA", "", gg_o), gsub("NA", "", gg_no), sep = " _ ")

    coeffs <- cbind(coeffs_ordered, coeffs_unordered)
    colnames(coeffs) <- c(pasted_ordered_terms, pasted_unordered_terms)

    log_fitness <- rep(0, ncoeff)
    log_fitness[match( these.genotypes.string, all.genotypes.string)] <- log(x$Fitness)
        
    
    solve(coeffs, log_fitness)
    
    
    
    
    num_genes <- unlist(lapply(splitted_genots, length))
    total_genes <- length(unique(unlist(splitted_genots)))
    total_genots <- factorial(total_genes)
    

    ## Here is where we do the heavy lifting:
    ## order rows by number of genes
    ## get coefficients of lower order
    ## do the math
    ## repeat

    if(nrow(x) < factorial(total_genots))
        warning("Number of genotypes < factorial(number of genes). ",
                "This is unlikely to do what you want, ",
                "and you'll have to do the math to understand the fitness ",
                "of the specified and the missing genotypes",
                " (we cannot guess your intentions).")
        
    x <- x[order(num_genes),  ]
    if(!isTRUE(all.equal(fwt, 1))) {
        message("Fitness of wildtype != 1. ",
                "Dividing all fitnesses by fitness of wildtype.")
        x[, 2] <- x[, 2]/fwt
    }
    s <- x[, 2] - 1
    names(s) <- x[, 1]
    return(s)
}




## A test
## B involved in both epi and order
## D only in order
## F, G, H only in epi
## A and E no interaction
foi3 <- allFitnessEffects(
    orderEffects = c("D>B" = -0.2, "B > D" = 0.3),
    epistasis = c("B:F" = 0.12, "F:H" = 0.14),
    noIntGenes = c("A" = 0.05),
    geneToModule = c("D"= "d1", "B" = "b1, b2",
                     "F" = "f1, f2", "H" = "h"))

generateAllGenotypes_minimal(foi3, max = 256)

## We obtain 255 = 256 - 1 (no wildtype)







m2 <- data.frame(Genotype = c("A > C _ E, F", "C > B _ G",
                              "_ G", "C _"),
                 Fitness = c(1, 2, 3, 4),
                 stringsAsFactors = FALSE)



m2s <- data.frame(Genotype = c("A > B _ c, d", " > B _ c",
                               "_ d", "B _"),
                  Fitness = c(1, 2, 3, 4),
                  stringsAsFactors = FALSE)
