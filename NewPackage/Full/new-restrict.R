# - Say that a user can use a "0" as a gene name, but that is BAD idea.
## - Modules and order effects can be kind of funny?
library(data.table)
library(Rcpp)
library(gtools) ## for permutations
## setwd("../../")

sourceCpp("new-restrict.cpp", verbose = TRUE)
## sourceCpp("t1.cpp", verbose = TRUE)


check.gm <- function(gm) {
    ## Yes, Root: we want no ambiguities
    if(gm[1] != "Root")
        stop("First value of a module table must be Root")
    if(names(gm)[1] != "Root")
        stop("First name of a module table must be Root")
    if(length(unique(names(gm))) != length(gm))
        stop("Number of unique module names different from length of vector")
    
}

gtm2 <- function(x) {
    data.frame(cbind(nice.vector.eo(x, ","), x))
}

nice.vector.eo <- function(z, sep) {
    ## with epistasis, maybe we want sorted?
    setdiff(unlist(lapply(strsplit(z, " "),
                                    function(u) strsplit(u, sep))), "")
}

gm.to.geneModuleL <- function(gm, one.to.one) {
    ## the table will be sorted by gene name
    check.gm(gm)
    ## the named vector with the mapping into the long geneModule df
    geneMod <- as.data.frame(rbindlist(lapply(gm, gtm2)))
    geneMod$Module <- names(gm)[geneMod[, 2]] ## reverse lookup table
    colnames(geneMod)[1] <- c("Gene")
    geneMod <- geneMod[, -2]
    geneMod$Gene <- as.character(geneMod$Gene)
    ## geneMod$Module <- as.character(geneMod$Module) ## already a char
    geneMod <- geneMod[c(1, order(geneMod$Gene[-1]) + 1), ] 
    geneMod$GeneNumID <- 0:(nrow(geneMod) - 1)

    ## this assumes sorted! and they need not be
    ## rlem <- rle(geneMod$Module)
    ## geneMod$ModuleNumID <- rep( 0:(length(rlem$values) - 1), rlem$lengths)
    if(one.to.one) {
        geneMod$ModuleNumID <- geneMod$GeneNumID
        if(!(identical(geneMod$Gene, geneMod$Module)))
            stop("Impossible: we are supposed to be in one-to-one for Module-Gene.")
    } else {
        idm <- seq.int(length(gm) - 1)
        idm <- c("Root" = 0L, idm)
        names(idm) <- names(gm)
        geneMod$ModuleNumID <- idm[geneMod[, "Module"]]
    }
    rownames(geneMod) <- 1:nrow(geneMod)
    geneMod   
}

geneModuleNull <- function(namesM) {
    v <- c("Root", setdiff(namesM, "Root"))
    names(v) <- v
    return(v)
}

list.of.deps <- function(x) {
    ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
    ##                 "SM" = 2, "semimonotone" = 2)
    lookupTypeDep <- c("MN" = "monotone",
                       "monotone" = "monotone",
                       "SM" = "semimonotone",
                       "semimonotone" = "semimonotone",
                       "XOR" = "xmpn",
                       "xmpn" = "xmpn",
                       "XMPN" = "xmpn")
    ## FIXME: check values of typeDep

    if(length(x) > 1) {
        if(length(unique(x$s))!= 1)
            stop("Not all s identical within a child")
        if(length(unique(x$sh))!= 1)
            stop("Not all sh identical within a child")
        if(length(unique(x$typeDep))!= 1)
            stop("Not all typeDep identical within a child")
        if(length(unique(x$child))!= 1)
            stop("child not unique")
    }
    return(list(
        child = unique(x$child),
        s = unique(x$s),
        sh = unique(x$sh),
        typeDep = lookupTypeDep[unique(x$typeDep)],
        parents = unlist(x$parent)))

}

to.long.rt <- function(rt, idm, verbosity = 0) {
    if(is.numeric(rt$parent))
        rt$parent <- as.character(rt$parent)
    if(!("Root" %in% rt$parent))
        stop("Root must be one parent node")
    if(is.numeric(rt$child))
        rt$child <- as.character(rt$child)
    ## rt$parent <- unlist(lapply(rt$parent, nice.string))
    ## rt$child <- unlist(lapply(rt$child, nice.string))
   
    srt <- rt[order(rt$child), ]

    ## Not relevant if we allow non-numeric names
    ## all.child.genes <- as.integer(
    ##     unlist(lapply(rt[, 2],
    ##                   function(x) strsplit(x, ","))))
    ## ## check all childs
    ## if(!identical(sort(unique(all.child.genes)),
    ##               seq.int(max(all.child.genes))))
    ##     stop("Not all children present")
    long.rt <- lapply(split(srt, srt$child), list.of.deps)

    ## geneModule <- gene.to.module(srt)
    ## idm <- seq.int(length(names(long.rt)))
    ## names(idm) <- names(long.rt)
    ## idm <- c("0" = 0L, idm)
    ## geneModule$ModuleNumID <- idm[geneModule[, "Module"]]

    ## idm is just a look up table for the id of the module
    ## idm <- unique(geneModule$ModuleNumID)
    ## names(idm) <- unique(geneModule$Module)
    
    ## add integer IDs
    addIntID <- function(z, idm) {
        z$childNumID <- idm[z$child]
        z$parentsNumID <- idm[z$parents]
        if( any(is.na(z$parentsNumID)) ||
           any(is.na(z$childNumID)) ) {
            stop(paste("An ID is NA:",
                       "Is a gene part of two different modules?",
                       "(That includes being by itself and part",
                       "of a module.)"))
            
        }
        ## sort the parents, always.
        o <- order(z$parentsNumID)
        z$parentsNumID <- z$parentsNumID[o]
        z$parents <- as.character(z$parents[o])
        return(z)
    }
    long.rt <- lapply(long.rt, function(x) addIntID(x, idm = idm))
   
    if(verbosity >= 4) {
        message(paste("Number of drivers: ",
                      length(unique(geneModule[, "Gene"]))))
        message(paste("Number of modules: ",
                      length(unique(geneModule[, "Module"]))))
    }
    return(long.rt)
    ## return(list(long.rt = long.rt, geneModule = geneModule))
}


epist.order.element <- function(x, y, sep) {
    list(ids = nice.vector.eo(x, sep = sep), s = y)
}


to.long.epist.order <- function(epor, sep) {
    if(is.vector(epor))
        long <- Map(function(x, y) epist.order.element(x, y, sep),
                       names(epor), epor)
    else if(is.data.frame(epor)) 
        long <- Map(function(x, y) epist.order.element(x, y, sep),
                    as.character(epor$ids),
                    epor$s)
    names(long) <- NULL
    return(long)
}

## addIntID.epist.order <- function(z, idm, sort) {
##     z$NumID <- idm[z$ids]
##     if(sort) {
##         ## essential for epistasis, but never do for order effects
##         o <- order(z$NumID)
##         z$NumID <- z$NumID[o]
##         z$ids <- z$ids[o]
##     }
##     return(z)
## }


addIntID.epist.order <- function(z, idm, sort, sign) {
    if( sort && (!sign))
        warning("sort is TRUE but sign is FALSE. You sure?")
    if((!sort) && sign)
        warning("sort is FALSE but sign is TRUE. You sure?")
    ## Adds the integer, but takes care of sign if sign = TRUE
    if(!sign) {
        z$NumID <- idm[z$ids]
        signs <- grep("^-", z$ids)
        if(length(signs))
            stop("You have a negative sign and you said sign = FALSE")
    } else {
        unsigned <- setdiff(unlist(lapply(z$ids,
                                          function(z) strsplit(z, "^-"))), "")
        NumID <- idm[unsigned]
        signs <- grep("^-", z$ids)
        if(length(signs)) {
            NumID[signs] <- -1 * NumID[signs]
        }
        z$NumID <- NumID
    }
    if(sort) {
        ## Essential for epistasis, but never do for order effects
        o <- order(z$NumID)
        z$NumID <- z$NumID[o]
        z$ids <- z$ids[o]
    }
    return(z)
}


checkRT <- function(mdeps) {
    if(ncol(mdeps) != 5)
        stop("mdeps must be of exactly 5 columns")
    if(!identical(colnames(mdeps), c("parent", "child", "s", "sh", "typeDep")))
        stop(paste("Column names of mdeps not of appropriate format. ",
                   "Should be parent, child, s, sh, typeDep"))
}





allFitnessEffects <- function(rT = NULL,
                              epistasis = NULL,
                              orderEffects = NULL,
                              noIntGenes = NULL,
                              geneToModule = NULL) {
    ## restrictions: the usual rt

    ## epistasis: as it says, with the ":"

    ## orderEffects: the ">"
    
    ## All of the above can be genes or can be modules (if you pass a
    ## geneToModule)

    ## rest: rest of genes, with fitness


    ## For epistasis and order effects we create the output object but
    ## missing the numeric ids of genes. With rT we do it in one go, as we
    ## already know the mapping of genes to numeric ids. We could do the
    ## same in epistasis and order, but we would be splitting twice
    ## (whereas for rT extracting the names is very simple).

    
    rtNames <- NULL
    epiNames <- NULL
    orNames <- NULL
    if(!is.null(rT))
        rtNames <- unique(c(rT$parent, rT$child))

    if(!is.null(epistasis)) {
        long.epistasis <- to.long.epist.order(epistasis, ":")
        ## epiNames <- unique(unlist(lapply(long.epistasis, function(x) x$ids)))
        ## deal with the possible negative signs
        epiNames <- setdiff(unique(
            unlist(lapply(long.epistasis,
                          function(x) lapply(x$ids,
                                             function(z) strsplit(z, "^-"))))),
                            "")
        
    } else {
        long.epistasis <- list()
    }
    if(!is.null(orderEffects)) {
        long.orderEffects <- to.long.epist.order(orderEffects, ">")
        orNames <- unique(unlist(lapply(long.orderEffects, function(x) x$ids)))
    } else {
        long.orderEffects <- list()
    }
    allModuleNames <- unique(c(rtNames, epiNames, orNames))
    if(is.null(geneToModule)) {
        gMOneToOne <- TRUE
        geneToModule <- geneModuleNull(allModuleNames)
    } else {
        gMOneToOne <- FALSE
        if(any(is.na(match(setdiff(names(geneToModule), "Root"), allModuleNames))))
            stop(paste("Some values in geneToModule not present in any of",
                       " rT, epistasis, or order effects"))
        if(any(is.na(match(allModuleNames, names(geneToModule)))))
            stop(paste("Some values in rT, epistasis, ",
                       "or order effects not in geneToModule"))
    }
    geneModule <- gm.to.geneModuleL(geneToModule, one.to.one = gMOneToOne)

    idm <- unique(geneModule$ModuleNumID)
    names(idm) <- unique(geneModule$Module)

    if(!is.null(rT)) {
        checkRT(rT)
        long.rt <- to.long.rt(rT, idm)
    } else {
        long.rt <- list() ## yes, we want an object of length 0
    }

    ## Append the numeric ids to epistasis and order
    if(!is.null(epistasis)) {
        long.epistasis <- lapply(long.epistasis,
                                 function(x)
                                     addIntID.epist.order(x, idm,
                                                          sort = TRUE,
                                                          sign = TRUE))
    }
    if(!is.null(orderEffects)) {
        long.orderEffects <- lapply(long.orderEffects,
                                    function(x)
                                        addIntID.epist.order(x, idm,
                                                             sort = FALSE,
                                                             sign = FALSE))
    }
    
    
    if(!is.null(noIntGenes)) {
        mg <- max(geneModule[, "GeneNumID"])
        gnum <- seq_along(noIntGenes) + mg
        if(!is.null(names(noIntGenes))) {
            ng <- names(noIntGenes)
            if(any(ng %in% geneModule[, "Gene"] ))
                stop("A gene in noIntGenes also present in the other terms")
        } else {
            ng <- gnum
        }
        geneNoInt <- data.frame(Gene = as.character(ng),
                                GeneNumID = gnum,
                                s = noIntGenes,
                                stringsAsFactors = FALSE)
    } else {
        geneNoInt <- data.frame()
    }
    if( (length(long.rt) + length(long.epistasis) + length(long.orderEffects) +
             nrow(geneNoInt)) == 0)
        stop("You have specified nothing!")
    out <- list(long.rt = long.rt,
                long.epistasis = long.epistasis,
                long.orderEffects = long.orderEffects,
                long.geneNoInt = geneNoInt,
                geneModule = geneModule,
                gMOneToOne = gMOneToOne)
    class(out) <- c("fitnessEffects")
    return(out)
}


rtAndGeneModule <- function(mdeps, gM = NULL) {
    ## To show a table of restrictions when there are modules. Do not use
    ## for anything else. Maybe as intermediate to plotting.
    
    ## Specify restriction table of modules and a mapping of modules to
    ## genes. gM is a named vector; names are modules, values are elements
    ## of each module.

    ## We do nothing important if gM is NULL except checks

    ## If there are modules, the table shows the individual genes.
    checkRT(mdeps)
    ## if(ncol(mdeps) != 5)
    ##     stop("mdeps must be of exactly 5 columns")
    ## if(!identical(colnames(mdeps), c("parent", "child", "s", "sh", "typeDep")))
    ##     stop(paste("Column names of mdeps not of appropriate format. ",
    ##                "Should be parent, child, s, sh, typeDep"))
    if(!is.null(gM)) {
        if(any(is.na(match(mdeps[ , 1], names(gM)))))
            stop("Some values in parent not from a known module")
        if(any(is.na(match(mdeps[ , 2], names(gM)))))
            stop("Some values in child not from a known module")
        if(any(is.na(match(names(gM), c(mdeps[, 1], mdeps[, 2])))))
            stop("Some values in module in neither parent or child")
        
        parent <- gM[mdeps[, 1]]
        child <- gM[mdeps[, 2]]
        df <- data.frame(parent = parent,
                         child = child,
                         s = mdeps$s,
                         sh = mdeps$sh,
                         typeDep = mdeps$typeDep,
                         stringsAsFactors = FALSE)
    } else {
        df <- mdeps
    }
    rownames(df) <- seq.int(nrow(df))
    return(df)
}

## wrap.test.rt <- function(rt, gM = NULL) {
##     ## FIXME add epistasis and orderEffects
##     lrt <- allFitnessEffects(rt, geneToModule = gM)
##     ## wrap_test_rt(lrt$long.rt)
##     wrap_test_rt(lrt$long.rt, lrt$geneModule)
## }


wrap.readFitnessEffects <- function(rt, epi, oe, ni, gm, echo = TRUE) {
    tt <- allFitnessEffects(rt, epi, oe, ni, gm)
    readFitnessEffects(tt, echo = echo)
    ## readFitnessEffects(tt$long.rt,
    ##                    tt$long.epistasis,
    ##                    tt$long.orderEffects,
    ##                    tt$long.geneNoInt,
    ##                    tt$geneModule,
    ##                    tt$gMOneToOne,
    ##                    echo = TRUE)
}

evalGenotype <- function(genotype, fitnessEffects, verbose = FALSE) {
    if(!is.integer(genotype)) {
        gm <- fitnessEffects$geneModule
        genotype <- gm$GeneNumID[match(genotype, gm$Gene)]
    }
    if(any(is.na(genotype)))
        stop("genotype contains NAs or genes not in fitnessEffects")
    if(!length(genotype))
        stop("genotypes must have at least one mutated gene")
    if(any(genotype < 0)) {
        stop(paste("genotypes cannot contain negative values.",
                   "If you see this message, you found a bug."))
    }
    evalRGenotype(genotype, fitnessEffects, verbose)
}

## For multiple genotypes, lapply the matching.
## Nope, I think unneeded
## internal.convert_genotypes <- function(genotypes, gm) {
##     genotypes <- lapply(lg, function(x) gm$GeneNumID[match(x, gm$Gene)])
## }


evalAllGenotypes <- function(fitnessEffects, order = TRUE, max = 256) {

    if(order)
        tot <- function(n) {sum(sapply(seq.int(n),
                                       function(x) choose(n, x) * factorial(x)))}
    else
        tot <- function(n) {2^n}
    nn <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
    tnn <- tot(nn)
    if(tnn > max) {
        m <- paste("There are ", tnn, "genotypes.")
        m <- paste(m, " This is larger than max.")
        m <- paste(m, "Adjust max and rerun if you want")
        stop(m)
    }
    if(order) {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) permutations(n = n, r = x))}
        f2 <- function(names) {
            n <- length(names)
            lapply(seq.int(n), function(x) permutations(n = n, r = x, v = names))
        }
    } else {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) combinations(n = n, r = x))}
        f2 <- function(names) {
            n <- length(names)
            lapply(seq.int(n), function(x) combinations(n = n, r = x, v = names))
        } 
    }
    genotNums <- f1(nn)
    names <- c(fitnessEffects$geneModule$Gene[-1],
               fitnessEffects$long.geneNoInt$Gene)
    genotNames <- f2(names)

    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    genotNums <- list.of.vectors(genotNums)
    genotNames <- unlist(lapply(list.of.vectors(genotNames),
                                function(z) paste(z,
                                                  collapse = if(order){" > "} else {", "})))
    
    allf <- vapply(genotNums,
                   function(x) evalRGenotype(x, fitnessEffects, FALSE),
                   1.1)
    df <- data.frame(Genotype = genotNames, Fitness = allf,
                     stringsAsFactors = FALSE)
    return(df)
}



print.genotToFitness <- function(x) {
    print(x$rt)

}



## FIXME
## - print.genotToFitness?
## - wrap.test.rt and evalGenotype
## - getting epistatis and order in C++
## - testing of examples in C++: several genotypes with evalGenotype
##      - see below for synthetic lethality, etc, examples
##      - add XOR example: it is done via epistasis code.



