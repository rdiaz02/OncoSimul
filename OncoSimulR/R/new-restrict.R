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





## Say which are drivers: populate the drv vector.
## // In R, the user says which are the drivers. If does not say anuthing,
## // the default (NULL) then drivers are all in poset, epist, restrict. The
## // user can pass a vector with the names of the genes (not modules). Allow
## // also for empty, so this is faster if not needed. And check that if we
## // use a stopping rule based on drivers that drv vectors is not empty.


## - Say that a user can use a "0" as a gene name, but that is BAD idea.
## - Modules and order effects can be kind of funny?

## Gene names can contain no spaces, commas, or ">"


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

## nice.vector.eo <- function(z, sep) {
##     ## with epistasis, maybe we want sorted?
##     setdiff(unlist(lapply(strsplit(z, " "),
##                                     function(u) strsplit(u, sep))), "")
## }

nice.vector.eo <- function(z, sep, rm.sign = FALSE) {
    ## with epistasis, maybe we want sorted?
    if(! rm.sign)
        return(setdiff(unlist(lapply(strsplit(z, " "),
                              function(u) strsplit(u, sep))), ""))
    else ## remove the " ", the -, and the sep
        return(setdiff(unlist(lapply(strsplit(z, " "), function(u)
            strsplit(unlist(strsplit(u, "-")), sep))), ""))
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
        ## Works, but does not give the most ordered module names. But
        ## keeps implicit order given by user.
        idm <- seq.int(length(gm) - 1)
        idm <- c("Root" = 0L, idm)
        names(idm) <- names(gm)
        ## This sorts by the names but is not optimal either
        ## idm1 <- seq.int(length(gm) - 1)
        ## idm <- c(0L, idm1)
        ## names(idm) <- c("Root", sort(setdiff(names(gm), "Root")))
        geneMod$ModuleNumID <- idm[geneMod[, "Module"]]
    }
    if(length(unique(geneMod$Gene)) != nrow(geneMod)) {
        stop("Are there identical gene names in different modules?")
    }
    ## I think this is unreacheable now. Caught earlier.
    if(length(unique(geneMod$GeneNumID)) != nrow(geneMod)) {
        stop("Are there identical gene names in different modules?")
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
                       "AND" = "monotone",
                       "monotone" = "monotone",
                       "CMPN" = "monotone",
                       "OR" = "semimonotone",
                       "SM" = "semimonotone",
                       "semimonotone" = "semimonotone",
                       "DMPN" = "semimonotone",
                       "XOR" = "xmpn",
                       "xmpn" = "xmpn",
                       "XMPN" = "xmpn",
                       "--"   = "--",
                       "-" = "--")
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
    typeDep <- lookupTypeDep[unique(x$typeDep)]
    if(any(is.na(typeDep)))
        stop("typeDep value incorrect. Check spelling.")
    return(list(
        child = unique(x$child),
        s = unique(x$s),
        sh = unique(x$sh),
        typeDep = typeDep,
        parents = unlist(x$parent)))

}

to.long.rt <- function(rt, idm) {
    ## We now do this inconditionally, so that we do not need to use the
    ## "stringsAsFactors = FALSE". This is now done before
    ## if(is.numeric(rt$parent))
    ##     rt$parent <- as.character(rt$parent)
    ## if(is.numeric(rt$child))
    ##     rt$child <- as.character(rt$child)
   
    
    if(!("Root" %in% rt$parent))
        stop("Root must be one parent node")

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
            ## I think this can no longer be reached ever. Caught before.
            stop(paste("An ID is NA:",
                       "Is a gene part of two different modules?",
                       "(That includes being by itself and part",
                       "of a module.)"))
            
        }
        ## These checks could be somewhere else instead of here.
        if(length(unique(z$parentsNumID)) != length(z$parentsNumID))
            stop("No node can have the same parent multiple times")
        if(length(unique(z$parents)) != length(z$parents))
            stop("No node can have the same parent multiple times")
        if(length(z$child) > 1)
            stop("Child nodes can have one, and only one, member")
        
        ## sort the parents, always.
        o <- order(z$parentsNumID)
        z$parentsNumID <- z$parentsNumID[o]
        z$parents <- as.character(z$parents[o])
        return(z)
    }
    long.rt <- lapply(long.rt, function(x) addIntID(x, idm = idm))
   
    ## if(verbosity >= 4) {
    ##     message(paste("Number of drivers: ",
    ##                   length(unique(geneModule[, "Gene"]))))
    ##     message(paste("Number of modules: ",
    ##                   length(unique(geneModule[, "Module"]))))
    ## }
    return(long.rt)
    ## return(list(long.rt = long.rt, geneModule = geneModule))
}


epist.order.element <- function(x, y, sep, rm.sign = FALSE) {
    list(ids = nice.vector.eo(x, sep = sep, rm.sign = rm.sign), s = y)
}

oe.to.df <- function(x) {
    ma <- matrix(ncol = 2, nrow = 1 + length(x) - 2)
    if(length(x) == 2) {
        ma[1, 1] <- x[1]
        ma[1, 2] <- x[2]
    } else {
        ma[, 1] <- x[-length(x)]
        ma[, 2] <- x[-1]
    }
    return(data.frame(ma, stringsAsFactors = FALSE))
}


epist.order.to.pairs.modules <- function(x, sep, rm.sign = TRUE) {
    ## We discard, do not even accept, the coefficient
    tmp <- epist.order.element(x, -99, sep = sep, rm.sign = rm.sign)$ids
    if(length(tmp) > 1) {
        ## if a single gene, as when we specify all genotypes, we do not
        ## want this
        if(sep == ":")
            return(data.frame(combinations(n = length(tmp), r = 2, v = tmp),
                              stringsAsFactors = FALSE))
        else if(sep == ">") {
            return(oe.to.df(tmp))
        }
    }
}

to.pairs.modules <- function(x, sep, rm.sign = TRUE) {
    df <- data.frame(rbindlist(
        lapply(names(x),
               function(z) epist.order.to.pairs.modules(z, sep))))
    if(nrow(df) == 0L) { ## if only single genes in epist, we get nothing here.
        return(df)
    } else {
        colnames(df) <- c("parent", "child")
        if(sep == ":")
            df$typeDep <- "epistasis"
        else if(sep == ">")
            df$typeDep <- "orderEffect"
        return(unique(df))
    }
}



to.long.epist.order <- function(epor, sep, rm.sign = FALSE) {
    ## just vectors for now
    long <- Map(function(x, y) epist.order.element(x, y, sep, rm.sign),
                names(epor), epor)
    ## if(is.vector(epor))
    ##     long <- Map(function(x, y) epist.order.element(x, y, sep, rm.sign),
    ##                    names(epor), epor)
    ## else if(is.data.frame(epor)) 
    ##     long <- Map(function(x, y) epist.order.element(x, y, sep, rm.sign),
    ##                 as.character(epor$ids),
    ##                 epor$s)
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
    if(length(unique(z$NumID)) != length(z$NumID))
        stop("No node can have the same id multiple times")
    if(length(unique(z$ids)) != length(z$ids))
        stop("No node can have the same id multiple times")
    return(z)
}


checkRT <- function(mdeps) {
    if(ncol(mdeps) != 5)
        stop("mdeps must be of exactly 5 columns")
    if(!identical(colnames(mdeps), c("parent", "child", "s", "sh", "typeDep")))
        stop(paste("Column names of mdeps not of appropriate format. ",
                   "Should be parent, child, s, sh, typeDep"))
}


getNamesID <- function(fp) {
    ## Return a lookup table for names based on numeric IDs
    idname <- c(fp$geneModule$GeneNumID,  fp$long.geneNoInt$GeneNumID)
    names(idname) <- c(fp$geneModule$Gene, fp$long.geneNoInt$Gene)
    return(idname[-1]) ## remove Root!!
}


getGeneIDNum <- function(geneModule, geneNoInt, drv, sort = TRUE) {
    ## It returns the genes, as NumID, in the given vector with names drv
    ## initMutant uses this, for simplicity, without sorting, but noInt
    ## are always sorted

    ## Also used for the drivers with sort = TRUE

    ## Yes, we must do it twice because we do not know before hand which
    ## is which. This makes sure no NA. Period.
    if(any(is.na( match(drv, c(geneModule$Gene, geneNoInt$Gene))))) {
        stop(paste("For driver or initMutant you have passed genes",
                   "not in the fitness table."))
    }
    
    indicesM <- as.vector(na.omit(match( drv, geneModule$Gene)))
    indicesI <- as.vector(na.omit(sort(match( drv, geneNoInt$Gene))))
    if(sort) {
        indicesM <- sort(indicesM)
    }
    return(c(
        geneModule$GeneNumID[indicesM],
        geneNoInt$GeneNumID[indicesI])
    )
}

allFitnessEffects <- function(rT = NULL,
                              epistasis = NULL,
                              orderEffects = NULL,
                              noIntGenes = NULL,
                              geneToModule = NULL,
                              drvNames = NULL,
                              keepInput = TRUE) {
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
    if(!is.null(rT)) {
        ## This is really ugly, but to prevent the stringsAsFactors I need it here:
        rT$parent <- as.character(rT$parent)
        rT$child <- as.character(rT$child)
        rT$typeDep <- as.character(rT$typeDep)
        rtNames <- unique(c(rT$parent, rT$child))
    }
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

    if((length(long.rt) + length(long.epistasis) + length(long.orderEffects)) > 1) {
        graphE <- fitnessEffectsToIgraph(rT, epistasis, orderEffects)
    } else {
        graphE <- NULL
    }

    if(!is.null(drvNames)) {
        drv <- getGeneIDNum(geneModule, geneNoInt, drvNames)
    } else {
        drv <- geneModule$GeneNumID[-1]
    }
    if(!keepInput) {
        rT <- epistasis <- orderEffects <- noIntGenes <- NULL
    }
    out <- list(long.rt = long.rt,
                long.epistasis = long.epistasis,
                long.orderEffects = long.orderEffects,
                long.geneNoInt = geneNoInt,
                geneModule = geneModule,
                gMOneToOne = gMOneToOne,
                geneToModule = geneToModule,
                graph = graphE,
                drv = drv,
                rT = rT,
                epistasis = epistasis,
                orderEffects = orderEffects,
                noIntGenes = noIntGenes                
                )
    class(out) <- c("fitnessEffects")
    return(out)
}


## No longer used
## rtAndGeneModule <- function(mdeps, gM = NULL) {
##     ## To show a table of restrictions when there are modules. Do not use
##     ## for anything else. Maybe as intermediate to plotting.
    
##     ## Specify restriction table of modules and a mapping of modules to
##     ## genes. gM is a named vector; names are modules, values are elements
##     ## of each module.

##     ## We do nothing important if gM is NULL except checks

##     ## If there are modules, the table shows the individual genes.
##     checkRT(mdeps)
##     ## if(ncol(mdeps) != 5)
##     ##     stop("mdeps must be of exactly 5 columns")
##     ## if(!identical(colnames(mdeps), c("parent", "child", "s", "sh", "typeDep")))
##     ##     stop(paste("Column names of mdeps not of appropriate format. ",
##     ##                "Should be parent, child, s, sh, typeDep"))
##     if(!is.null(gM)) {
##         if(any(is.na(match(mdeps[ , 1], names(gM)))))
##             stop("Some values in parent not from a known module")
##         if(any(is.na(match(mdeps[ , 2], names(gM)))))
##             stop("Some values in child not from a known module")
##         if(any(is.na(match(names(gM), c(mdeps[, 1], mdeps[, 2])))))
##             stop("Some values in module in neither parent or child")
        
##         parent <- gM[mdeps[, 1]]
##         child <- gM[mdeps[, 2]]
##         df <- data.frame(parent = parent,
##                          child = child,
##                          s = mdeps$s,
##                          sh = mdeps$sh,
##                          typeDep = mdeps$typeDep,
##                          stringsAsFactors = FALSE)
##     } else {
##         df <- mdeps
##     }
##     rownames(df) <- seq.int(nrow(df))
##     return(df)
## }

## wrap.test.rt <- function(rt, gM = NULL) {
##     ## FIXME add epistasis and orderEffects
##     lrt <- allFitnessEffects(rt, geneToModule = gM)
##     ## wrap_test_rt(lrt$long.rt)
##     wrap_test_rt(lrt$long.rt, lrt$geneModule)
## }

## No longer used
## wrap.readFitnessEffects <- function(rt, epi, oe, ni, gm, echo = TRUE) {
##     tt <- allFitnessEffects(rt, epi, oe, ni, gm)
##     readFitnessEffects(tt, echo = echo)
##     ## readFitnessEffects(tt$long.rt,
##     ##                    tt$long.epistasis,
##     ##                    tt$long.orderEffects,
##     ##                    tt$long.geneNoInt,
##     ##                    tt$geneModule,
##     ##                    tt$gMOneToOne,
##     ##                    echo = TRUE)
## }

evalGenotype <- function(genotype, fitnessEffects,
                         verbose = FALSE,
                         echo = FALSE,
                         model = "") {
    ## genotype can be a vector of integers, that are the exact same in
    ## the table of fitnessEffects or a vector of strings, or a vector (a
    ## string) with genes separated by "," or ">"
    
    if(echo)
        cat(paste("Genotype: ", genotype))
    if(!is.integer(genotype)) {
        if(length(grep(">", genotype))) {
            genotype <- nice.vector.eo(genotype, ">")
        } else if(length(grep(",", genotype))) {
            genotype <- nice.vector.eo(genotype, ",")
        }
        all.g.nums <- c(fitnessEffects$geneModule$GeneNumID,
                        fitnessEffects$long.geneNoInt$GeneNumID)
        all.g.names <- c(fitnessEffects$geneModule$Gene,
                         fitnessEffects$long.geneNoInt$Gene)
        genotype <- all.g.nums[match(genotype, all.g.names)]
    }
    if(any(is.na(genotype)))
        stop("genotype contains NAs or genes not in fitnessEffects")
    if(!length(genotype))
        stop("genotypes must have at least one mutated gene")
    if(any(genotype < 0)) {
        stop(paste("genotypes cannot contain negative values.",
                   "If you see this message, you found a bug."))
    }
    if(model %in% c("Bozic", "bozic1", "bozic2") )
        prodNeg <- TRUE
    else
        prodNeg <- FALSE
    ff <- evalRGenotype(genotype, fitnessEffects, verbose, prodNeg)


    if(echo) {
        if(!prodNeg)
            cat(" Fitness: ", ff, "\n")
        else
            cat(" Death rate: ", ff, "\n")
    } ## else {
    ##     return(ff)
    ## }
    return(ff)
}

## For multiple genotypes, lapply the matching.
## Nope, I think unneeded
## internal.convert_genotypes <- function(genotypes, gm) {
##     genotypes <- lapply(lg, function(x) gm$GeneNumID[match(x, gm$Gene)])
## }


evalAllGenotypes <- function(fitnessEffects, order = TRUE, max = 256,
                             addwt = FALSE,
                             model = "") {
    
    if(order)
        tot <- function(n) {sum(sapply(seq.int(n),
                                       function(x) choose(n, x) * factorial(x)))}
    else
        tot <- function(n) {2^n}
    nn <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
    tnn <- tot(nn)
    if(tnn > max) {
        m <- paste("There are", tnn, "genotypes.")
        m <- paste(m, "This is larger than max.")
        m <- paste(m, "Adjust max and rerun if you want")
        stop(m)
    }
    if(order) {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) permutations(n = n, r = x))
        }
    } else {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) combinations(n = n, r = x))}
        
    }
    genotNums <- f1(nn)
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    genotNums <- list.of.vectors(genotNums)
    names <- c(fitnessEffects$geneModule$Gene[-1],
               fitnessEffects$long.geneNoInt$Gene)

    genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]),
                         function(z)
                             paste(z,
                                   collapse = if(order){" > "} else {", "} )))
    ## This ain't the best, as we repeatedly read and convert
    ## fitnessEffects.  If this were slow, prepare C++ function that takes
    ## vectors of genotypes and uses same fitnessEffects.
    if(model %in% c("Bozic", "bozic1", "bozic2") )
        prodNeg <- TRUE
    else
        prodNeg <- FALSE
    allf <- vapply(genotNums,
                   function(x) evalRGenotype(x, fitnessEffects, FALSE, prodNeg),
                   1.1)
    df <- data.frame(Genotype = genotNames, Fitness = allf,
                     stringsAsFactors = FALSE)
    if(addwt)
        df <- rbind(data.frame(Genotype = "WT", Fitness = 1,
                               stringsAsFactors = FALSE), df)
    if(prodNeg)
        colnames(df)[match("Fitness", colnames(df))] <- "Death_rate"
    return(df)
}


fitnessEffectsToIgraph <- function(rT, epistasis, orderEffects) {

    df0 <- df1 <- df2 <- data.frame()
             
    if(!is.null(rT)) {
        df0 <- rT[, c("parent", "child", "typeDep")]
    }
    if(!is.null(epistasis)) {
        df1 <- to.pairs.modules(epistasis, ":")
    }
    if(!is.null(orderEffects)) {
        df2 <- to.pairs.modules(orderEffects, ">")
    }
    df <- rbind(df0, df1, df2)
    g1 <- graph.data.frame(df)
    E(g1)$color <- "black"
    E(g1)$color[E(g1)$typeDep == "SM"] <- "blue"
    E(g1)$color[E(g1)$typeDep == "XMPN"] <- "red"
    E(g1)$color[E(g1)$typeDep == "epistasis"] <- "orange"
    E(g1)$color[E(g1)$typeDep == "orderEffect"] <- "orange"
    E(g1)$lty <- 1
    E(g1)$lty[E(g1)$typeDep == "epistasis"] <- 2
    E(g1)$lty[E(g1)$typeDep == "orderEffect"] <- 3
    E(g1)$arrow.mode <- 2
    E(g1)$arrow.mode[E(g1)$typeDep == "epistasis"] <- 0
    E(g1)$arrow.mode[E(g1)$typeDep == "orderEffect"] <- 2
    return(g1)
}


plot.fitnessEffects <- function(x, type = "graphNEL",
                                layout = NULL,
                                expandModules = FALSE,
                                autofit = FALSE,
                                scale_char = ifelse(type == "graphNEL", 1/10, 5),
                                return_g = FALSE,
                                lwdf = 1, ## multiplier for lwd for graphNEL
                                ...) {
    ## some other layouts I find OK
    ## layout.circle
    ## layout.reingold.tilford if really a tree
    ## o.w. it will use the default
    g <- x$graph
    
    if(type == "igraph") {
        if(expandModules && (!x$gMOneToOne)) {
            ## vlabels <- fe$geneToModule[vertex.attributes(g)$name]
            vlabels <- x$geneToModule[V(g)$name]
            V(g)$label <- vlabels
        }
        if(autofit) {
            plot(0, type = "n", axes = FALSE, ann = FALSE)
            ## ideas from http://stackoverflow.com/questions/14472079/match-vertex-size-to-label-size-in-igraph
            ## vsize <- (strwidth(V(g)$label) + strwidth("oo")) * 200
            ## but this is a kludge.
            vsize <- (nchar(V(g)$label) + 1) * scale_char
            plot.igraph(g, vertex.size = vsize, vertex.shape = "rectangle",
                        layout = layout)
        } else {
            plot.igraph(g, layout = layout)
        }
        if(return_g) return(g)
    }
    else if (type == "graphNEL") {
        g1 <- igraph.to.graphNEL(g)
        c1 <- unlist(lapply(edgeData(g1), function(x) x$color))
        names(c1) <- sub("|", "~", names(c1), fixed = TRUE)
        s1 <- unlist(lapply(edgeData(g1), function(x) x$lty))
        names(s1) <- names(c1)
        a1 <- unlist(lapply(edgeData(g1), function(x) max(x$arrow.mode - 1, 0)))
        names(a1) <- names(c1)
        lwd <- s1
        lwd[lwd == 2] <- 2 ## o.w. too thin
        lwd[lwd == 3] <- 2 ## o.w. too thin
        lwd <- lwd * lwdf
        nAttrs <- list()
        if(expandModules && (!x$gMOneToOne)) {
            nnodes <- x$geneToModule[nodes(g1)]
            names(nnodes) <- nodes(g1)
            nAttrs$label <- nnodes
        }
        if(autofit) {
            nAttrs$width <- (nchar(nAttrs$label) + 1) * scale_char
            names(nAttrs$width) <- names(nAttrs$label)
            plot(g1, edgeAttrs = list(arrowsize = a1, lty = s1, lwd = lwd,
                         color = c1), attrs=list(node=list(shape = "rectangle")),
                 nodeAttrs = nAttrs)
            
        } else {
            plot(g1, edgeAttrs = list(arrowsize = a1, lty = s1, lwd = lwd,
                         color = c1),
                 nodeAttrs = nAttrs)
        }
        if(return_g) return(g1)
    } else {
        stop("plot type not recognized")
    }
}

## plot.fitnessEffects <- plotFitnessEffects

## FIXME: in the help: say we cannot properly show 3- or higher order
## interactions.



nr_oncoSimul.internal <- function(rFE, 
                                  birth, 
                                  death,
                                  mu,
                                  initSize,
                                  sampleEvery,
                                  detectionSize,
                                  finalTime,
                                  initSize_species,
                                  initSize_iter,
                                  seed,
                                  verbosity,
                                  speciesFS,
                                  ratioForce,
                                  typeFitness,
                                  max.memory,
                                  mutationPropGrowth,
                                  initMutant,
                                  max.wall.time,
                                  keepEvery,
                                  alpha,
                                  K,
                                  detectionDrivers,
                                  onlyCancer,
                                  errorHitWallTime,
                                  max.num.tries,
                                  errorHitMaxTries,
                                  minDetectDrvCloneSz,
                                  extraTime,
                                  keepPhylog) {
    if(!inherits(rFE, "fitnessEffects"))
        stop(paste("rFE must be an object of class fitnessEffects",
                   "as created, for instance, with function",
                   "allFitnessEffects"))
    if(!is.null(initMutant)) {
       if(length(grep(">", initMutant))) {
            initMutant <- nice.vector.eo(initMutant, ">")
        } else if(length(grep(",", initMutant))) {
            initMutant <- nice.vector.eo(initMutant, ",")
        }
        initMutant <- getGeneIDNum(rFE$geneModule,
                             rFE$long.geneNoInt,
                                   initMutant,
                                   FALSE)
    } else {
        initMutant <- vector(mode = "integer")
    }
    ## these are never user options
    ## if(initSize_species < 10) {
    ##     warning("initSize_species too small?")
    ## }
    ## if(initSize_iter < 100) {
    ##     warning("initSize_iter too small?")
    ## }

    if(typeFitness %in% c("bozic1", "bozic2")) {
        thesh <- unlist(lapply(rFE$long.rt, function(x) x$sh))
        ## thes <- unlist(lapply(rFE$long.rt, function(x) x$s))
        thes <- unlist(c(lapply(rFE$long.rt, function(x) x$s),
                         lapply(rFE$long.epistasis, function(x) x$s),
                         lapply(rFE$long.orderEffects, function(x) x$s),
                         rFE$long.geneNoInt$s))
        
        if(any(thes > 1 )) {
            m <- paste("You are using a Bozic model with",
                       "the new restriction specification, and you have",
                       "at least one s > 1.",
                       "But that is the same as setting that to 1:",
                       "we obviously cannot allow negative death rates,",
                       "nor problems derived from multiplying odd or even",
                       "numbers of negative numbers.",
                       "You can set those > 1 to 1, but then you will",
                       "eventually get the simulations to abort when the",
                       "death rate becomes 0.")
            stop(m)
        }
        if(any( (thesh <= -1) & (thesh > -Inf))) {
            m <- paste("You are using a Bozic model with",
                       "the new restriction specification, and you have",
                       "at least one sh <= -1. Maybe you mean -Inf?")
            warning(m)
        }
        if(any(thes == 1 )) {
            m <- paste("You are using a Bozic model with",
                       "the new restriction specification, and you have",
                       "at least one s of 1. If that gene is mutated, ",
                       "this will lead to a death rate of 0",
                       "and the simulations will abort when you get a non finite",
                       "value.")
            warning(m)
        }
    }
    ## call <- match.call()
    return(c(
        nr_BNB_Algo5(rFE = rFE,
                     mu = mu,
                 death = death,
                 initSize = initSize,
                 sampleEvery = sampleEvery,
                 detectionSize = detectionSize,
                 finalTime = finalTime,
                 initSp = initSize_species,
                 initIt = initSize_iter,
                 seed = seed,
                 verbosity = verbosity,
                 speciesFS = speciesFS,
                 ratioForce = ratioForce,
                 typeFitness_ = typeFitness,
                 maxram = max.memory,
                 mutationPropGrowth = as.integer(mutationPropGrowth),
                 initMutant_ = initMutant, 
                 maxWallTime = max.wall.time,
                 keepEvery = keepEvery,
                 alpha = alpha,
                 K = K,
                 detectionDrivers = detectionDrivers,
                 onlyCancer = onlyCancer,
                 errorHitWallTime = errorHitWallTime,
                 maxNumTries = max.num.tries,
                 errorHitMaxTries = errorHitMaxTries,
                 minDetectDrvCloneSz = minDetectDrvCloneSz,
                     extraTime = extraTime,
                     keepPhylog),
        Drivers = list(rFE$drv), ## but when doing pops, these will be repeated
        geneNames = list(names(getNamesID(rFE)))
    ))
}
















