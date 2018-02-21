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
    ## Yes, Root: we want no ambiguities.

    ## Actually, that sucks. So we do not require it, but check for
    ## consistency.
    
    if(any(gm == "Root") && (gm[1] != "Root") )
        stop("If Root is in the module table, it must be the first")

    if(any(names(gm) == "Root") && (names(gm)[1] != "Root") )
        stop("If the name Root is in the module table, it must be the first")

    if( (names(gm)[1] == "Root") && (gm[1] != "Root") )
        stop("The name Root is in the module table, but is not of Root")

    if( (gm[1] == "Root") && (names(gm)[1] != "Root") )
        stop("Root is in the module table, but with a different name")

    ## if(gm[1] != "Root")
    ##     stop("First value of a module table must be Root")
    ## if(names(gm)[1] != "Root")
    ##     stop("First name of a module table must be Root")
    if(length(unique(names(gm))) != length(gm))
        stop("Number of unique module names different from length of vector")

    if(gm[1] != "Root")
        gm <- c("Root" = "Root", gm)
    return(gm)
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
    gm <- check.gm(gm)
   
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



allFitnessORMutatorEffects <- function(rT = NULL,
                                       epistasis = NULL,
                                       orderEffects = NULL,
                                       noIntGenes = NULL,
                                       geneToModule = NULL,
                                       drvNames = NULL,
                                       keepInput = TRUE,
                                       ## refFE = NULL,
                                       calledBy = NULL) {
    ## From allFitnessEffects. Generalized so we deal with Fitness
    ## and mutator.
    
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

    ## called appropriately?
    if( !(calledBy %in% c("allFitnessEffects", "allMutatorEffects") ))
        stop("How did you call this function?. Bug.")
    
    if(calledBy == "allMutatorEffects") {
        ## very paranoid check
        if( !is.null(rT) || !is.null(orderEffects) || !is.null(drvNames))
            stop("allMutatorEffects called with forbidden arguments.",
                 "Is this an attempt to subvert the function?")
    }
    
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
        if(inherits(noIntGenes, "character")) {
            wm <- paste("noIntGenes is a character vector.",
                        "This is probably not what you want, and will",
                        "likely result in an error downstream.",
                        "You can get messages like",
                        " 'not compatible with requested type', and others.",
                        "We are stopping.")
            stop(wm)
        }
            
        mg <- max(geneModule[, "GeneNumID"])
        gnum <- seq_along(noIntGenes) + mg
        if(!is.null(names(noIntGenes))) {
            ng <- names(noIntGenes)
            if( grepl(",", ng, fixed = TRUE) || grepl(">", ng, fixed = TRUE)
                || grepl(":", ng, fixed = TRUE))
                stop("The name of some noIntGenes contain a ',' or a '>' or a ':'")
            if(any(ng %in% geneModule[, "Gene"] ))
                stop("A gene in noIntGenes also present in the other terms")
            if(any(duplicated(ng)))
                stop("Duplicated gene names in geneNoInt")
            if(any(is.na(ng)))
                stop("In noIntGenes some genes have names, some don't.",
                     " Name all of them, or name none of them.")
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

    if(calledBy == "allFitnessEffects") {
        if((length(long.rt) + length(long.epistasis) + length(long.orderEffects)) > 1) {
            graphE <- fitnessEffectsToIgraph(rT, epistasis, orderEffects)
        } else {
            graphE <- NULL
        }
    } else {
        graphE <- NULL
    }
    if(!is.null(drvNames)) {
        drv <- unique(getGeneIDNum(geneModule, geneNoInt, drvNames))
        ## drivers should never be in the geneNoInt; Why!!!???
        ## Catch the problem. This is an overkill,
        ## so since we catch the issue, we could leave the geneNoInt. But
        ## that should not be there in this call.
        ## if(any(drvNames %in% geneNoInt$Gene)) {
        ##     stop(paste("At least one gene in drvNames is a geneNoInt gene.",
        ##                "That is not allowed.",
        ##                "If that gene is a driver, pass it as gene in the epistasis",
        ##                "component."))
        ## }
        ## drv <- getGeneIDNum(geneModule, NULL, drvNames)
    } else {
        ## we used to have this default
        ## drv <- geneModule$GeneNumID[-1]
        drv <- vector(mode = "integer", length = 0L)
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
    if(calledBy == "allFitnessEffects") {
        class(out) <- c("fitnessEffects")
    } else if(calledBy == "allMutatorEffects") {
        class(out) <- c("mutatorEffects")
    }
    return(out)
}


allFitnessEffects <- function(rT = NULL,
                              epistasis = NULL,
                              orderEffects = NULL,
                              noIntGenes = NULL,
                              geneToModule = NULL,
                              drvNames = NULL,
                              genotFitness = NULL,
                              keepInput = TRUE) {

    if(!is.null(genotFitness)) {
        if(!is.null(rT) || !is.null(epistasis) ||
           !is.null(orderEffects) || !is.null(noIntGenes) ||
           !is.null(geneToModule)) {
            stop("You have a non-null genotFitness.",
                 " If you pass the complete genotype to fitness mapping",
                 " you cannot pass any of rT, epistasis, orderEffects",
                 " noIntGenes or geneToModule.")
        }
        epistasis <- from_genotype_fitness(genotFitness)
    }
    allFitnessORMutatorEffects(
        rT = rT,
        epistasis = epistasis,
        orderEffects = orderEffects,
        noIntGenes = noIntGenes,
        geneToModule = geneToModule,
        drvNames = drvNames,
        keepInput = keepInput,
        calledBy = "allFitnessEffects")
}


## allFitnessEffects <- function(rT = NULL,
##                               epistasis = NULL,
##                               orderEffects = NULL,
##                               noIntGenes = NULL,
##                               geneToModule = NULL,
##                               drvNames = NULL,
##                               keepInput = TRUE) {
##     ## restrictions: the usual rt

##     ## epistasis: as it says, with the ":"

##     ## orderEffects: the ">"
    
##     ## All of the above can be genes or can be modules (if you pass a
##     ## geneToModule)

##     ## rest: rest of genes, with fitness


##     ## For epistasis and order effects we create the output object but
##     ## missing the numeric ids of genes. With rT we do it in one go, as we
##     ## already know the mapping of genes to numeric ids. We could do the
##     ## same in epistasis and order, but we would be splitting twice
##     ## (whereas for rT extracting the names is very simple).

    
##     rtNames <- NULL
##     epiNames <- NULL
##     orNames <- NULL
##     if(!is.null(rT)) {
##         ## This is really ugly, but to prevent the stringsAsFactors I need it here:
##         rT$parent <- as.character(rT$parent)
##         rT$child <- as.character(rT$child)
##         rT$typeDep <- as.character(rT$typeDep)
##         rtNames <- unique(c(rT$parent, rT$child))
##     }
##     if(!is.null(epistasis)) {
##         long.epistasis <- to.long.epist.order(epistasis, ":")
##         ## epiNames <- unique(unlist(lapply(long.epistasis, function(x) x$ids)))
##         ## deal with the possible negative signs
##         epiNames <- setdiff(unique(
##             unlist(lapply(long.epistasis,
##                           function(x) lapply(x$ids,
##                                              function(z) strsplit(z, "^-"))))),
##                             "")
        
##     } else {
##         long.epistasis <- list()
##     }
##     if(!is.null(orderEffects)) {
##         long.orderEffects <- to.long.epist.order(orderEffects, ">")
##         orNames <- unique(unlist(lapply(long.orderEffects, function(x) x$ids)))
##     } else {
##         long.orderEffects <- list()
##     }
##     allModuleNames <- unique(c(rtNames, epiNames, orNames))
##     if(is.null(geneToModule)) {
##         gMOneToOne <- TRUE
##         geneToModule <- geneModuleNull(allModuleNames)
##     } else {
##         gMOneToOne <- FALSE
##         if(any(is.na(match(setdiff(names(geneToModule), "Root"), allModuleNames))))
##             stop(paste("Some values in geneToModule not present in any of",
##                        " rT, epistasis, or order effects"))
##         if(any(is.na(match(allModuleNames, names(geneToModule)))))
##             stop(paste("Some values in rT, epistasis, ",
##                        "or order effects not in geneToModule"))
##     }
##     geneModule <- gm.to.geneModuleL(geneToModule, one.to.one = gMOneToOne)

##     idm <- unique(geneModule$ModuleNumID)
##     names(idm) <- unique(geneModule$Module)

##     if(!is.null(rT)) {
##         checkRT(rT)
##         long.rt <- to.long.rt(rT, idm)
##     } else {
##         long.rt <- list() ## yes, we want an object of length 0
##     }

##     ## Append the numeric ids to epistasis and order
##     if(!is.null(epistasis)) {
##         long.epistasis <- lapply(long.epistasis,
##                                  function(x)
##                                      addIntID.epist.order(x, idm,
##                                                           sort = TRUE,
##                                                           sign = TRUE))
##     }
##     if(!is.null(orderEffects)) {
##         long.orderEffects <- lapply(long.orderEffects,
##                                     function(x)
##                                         addIntID.epist.order(x, idm,
##                                                              sort = FALSE,
##                                                              sign = FALSE))
##     }
    
##     if(!is.null(noIntGenes)) {
##         mg <- max(geneModule[, "GeneNumID"])
##         gnum <- seq_along(noIntGenes) + mg
##         if(!is.null(names(noIntGenes))) {
##             ng <- names(noIntGenes)
##             if(any(ng %in% geneModule[, "Gene"] ))
##                 stop("A gene in noIntGenes also present in the other terms")
##         } else {
##             ng <- gnum
##         }
##         geneNoInt <- data.frame(Gene = as.character(ng),
##                                 GeneNumID = gnum,
##                                 s = noIntGenes,
##                                 stringsAsFactors = FALSE)
##     } else {
##         geneNoInt <- data.frame()
##     }
##     if( (length(long.rt) + length(long.epistasis) + length(long.orderEffects) +
##              nrow(geneNoInt)) == 0)
##         stop("You have specified nothing!")

##     if((length(long.rt) + length(long.epistasis) + length(long.orderEffects)) > 1) {
##         graphE <- fitnessEffectsToIgraph(rT, epistasis, orderEffects)
##     } else {
##         graphE <- NULL
##     }

##     if(!is.null(drvNames)) {
##         drv <- getGeneIDNum(geneModule, geneNoInt, drvNames)
##     } else {
##         drv <- geneModule$GeneNumID[-1]
##     }
##     if(!keepInput) {
##         rT <- epistasis <- orderEffects <- noIntGenes <- NULL
##     }
##     out <- list(long.rt = long.rt,
##                 long.epistasis = long.epistasis,
##                 long.orderEffects = long.orderEffects,
##                 long.geneNoInt = geneNoInt,
##                 geneModule = geneModule,
##                 gMOneToOne = gMOneToOne,
##                 geneToModule = geneToModule,
##                 graph = graphE,
##                 drv = drv,
##                 rT = rT,
##                 epistasis = epistasis,
##                 orderEffects = orderEffects,
##                 noIntGenes = noIntGenes                
##                 )
##     class(out) <- c("fitnessEffects")
##     return(out)
## }


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



evalGenotypeORMut <- function(genotype,
                              fmEffects,
                              verbose = FALSE,
                              echo = FALSE,
                              model = "",
                              calledBy_= NULL) {
    ## genotype can be a vector of integers, that are the exact same in
    ## the table of fmEffects or a vector of strings, or a vector (a
    ## string) with genes separated by "," or ">"

    if( !(calledBy_ %in% c("evalGenotype", "evalGenotypeMut") ))
        stop("How did you call this function?. Bug.")
    
    if(echo)
        cat(paste("Genotype: ", genotype))
    if(!is.integer(genotype)) {
        if(length(grep(">", genotype))) {
            genotype <- nice.vector.eo(genotype, ">")
        } else if(length(grep(",", genotype))) {
            genotype <- nice.vector.eo(genotype, ",")
        }
        all.g.nums <- c(fmEffects$geneModule$GeneNumID,
                        fmEffects$long.geneNoInt$GeneNumID)
        all.g.names <- c(fmEffects$geneModule$Gene,
                         fmEffects$long.geneNoInt$Gene)
        genotype <- all.g.nums[match(genotype, all.g.names)]
    }
    if(any(is.na(genotype)))
        stop("genotype contains NAs or genes not in fitnessEffects/mutatorEffects")
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
    ff <- evalRGenotype(rG = genotype,
                        rFE = fmEffects,
                        verbose = verbose,
                        prodNeg = prodNeg,
                        calledBy_ = calledBy_)

    
    if(echo) {
        if(calledBy_ == "evalGenotype") {
            if(!prodNeg)
                cat(" Fitness: ", ff, "\n")
            else
                cat(" Death rate: ", ff, "\n")
        } else if(calledBy_ == "evalGenotypeMut") {
            cat(" Mutation rate product :", ff, "\n")
        }
        
    } 
    return(ff)
}

evalGenotype <- function(genotype, fitnessEffects,
                         verbose = FALSE,
                         echo = FALSE,
                         model = "") {
    if(inherits(fitnessEffects, "mutatorEffects"))
         stop("You are trying to get the fitness of a mutator specification. ",
             "You did not pass an object of class fitnessEffects.")

    evalGenotypeORMut(genotype = genotype,
                       fmEffects = fitnessEffects,
                       verbose = verbose,
                       echo = echo,
                       model  = model ,
                       calledBy_= "evalGenotype"
                       )
}


evalGenotypeFitAndMut <- function(genotype,
                                  fitnessEffects,
                                  mutatorEffects,
                                  verbose = FALSE,
                                  echo = FALSE,
                                  model = "") {
    prodNeg <- FALSE
    ## Next is from evalGenotypeAndMut
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

    full2mutator_ <- matchGeneIDs(mutatorEffects,
                                  fitnessEffects)$Reduced
    if(model %in% c("Bozic", "bozic1", "bozic2") )
        prodNeg <- TRUE
    else
        prodNeg <- FALSE
    evalRGenotypeAndMut(genotype, fitnessEffects,
                        mutatorEffects, full2mutator_,
                        verbose = verbose,
                        prodNeg = prodNeg)
}

## evalGenotype <- function(genotype, fitnessEffects,
##                          verbose = FALSE,
##                          echo = FALSE,
##                          model = "") {
##     ## genotype can be a vector of integers, that are the exact same in
##     ## the table of fitnessEffects or a vector of strings, or a vector (a
##     ## string) with genes separated by "," or ">"
    
##     if(echo)
##         cat(paste("Genotype: ", genotype))
##     if(!is.integer(genotype)) {
##         if(length(grep(">", genotype))) {
##             genotype <- nice.vector.eo(genotype, ">")
##         } else if(length(grep(",", genotype))) {
##             genotype <- nice.vector.eo(genotype, ",")
##         }
##         all.g.nums <- c(fitnessEffects$geneModule$GeneNumID,
##                         fitnessEffects$long.geneNoInt$GeneNumID)
##         all.g.names <- c(fitnessEffects$geneModule$Gene,
##                          fitnessEffects$long.geneNoInt$Gene)
##         genotype <- all.g.nums[match(genotype, all.g.names)]
##     }
##     if(any(is.na(genotype)))
##         stop("genotype contains NAs or genes not in fitnessEffects")
##     if(!length(genotype))
##         stop("genotypes must have at least one mutated gene")
##     if(any(genotype < 0)) {
##         stop(paste("genotypes cannot contain negative values.",
##                    "If you see this message, you found a bug."))
##     }
##     if(model %in% c("Bozic", "bozic1", "bozic2") )
##         prodNeg <- TRUE
##     else
##         prodNeg <- FALSE
##     ff <- evalRGenotype(genotype, fitnessEffects, verbose, prodNeg)


##     if(echo) {
##         if(!prodNeg)
##             cat(" Fitness: ", ff, "\n")
##         else
##             cat(" Death rate: ", ff, "\n")
##     } ## else {
##     ##     return(ff)
##     ## }
##     return(ff)
## }

## For multiple genotypes, lapply the matching.
## Nope, I think unneeded
## internal.convert_genotypes <- function(genotypes, gm) {
##     genotypes <- lapply(lg, function(x) gm$GeneNumID[match(x, gm$Gene)])
## }

## I am here: simplify this

evalAllGenotypesORMut <- function(fmEffects,
                                  order = FALSE, max = 256,
                             addwt = FALSE,
                             model = "",
                             calledBy_ = "") {
##                             minimal = FALSE) {
    if( !(calledBy_ %in% c("evalGenotype", "evalGenotypeMut") ))
        stop("How did you call this function?. Bug.")
    
    if( (calledBy_ == "evalGenotype" ) &&
        (!inherits(fmEffects, "fitnessEffects")))
        stop("You are trying to get the fitness of a mutator specification. ",
             "You did not pass an object of class fitnessEffects.")
    if( (calledBy_ == "evalGenotypeMut" ) &&
        (!inherits(fmEffects, "mutatorEffects")))
        stop("You are trying to get the mutator effects of a fitness specification. ",
             "You did not pass an object of class mutatorEffects.")

    ## if(!minimal)
    allg <- generateAllGenotypes(fitnessEffects = fmEffects,
                                 order = order, max = max)
    ## else
        ## allg <- generateAllGenotypes_minimal(fitnessEffects = fmEffects,
        ##                                      max = max)
    ## if(order)
    ##     tot <- function(n) {sum(sapply(seq.int(n),
    ##                                    function(x) choose(n, x) * factorial(x)))}
    ## else
    ##     tot <- function(n) {2^n}
    ## nn <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
    ## tnn <- tot(nn)
    ## if(tnn > max) {
    ##     m <- paste("There are", tnn, "genotypes.")
    ##     m <- paste(m, "This is larger than max.")
    ##     m <- paste(m, "Adjust max and rerun if you want")
    ##     stop(m)
    ## }
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
    ##                      function(z)
    ##                          paste(z,
    ##                                collapse = if(order){" > "} else {", "} )))
    ## This ain't the best, as we repeatedly read and convert
    ## fitnessEffects.  If this were slow, prepare C++ function that takes
    ## vectors of genotypes and uses same fitnessEffects.


    if(model %in% c("Bozic", "bozic1", "bozic2") )
        prodNeg <- TRUE
    else
        prodNeg <- FALSE
    allf <- vapply(allg$genotNums,
                   function(x) evalRGenotype(x, fmEffects, FALSE,
                                             prodNeg,
                                             calledBy_),
                   1.1)
    df <- data.frame(Genotype = allg$genotNames, Fitness = allf,
                     stringsAsFactors = FALSE)
    ## Why am I doing this? I am not computing the genotype.  I test the
    ## evaluation of the empty genotype in the tests. But its evaluation
    ## generates warnings that are not needed. And by decree (even in the
    ## C++) this thing has a fitness of 1, a mutator effect of 1 since
    ## there are no terms.
    
    if(addwt)
        df <- rbind(data.frame(Genotype = "WT", Fitness = 1,
                               stringsAsFactors = FALSE), df)
    if(calledBy_ == "evalGenotype") { 
        if(prodNeg)
            colnames(df)[match("Fitness", colnames(df))] <- "Death_rate"
        class(df) <- c("evalAllGenotypes", class(df))
    } else if (calledBy_ == "evalGenotypeMut") {
        colnames(df)[match("Fitness", colnames(df))] <- "MutatorFactor"
        class(df) <- c("evalAllGenotypesMut", class(df))
    }
    return(df)
}

evalAllGenotypes <- function(fitnessEffects, order = FALSE, max = 256,
                             addwt = FALSE,
                             model = "") {
    evalAllGenotypesORMut(
        fmEffects = fitnessEffects,
        order = order,
        max = max,
        addwt = addwt,
        model = model,
        calledBy_= "evalGenotype"
    )
}

generateAllGenotypes <- function(fitnessEffects, order = TRUE, max = 256) {
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
    ## With mutator, the ids of genes need not go from 1:n
    vid <- allNamedGenes(fitnessEffects)$GeneNumID
    if(order) {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vid))
        }
    } else {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vid))}
        
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
    return(list(genotNums = genotNums,
                genotNames = genotNames
                ))
}

evalAllGenotypesFitAndMut <- function(fitnessEffects, mutatorEffects,
                                   order = FALSE, max = 256,
                                   addwt = FALSE,
                                   model = "" ){
##                                   minimal = FALSE) {
    ## if(!minimal)
    allg <- generateAllGenotypes(fitnessEffects = fitnessEffects,
                                 order = order, max = max)
    ## else
        ## allg <- generateAllGenotypes_minimal(fitnessEffects = fitnessEffects,
        ##                                      max = max)
    
    if(model %in% c("Bozic", "bozic1", "bozic2") ) {
        prodNeg <- TRUE
    } else {
        prodNeg <- FALSE
    }
    
    full2mutator_ <- matchGeneIDs(mutatorEffects,
                                  fitnessEffects)$Reduced
    allf <- t(vapply(allg$genotNums,
                   function(x) evalRGenotypeAndMut(x,
                                                   rFE = fitnessEffects,
                                                   muEF= mutatorEffects,
                                                   full2mutator_ = full2mutator_,
                                                   verbose = FALSE,
                                                   prodNeg = prodNeg),
                   c(1.1, 2.2)))

    df <- data.frame(Genotype = allg$genotNames, Fitness = allf[, 1],
                     MutatorFactor = allf[, 2],
                     stringsAsFactors = FALSE)
    if(addwt)
        df <- rbind(data.frame(Genotype = "WT", Fitness = 1,
                               MutatorFactor = 1,
                               stringsAsFactors = FALSE), df)
    if(prodNeg)
        colnames(df)[match("Fitness", colnames(df))] <- "Death_rate"
    class(df) <- c("evalAllGenotypesFitAndMut", class(df))
    return(df)
}






## evalAllGenotypes <- function(fitnessEffects, order = TRUE, max = 256,
##                              addwt = FALSE,
##                              model = "") {
    
##     if(order)
##         tot <- function(n) {sum(sapply(seq.int(n),
##                                        function(x) choose(n, x) * factorial(x)))}
##     else
##         tot <- function(n) {2^n}
##     nn <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
##     tnn <- tot(nn)
##     if(tnn > max) {
##         m <- paste("There are", tnn, "genotypes.")
##         m <- paste(m, "This is larger than max.")
##         m <- paste(m, "Adjust max and rerun if you want")
##         stop(m)
##     }
##     if(order) {
##         f1 <- function(n) {
##             lapply(seq.int(n), function(x) permutations(n = n, r = x))
##         }
##     } else {
##         f1 <- function(n) {
##             lapply(seq.int(n), function(x) combinations(n = n, r = x))}
        
##     }
##     genotNums <- f1(nn)
##     list.of.vectors <- function(y) {
##         ## there's got to be a simpler way
##         lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
##                function(m) m[[1]])
##     }
##     genotNums <- list.of.vectors(genotNums)
##     names <- c(fitnessEffects$geneModule$Gene[-1],
##                fitnessEffects$long.geneNoInt$Gene)

##     genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]),
##                          function(z)
##                              paste(z,
##                                    collapse = if(order){" > "} else {", "} )))
##     ## This ain't the best, as we repeatedly read and convert
##     ## fitnessEffects.  If this were slow, prepare C++ function that takes
##     ## vectors of genotypes and uses same fitnessEffects.
##     if(model %in% c("Bozic", "bozic1", "bozic2") )
##         prodNeg <- TRUE
##     else
##         prodNeg <- FALSE
##     allf <- vapply(genotNums,
##                    function(x) evalRGenotype(x, fitnessEffects, FALSE, prodNeg),
##                    1.1)
##     df <- data.frame(Genotype = genotNames, Fitness = allf,
##                      stringsAsFactors = FALSE)
##     if(addwt)
##         df <- rbind(data.frame(Genotype = "WT", Fitness = 1,
##                                stringsAsFactors = FALSE), df)
##     if(prodNeg)
##         colnames(df)[match("Fitness", colnames(df))] <- "Death_rate"
##     return(df)
## }

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
    ## for special case of simple epi effects
    if(nrow(df) == 0) return(NULL)
    
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
                                  K,
                                  detectionDrivers,
                                  onlyCancer,
                                  errorHitWallTime,
                                  max.num.tries,
                                  errorHitMaxTries,
                                  minDetectDrvCloneSz,
                                  extraTime,
                                  keepPhylog,
                                  detectionProb,
                                  MMUEF = NULL ## avoid partial matching, and set default
                                  ) {
    if(!inherits(rFE, "fitnessEffects"))
        stop(paste("rFE must be an object of class fitnessEffects",
                   "as created, for instance, with function",
                   "allFitnessEffects"))

    if(countGenesFe(rFE) < 2) {
        stop("There must be at least two genes (loci) in the fitness effects.",
             "If you only care about a degenerate case with just one,",
             "you can enter a second gene (locus)",
             "with fitness effect of zero.")
    }
    
    if( (length(mu) == 1) && !(is.null(names(mu)))) {
        stop("A length 1 mutation, but named. ",
             "This is ambiguous. ",
             "If you want per-gene mutation rates, each gene",
             "must have its entry in the mu vector. ",
             "(And regardless of per-gene mutation rates ",
             " there must be at least two gene/loci).")
    }

    namedGenes <- allNamedGenes(rFE)

    if( length(mu) > 1) {
        if(is.null(names(mu)))
            stop("When using per-gene mutation rates the ",
                 "mu vector must be named ",
                 "(and if you have noIntGenes, those must have names).")
        if(length(mu) != countGenesFe(rFE))
            stop("When using per-gene mutation rates, ",
                 "there must be the same number of genes in the ",
                 "mu vector and the fitness effects object.")
        if(!identical(sort(namedGenes$Gene),
                      sort(names(mu))))
            stop("When using per-gene mutation rates, ",
                 "names of genes must match those in the ",
                 "restriction table.")
        mu <- mu[order(match(names(mu), namedGenes$Gene))]
        ## Hyperparanoid check. Should never, ever, happen.
        if(!identical(names(mu), namedGenes$Gene))
            stop("Names of re-ordered mu do not match names of genes")
        minmu <- 1e-40
        if(any(mu <= minmu))
            stop(paste("At least one per-gene mutation rate is negative",
                       "or less than", minmu,". Remember that the per-base",
                       "mutation rate in the human genome is about 1e-11 to 1e-9."))
    }
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
       if(length(initMutant) >= countGenesFe(rFE)) {
        stop("For initMutant you passed as many, or more genes, mutated ",
             "than the number of genes in the genotype (fitness effects).")
    }
       
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

    if(!is.null(MMUEF)) {
        if(!inherits(MMUEF, "mutatorEffects"))
            stop("muEF must be a mutatorEffects object")
        full2mutator_ <- matchGeneIDs(MMUEF, rFE)
        if(any(is.na(full2mutator_$Full)))
            stop("Genes in mutatorEffects not present in fitnessEffects")
        if(any(is.na(full2mutator_)))
            stop("full2mutator failed for unknown reasons")
        full2mutator_ <- full2mutator_$Reduced
    } else {
        full2mutator_ <- vector(mode = "numeric", length = 0)
        ## muEF <- emptyFitnessEffects()
    }

    dpr <- detectionProbCheckParse(detectionProb, initSize, verbosity)
    ## if( !is.null(cPDetect) && (sum(!is.null(p2), !is.null(n2)) >= 1 ))
    ##     stop("Specify only cPDetect xor both of p2 and n2")
    ## if( (is.null(p2) + is.null(n2)) == 1 )
    ##     stop("If you pass one of n2 or p2, you must also pass the other. ",
    ##          "Otherwise, we would not know what to do.")
    ## stopifnot(PDBaseline >= 0)
    ## stopifnot(n2 > PDBaseline)
    ## stopifnot(p2 < 1)
    ## stopifnot(p2 > 0)
    ## if( is.null(cPDetect) ) cPDetect <- -9
    ## if( is.null(p2)) p2 <- 9
    ## if( is.null(n2)) n2 <- -9

    ## call <- match.call()
    if(typeFitness == 'discreteModel'){
        return(c(
            discreteModel(rFE = rFE,
                    mu = mu,
                    popIni = initSize,
                    tMax = finalTime,
                    seed = seed,
                    sampleEvery = sampleEvery,
                    keepEvery = keepEvery,
                    verbosity = verbosity),
            Drivers = list(rFE$drv), ## but when doing pops, these will be repeated
            geneNames = list(names(getNamesID(rFE)))
        ))
    }else{
        return(c(
            nr_BNB_Algo5(rFE = rFE,
                     mu_ = mu,
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
                     K = K,
                     detectionDrivers = detectionDrivers,
                     onlyCancer = onlyCancer,
                     errorHitWallTime = errorHitWallTime,
                     maxNumTries = max.num.tries,
                     errorHitMaxTries = errorHitMaxTries,
                     minDetectDrvCloneSz = minDetectDrvCloneSz,
                     extraTime = extraTime,
                     keepPhylog = keepPhylog,
                     MMUEF = MMUEF,
                     full2mutator_ = full2mutator_,
                     n2 = dpr["n2"],
                     p2 = dpr["p2"],
                     PDBaseline = dpr["PDBaseline"],
                     cPDetect_i= dpr["cPDetect"],
                     checkSizePEvery = dpr["checkSizePEvery"]),
            Drivers = list(rFE$drv), ## but when doing pops, these will be repeated
            geneNames = list(names(getNamesID(rFE)))
        ))

    }

}


countGenesFe <- function(fe) {
    ## recall geneModule has Root always
    nrow(fe$geneModule) + nrow(fe$long.geneNoInt) - 1
}

allNamedGenes <- function(fe){
    ## Returns a data frame with genes and their names and verifies all
    ## genes have names.

    ## Root should always be first, but just in case
    ## avoid assuming it

    ## This does is not a good idea as it assumes the user did not use
    ## "keepInput = FALSE".
    ## lni <- length(fe$noIntGenes)
    ## ## FIXME:test
    ## if((lni > 0) &&
    ##        (is.null(names(fe$noIntGenes))))
    ##         stop("When using per-gene mutation rates the ",
    ##              "no interaction genes must be named ",
    ##              "(i.e., the noIntGenes vector must have names).")
    
    v1 <- fe$geneModule[, c("Gene", "GeneNumID")]
    if(nrow(fe$long.geneNoInt)) {
        v1 <- rbind(v1,
                    fe$long.geneNoInt[, c("Gene", "GeneNumID")])
    }
    v1 <- v1[-which(v1[, "Gene"] == "Root"), ]
    rownames(v1) <- NULL
    return(v1)
}


get.gene.counts <- function(x) {
                                        # , timeSample = "last",
                                        # typeSample = "whole") {
    ## From get.mut.vector. Used for now just for testing
    timeSample <- "last"
    the.time <- get.the.time.for.sample(x, timeSample, NULL)
    if(the.time < 0) {
        counts <- rep(0, length(x$geneNames))
        names(counts) <- x$geneNames
        freq <- rep(NA, length(x$geneNames))
        names(freq) <- x$geneNames
        return(list(counts = counts,
                    freq = freq,
                    popSize = 0))
    } 
    pop <- x$pops.by.time[the.time, -1]
    if(all(pop == 0)) {
        stop("You found a bug: this should never happen")
    }
    ## if(typeSample %in% c("wholeTumor", "whole")) {
    popSize <- x$PerSampleStats[the.time, 1]
    counts <- as.vector(tcrossprod(pop, x$Genotypes))
    names(counts) <- x$geneNames    
    return(list(counts = counts,
                freq = counts/popSize,
                popSize = popSize))
    ## return( (tcrossprod(pop,
    ##                     x$Genotypes)/popSize) )
    ## } else if (typeSample %in%  c("singleCell", "single")) {

    ##       return(x$Genotypes[, sample(seq_along(pop), 1, prob = pop)])
    ##   } else {
    ##         stop("Unknown typeSample option")
    ##     }
}



geneCounts <- function(x) {
    if(inherits(x, "oncosimulpop")) {
        return( do.call(rbind,
                        lapply(x,
                               function(z)
                               get.gene.counts(z)$counts))
               )
    } else {
        return( get.gene.counts(x)$counts)
    }
}





## geneNumIDReset <- function(x, ref){
##     ## Set GeneNumID of a fitnessEffect object, x, using ref as the
##     ## reference fitness effect object.
##     ## Check also if all in x are in ref.

##     gg <- allNamedGenes(ref)
##     gnid <- gg$GeneNumID
##     names(gnid) <- gg$Gene
##     ## FIXME: this later and conditional on what is in thee
##     gnid <- c("Root" = 0, gnid)
    
##     if(!all(x$geneModule$Gene %in% names(gnid) ))
##         stop("Some genes not in reference fitnessEffects (rebasing geneModule)")
##     x$geneModule$GeneNumID <- gnid[geneModule$Gene]

##     ## and then go over all the lists in the x object. 
    
##     if(nrow(x$long.geneNoInt)) {
##         ## now, mapping for the noInt if this is mutator
##         if(!all(x$long.geneNoInt$Gene %in% names(gnid) ))
##             stop("Some genes not in reference fitnessEffects (rebasing geneNoInt)")
##         x$long.geneNoInt$GeneNumID <- gnid[long.geneNoInt$Gene]
##     }
## }


matchGeneIDs <- function(x, refFE) {
    n1 <- allNamedGenes(x)
    n2 <- allNamedGenes(refFE)
    colnames(n1)[2] <- "Reduced"
    colnames(n2)[2] <- "Full"
    ## Non standard evaluation thing and a note being generated in check.
    ## See also http://stackoverflow.com/questions/33299845/when-writing-functions-for-an-r-package-should-i-avoid-non-standard-evaluation
    ## But does not work well with replace. So use the NULL trick
    Reduced <- NULL
    dplyr::full_join(n2, n1, by = "Gene") %>%
        mutate(Reduced = replace(Reduced, is.na(Reduced), -9))
}

    
detectionProbCheckParse <- function(x, initSize, verbosity) {
    default_p2 <- 0.1
    default_n2 <- 2 * initSize
    default_PDBaseline <- 1.2 * initSize
    default_checkSizePEvery <- 20
    ## No default cPDetect. That is done from p2 and n2 in C++.
    
    if((length(x) == 1) && (is.na(x))) {
        y <- vector()
        y["cPDetect"] <- -9
        y["p2"] <- 9
        y["n2"] <- -9
        y["PDBaseline"] <- -9
        y["checkSizePEvery"] <- Inf
        return(y)
    } else if((length(x) == 1) && (x == "default")) {
        ## Populate with defaults
        y <- vector()
        y["p2"] <- default_p2
        y["n2"] <- default_n2
        y["PDBaseline"] <- default_PDBaseline
        y["checkSizePEvery"] <- default_checkSizePEvery
        x <- y
    }

    if(length(x) >= 1) {
        if( !all(names(x) %in% c("cPDetect", "p2", "n2", "PDBaseline", "checkSizePEvery")))
            stop("Names of some components of detectionProb are not recognized.",
                 " Check for possible typos.")
    }
   
    ## This ain't conditional. If not given, always check
    if( !is.na(x["cPDetect"]) && (sum(!is.na(x["p2"]), !is.na(x["n2"])) >= 1 ))
        stop("Specify only cPDetect xor both of p2 and n2")
    if( (is.na(x["p2"]) + is.na(x["n2"])) == 1 )
        stop("If you pass one of n2 or p2, you must also pass the other. ",
             "Otherwise, we would not know what to do.")

    if(is.na(x["PDBaseline"])) {
        x["PDBaseline"] <- default_PDBaseline
        if(verbosity > -1)
            message("\n  PDBaseline set to default value of ", default_PDBaseline, "\n")
        }
    if(is.na(x["checkSizePEvery"])) {
        x["checkSizePEvery"] <- default_checkSizePEvery
        if(verbosity > -1)
            message("\n  checkSizePEvery set to default value of ",
                default_checkSizePEvery, "\n")
        }

    ## If we get here, we checked consistency of whether cPDetect or p2/n2.
    ## Fill up with defaults the missing values

    if(is.na(x["cPDetect"])) {
        if(is.na(x["p2"])) {
            if(!is.na(x["n2"])) stop("Eh? no p2 but n2? Bug")
            x["n2"] <- default_n2
            x["p2"] <- default_p2
            if(verbosity > -1)
                message("\n  n2, p2 set to default values of ",
                    default_n2, ", ", default_p2, "\n")
        }
    }
    
    if(x["checkSizePEvery"] <= 0)
        stop("checkSizePEvery <= 0")
    if(x["PDBaseline"] < 0)
        stop("PDBaseline < 0")
    
    if(!is.na(x["n2"])) {
        if(x["n2"] <= x["PDBaseline"])
            stop("n2 <= PDBaseline")
        if(x["p2"] >= 1) stop("p2 >= 1")
        if(x["p2"] <= 0) stop("p2 <= 0")
        x["cPDetect"] <- -9
    } else { ## only if x["cPDetect"] is not NA
        if(is.na(x["cPDetect"])) stop("eh? you found a bug")## paranoia
        x["n2"] <- -9
        x["p2"] <- -9
        if(verbosity > -1)
            message("\n Using user-specified cPDetect as n2, p2 not given.\n")
    }
    return(x)
}


## emptyFitnessEffects <- function() {
##     list(long.rt = list(),
##          long.epistasis = list(),
##          long.orderEffects = list(),
##          long.geneNoInt = list(),
##          geneModule = list(),
##          gMOneToOne = TRUE,
##          geneToModule = list(),
##          graph = NULL,
##          drv = vector(mode = "integer", length = 0)
##          )
## }



### Later, for all the effects, we will do some kind of dplyr match?

### a, b, c, in fitness, only a, c in mut.
### fitness table for a,b,c
### each row name transformed removing b (so leaving only present)
### each row transformed matched to row in mut table.

## t1 <- data.frame(v1 = c("a,b", "a,c", "b"), v2 = c("b", "c", "b"), v3 = 1:3, stringsAsFactors = FALSE)
## t2 <- data.frame(v2 = c("b", "c"), v4 = c(11, 12), stringsAsFactors = FALSE)
## full_join(t1, t2, by = "v2")



## FIXME

## new exit code
## check:
## baseline < n2
## p2 < 1

## baseline defaults to init size + .1
## use c < -9, n2 < -9, p2 < -9 for no values
## and check one of c or p2 and n2 are valid if using exit
