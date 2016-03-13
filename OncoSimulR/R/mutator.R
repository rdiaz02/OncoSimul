## Any gene not specified has a mutator effect of 1.
## So noIntGenes constructed from this noIntGenes plus the list
## of genes of the full fitness effect.
## Or not really? Better not. No need to fill things up.

## But .. what about the genome?

## Any mutator gene must have been in fitness, even if it has fitness
## effect of 0. Yes, because otherwise the mapping of gene names to
## numerical ids becomes a pita.



allMutatorEffects <- function(fe,
                              epistasis = NULL,
                              noIntGenes = NULL,
                              geneToModule = NULL,
                              keepInput = TRUE) {
    ## From allFitnesseffects, removing stuff and setting rT and
    ## orderEffects to NULL
    
    rT <- orderEffects <- NULL
    
    
    rtNames <- NULL
    epiNames <- NULL
    orNames <- NULL
    ## if(!is.null(rT)) {
    ##     ## This is really ugly, but to prevent the stringsAsFactors I need it here:
    ##     rT$parent <- as.character(rT$parent)
    ##     rT$child <- as.character(rT$child)
    ##     rT$typeDep <- as.character(rT$typeDep)
    ##     rtNames <- unique(c(rT$parent, rT$child))
    ## }
    
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
    ## if(!is.null(orderEffects)) {
    ##     long.orderEffects <- to.long.epist.order(orderEffects, ">")
    ##     orNames <- unique(unlist(lapply(long.orderEffects, function(x) x$ids)))
    ## } else {
    ##     long.orderEffects <- list()
    ## }
    
    ## allModuleNames <- unique(c(rtNames, epiNames, orNames))

    allModuleNames <- unique(epiNames)

    if(is.null(geneToModule)) {
        gMOneToOne <- TRUE
        geneToModule <- geneModuleNull(allModuleNames)
    } else {
        ## add Root here directly.
        if(geneToModule[1] != "Root")
            geneToModule <- c("Root" = "Root", geneToModule)
        gMOneToOne <- FALSE
        if(any(is.na(match(setdiff(names(geneToModule), "Root"), allModuleNames))))
            stop(paste("Some values in geneToModule not present in any of",
                       " rT, epistasis, or order effects"))
        if(any(is.na(match(allModuleNames, names(geneToModule)))))
            stop(paste("Some values in rT, epistasis, ",
                       "or order effects not in geneToModule"))
    }
    geneModule <- gm.to.geneModuleL(geneToModule, one.to.one = gMOneToOne)

    ## Set the mapping gene to GeneNumID here, using same mapping as in
    ## fitness effects.
    gg <- allNamedGenes(fea)
    gnid <- gg$GeneNumID
    names(gnid) <- gg$Gene
    gnid <- c("Root" = 0, gnid)
    geneModule$GeneNumID <- gnid[geneModule$Gene]
    
    ## FIXME: Remove Root? Probably not: expected in C++?
    idm <- unique(geneModule$ModuleNumID)
    names(idm) <- unique(geneModule$Module)

    ## if(!is.null(rT)) {
    ##     checkRT(rT)
    ##     long.rt <- to.long.rt(rT, idm)
    ## } else {
    ##     long.rt <- list() ## yes, we want an object of length 0
    ## }

    ## Append the numeric ids to epistasis and order
    if(!is.null(epistasis)) {
        long.epistasis <- lapply(long.epistasis,
                                 function(x)
                                     addIntID.epist.order(x, idm,
                                                          sort = TRUE,
                                                          sign = TRUE))
    }
    ## if(!is.null(orderEffects)) {
    ##     long.orderEffects <- lapply(long.orderEffects,
    ##                                 function(x)
    ##                                     addIntID.epist.order(x, idm,
    ##                                                          sort = FALSE,
    ##                                                          sign = FALSE))
    ## }
    
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
    if( (length(long.epistasis) + 
             nrow(geneNoInt)) == 0)
        stop("You have specified nothing for mutator effects!")
    
    ## if( (length(long.rt) + length(long.epistasis) + length(long.orderEffects) +
    ##          nrow(geneNoInt)) == 0)
    ##     stop("You have specified nothing!")

    ## if((length(long.rt) + length(long.epistasis) + length(long.orderEffects)) > 1) {
    ##     graphE <- fitnessEffectsToIgraph(rT, epistasis, orderEffects)
    ## } else {
    ##     graphE <- NULL
    ## }

    ## if(!is.null(drvNames)) {
    ##     drv <- getGeneIDNum(geneModule, geneNoInt, drvNames)
    ## } else {
    ##     drv <- geneModule$GeneNumID[-1]
    ## }

    if(!keepInput) {
        epistasis <- noIntGenes <- NULL
    }
    
    ## if(!keepInput) {
    ##     rT <- epistasis <- orderEffects <- noIntGenes <- NULL
    ## }
    out <- list(
        long.rt = NULL,
        long.epistasis = long.epistasis,
        long.orderEffects = NULL,
        long.geneNoInt = geneNoInt,
        geneModule = geneModule,
        gMOneToOne = gMOneToOne,
        geneToModule = geneToModule,
        ## graph = NULL,
        ##                drv = drv,
        ## rT = rT,
        epistasis = epistasis,
        ## orderEffects = orderEffects,
        noIntGenes = noIntGenes                
    )
    class(out) <- c("mutatorEffects")
    return(out)
}



evalGenotypeMut <- function(genotype, fitnessEffects,
                         verbose = FALSE,
                         echo = FALSE) {
    ## From evalGenotype

    
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
    ff <- evalRGenotypeMut(genotype, fitnessEffects, verbose)


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


