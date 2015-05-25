## - Say that a user can use a "0" as a gene name, but that is BAD idea.
## - Modules and order effects can be kind of funny?
library(data.table)
library(Rcpp)
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

gm.to.geneModuleL <- function(gm) {
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
    idm <- seq.int(length(gm) - 1)
    idm <- c("Root" = 0L, idm)
    names(idm) <- names(gm)
    geneMod$ModuleNumID <- idm[geneMod[, "Module"]]
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
        z$parents <- z$parents[o]
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
    geneModule <- gm.to.geneModuleL(geneToModule)

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
        geneNoInt <- data.frame(Gene = ng,
                                GeneNumID = gnum,
                                s = noIntGenes)
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
    ## class(out) <- c("fitnessEffects")
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

epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))
oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)



epineg <- c("-a:d" = 0.2, "a:d" = 0.3, "d:c" = 0.3)
epineg2 <- c("-a:d" = 0.2, "b:c" = 0.3)
epineg3 <- c("a:-d" = 0.2, "b:c" = 0.3)


allFitnessEffects(epistasis = epineg, geneToModule = NULL)

allFitnessEffects(epistasis = epineg2, geneToModule = NULL)

allFitnessEffects(epistasis = epineg3, geneToModule = NULL)

gme <- c("Root" = "Root", "a" = "1, 2", "d" = 3, "c" = 4)

gme2 <- c("Root" = "Root", "a" = "1, 2", "d" = 3, "b" = "5, 6", "c" = 4)

allFitnessEffects(epistasis = epineg, geneToModule = gme)
allFitnessEffects(epistasis = epineg2, geneToModule = gme2)

allFitnessEffects(m0, epistasis = epineg, geneToModule = gM3)
allFitnessEffects(m0, epistasis = epineg2, geneToModule = gM3)



oa <- allFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2), gM3)

oa2 <- allFitnessEffects(m0, epistm1,
                        oeffects1, runif(1000), gM3)


benchmark(wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                                  gM3, echo = FALSE),
          replications = 100)

benchmark(allFitnessEffects(m0, epistm1,
                            oeffects1, c(0.1, 0.1, 0.2), gM3),
          replications = 100)

benchmark(readFitnessEffects(oa2, echo = FALSE),
          replications = 10000)
benchmark(readFitnessEffects(oa, echo = FALSE),
          replications = 10000)



microbenchmark(readFitnessEffects(oa, echo = FALSE), times = 1000)

microbenchmark(allFitnessEffects(m0, epistm1,
                            oeffects1, c(0.1, 0.1, 0.2), gM3),
          times = 100)


evalGenotype <- function(genotype, fitnessEffects, verbose = FALSE) {
    if(!is.integer(genotype)) {
        gm <- fitnessEffects$geneModule
        genotype <- gm$GeneNumID[match(genotype, gm$Gene)]
    }
    evalRGenotype(genotype, fitnessEffects, verbose)
}


evalGenotype(c("d8", "2", "6"), oa, verbose = TRUE)


## For multiple genotypes, lapply the matching.
internal.convert_genotypes <- function(genotypes, gm) {
    genotypes <- lapply(lg, function(x) gm$GeneNumID[match(x, gm$Gene)])
}


## evalGenotype <- function(genotype, fitnessEffects, FALSE) {
##     eval_Genotype(genotype, fitnessEffects,  genotype)
## }



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



### examples here

m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)




gM2 <- c("Root" = "Root", "a" = "1, 2", "b2" = "3, 4, 5", "b" = "8",
         "c" = "7")

to.long.rt(m0, gm.to.geneModuleL(gM))

to.long.rt(m0, gm.to.geneModuleL(gM2))

m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)


gM <- c("Root" = "Root", "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")
gm.to.geneModuleL(gM)


allFitnessEffects(m0)

allFitnessEffects(m0, noIntGenes = c(0.1, 0, 0.2))

allFitnessEffects(m0, noIntGenes = c("u" = 0.1, "v" = 0, "mm" = 0.2))


rtAndGeneModule(m0, gM)
rtAndGeneModule(m0)


epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))
oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)

allFitnessEffects(m0, epistasis = epistm1,
                  orderEffects = oeffects1,
                  noIntGenes = c(0.1, 0, 0.2))

wrap.readFitnessEffects(NULL,
                        NULL,
                        NULL,
                        c(0.1, 0.1),
                        NULL)

wrap.readFitnessEffects(m0,
                        NULL,
                        NULL,
                        c(0.1, 0.1),
                        NULL)

wrap.readFitnessEffects(m0,
                        NULL,
                        NULL,
                        NULL,
                        NULL)

wrap.readFitnessEffects(NULL,
                        NULL,
                        NULL,
                        NULL,
                        NULL)


wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        NULL)

wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        gM2)

gM3 <- c("Root" = "Root", "d" = "d9, d8",
         "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")

m00 <- data.frame(parent = c("Root", "c"),
                 child  = c("c", "b"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)


m000 <- data.frame(parent = c("Root", "d", "c", "b"),
                 child  = c("c", "a", "a", "d"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

## good to check geneModule table is OK
gM4 <- c("Root" = "Root", "d" = "d9, d8, a, z9",
         "a" = "z700, u43, 78", "b" = "2, 3, 4, 5", "c" = "1, b, 6")

wrap.readFitnessEffects(m000, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        gM4)

## a complex one with originally disordered parents and children and weird
## dep on 0 and others.
m4 <- data.frame(parent = c(rep("Root", 7), "d", "c", "b", "g", "h", "d", "c", "c", "e"),
                 child  = c(letters[2:8], rep("a", 5), rep("b", 2), rep("g", 2)),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

wrap.readFitnessEffects(m4, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2), NULL)


wrap.readFitnessEffects(m4, epineg,
                        oeffects1, c(0.1, 0.1, 0.2), NULL)



benchmark(wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                                  gM3, echo = FALSE),
          replications = 100)
## FIXME make sure to test with 0 size elements: rT, epist, order





## We do something somewhat silly: we accept as input a set of
## restrictions and then find modules, etc. But it is probably better to
## simply take restrictions and modules separately, which minimizes the
## checking. And simplifies the code. But then, it is already there, and
## is more flexible now. But epistasis and order only as modules.







## FIXME add formal tests. Verify all output carefully.

## FIXME a converter from old posets to the new format.  do this INSIDE
## the code.  If you pass a single two-column poset and s and sh, the new
## format is created.

## FIXME I think we must enforce that 0 always show up. And it is the
## root? Or not, but enforce a policy


## FIXME: for sanity, and when dealing with the rT, etc, and to avoid the stingsAsFactors do
## mm[, "parent"] <- as.character(mm[, "parent"])
## and ditto for child

## synthetic lethality:
## specified as a row with colons: A:B
## change synbol; a box? plain text? getDefaultAttrs()
## plot(g1, nodeAttrs = list(shape = c(a = "box", b = "plaintext")))
## note how simple is to specify the shape of some node.
## see the "Hot to plot a graph using Rgraphviz"

## if using igraph, more shapes.
## http://stackoverflow.com/questions/7429162/how-to-use-igraph-vertex-shape-functionality
## How to get a tree?
## http://stackoverflow.com/questions/18270370/plot-tree-with-graph-tree-function-from-igraph

## but then ... maybe igraph is really the way to go for simuls? Or I
## should use both.


## How to specify order differences: e.g., A first, B second different
## from B first?

## Use ">" in restriction table

## In graphs, use a double-headed arrow, if possible with inverted arrow
## head, and a different colour (say, gray)



library(data.table)
rt2 <- data.frame(parent = c(
                      0, 0, 0,
                      1,
                      "2, 3",
                      4,
                      "5, 6, 7",
                      "5, 6, 7",
                      4,
                      4
                      ),
           child = c(
               1,
               7,
               3,
               2,
               4,
               5,
               8,
               9,
               9,
               6),
                  s = 0.1,
                  sh = 0.05,
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

rt3 <- data.frame(parent = c(
                      0, 0, 0,
                      1,
                      "2, 3",
                      4,
                      "5, 6, 7",
                      "5, 6, 7",
                      4,
                      4,
                      10,
                      10,
                      10,
                      10,
                      10,
                      10,
                      10,
                      10,
                      11,
                      12,
                      "13, 14, 15",
                      "13, 14, 15",                      
                      "16, 17, 18",
                      2,
                      4,
                      3
                      ),
           child = c(
               1,
               7,
               3,
               2,
               4,
               5,
               8,
               9,
               9,
               6,
               11,
               12,
               13,
               14,
               15,
               16,
               17,
               18,
               16,
               19,
               18,
               19,
               19,
               19,
               19,
               10
               ),
                  s = 0.1,
                  sh = 0.05,
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

## rt4 and rt5 are for checking. The values in sh would have no effect here.
rt4 <- data.frame(parent = c(
                      0, 0
                      ),
           child = c(
               1,
               2),
                  s = 0.1,
                  sh = c(0.05, -Inf),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

rt5 <- data.frame(parent = c(
                      0, 0
                      ),
           child = c(
               1,
               2),
                  s = c(-Inf, Inf),
                  sh = c(Inf, -Inf),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)


rt6 <- data.frame(parent = c(
                      0, 0
                      ),
           child = c(
               1,
               2),
                  s = c(0.1, 0.3),
                  sh = c(99, 99),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

rt7 <- data.frame(parent = c(
                      0, 0, 1, 2
                      ),
           child = c(
               1,
               2,
               3,
               4),
                  s = c(0.1, 0.2, 0.3, 0.4),
                  sh = c(99, 99, -0.05, -Inf),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)



rt8 <- data.frame(
    parent = c(
        0,
        0,
        "1,2",
        "1,2",
        5
        ),
    child = c(
        "1, 2",
        5,
        3,
        4,
        4),
    s = c(0.12,  0.5, 0.3, 0.4, 0.4),
    sh = c(99,  99, -0.03, -0.04, -0.04),
    typeDep = "MN",
    stringsAsFactors = FALSE)

rt8.sm <- data.frame(
    parent = c(
        0,
        0,
        "1,2",
        "1,2",
        5
        ),
    child = c(
        "1, 2",
        5,
        3,
        4,
        4),
    s = c(0.12,  0.5, 0.3, 0.4, 0.4),
    sh = c(99,  99, -0.03, -0.04, -0.04),
    typeDep = c(rep("MN",3), "SM", "SM"),
    stringsAsFactors = FALSE)

rt8.b <- data.frame(
    parent = c(
        0,
        0,
        0,
        "1,2",
        "1,2",
        5
        ),
    child = c(
        1,
        2,
        5,
        3,
        4,
        4),
    s = c(0.1, 0.2, 0.5, 0.3, 0.4, 0.4),
    sh = c(99, 99, 99, -0.03, -0.04, -0.04),
    typeDep = "SM",
    stringsAsFactors = FALSE)



rt9 <- data.frame(
    parent = c(
        0,
        0,
        "1,2",
        "1,2"
        ),
    child = c(
        1,
        2,
        3,
        4),
    s = c(0.1, 0.1, 0.1, 0.1),
    sh = c(-1, -2, -1, -2),
    typeDep = "MN",
    stringsAsFactors = FALSE)


rt11 <- data.frame(
    parent = c(
        0,
        0,
        0,
        1,
        2,
        "3,4",
        "3, 4",
        7
        ),
    child = c(
        1,
        2,
        7,
        "3, 4",
        "3,4",
        5,
        6,
        6),
    s = c(0.1, 0.2, 0.7, 0.34, 0.34, 0.5, 0.6, 0.6),
    sh = c(-1, -2, -7, -34, -34, -5, -6, -6),
    typeDep = c(rep("MN", 3), "SM", "SM", rep("MN", 3)),
    stringsAsFactors = FALSE)


rt11s <- data.frame(
    parent = c(
        0,
        0,
        0,
        "myc",
        "ras",
        "PT,alpa",
        "PT, alpa",
        "u1"
        ),
    child = c(
        "myc",
        "ras",
        "u1",
        "PT, alpa",
        "PT,alpa",
        "AC",
        "BG",
        "BG"),
    s = c(0.1, 0.2, 0.7, 0.34, 0.34, 0.5, 0.6, 0.6),
    sh = c(-1, -2, -7, -34, -34, -5, -6, -6),
    typeDep = c(rep("MN", 3), "SM", "SM", rep("MN", 3)),
    stringsAsFactors = FALSE)



rt12 <- data.frame(
    parent = c(
        0,
        0,
        0,
        "M",
        "B",
        "3,4",
        "3, 4",
        "C"
        ),
    child = c(
        "M",
        "B",
        "C",
        "3, 4",
        "3,4",
        "D",
        "E",
        "E"),
    s = c(0.1, 0.2, 0.7, 0.34, 0.34, 0.5, 0.6, 0.6),
    sh = c(-1, -2, -7, -34, -34, -5, -6, -6),
    typeDep = "MN",
    stringsAsFactors = FALSE)

## FIXME
## do I really want "as.integer"
## If I don't, i can use arbitrary things

nice.string <- function(z) {
    ## Spaces are ALWAYS ignored. Do not do silly things like having
    ## spaces in names

    ## why unique? and why sort? Because we want to recognize that, say, a
    ## set of parents as 1,2 is the same as 2,1. That will not happen with
    ## epistasis and order effects.
    paste(setdiff(sort(unique(unlist(lapply(strsplit(z, " "),
                                            function(u) strsplit(u, ","))))),
                  ""),
          collapse = ", ")
}


gtm <- function(x) {
    data.frame(cbind(unlist(strsplit(x, ", ")), x))
}



gene.to.module <- function(rt) {
    ##    all.modules <- unique(unlist(lapply(c(rt$parent, rt$child), nice.string)))
    all.modules <- unique(unlist(c(rt$parent, rt$child)))
    geneMod <- as.data.frame(rbindlist(lapply(all.modules, gtm)))
    colnames(geneMod) <- c("Gene", "Module")
    geneMod$Gene <- as.character(geneMod$Gene)
    geneMod$Module <- as.character(geneMod$Module)
    geneMod <- geneMod[order(geneMod$Gene), ]
    geneMod$GeneNumID <- 0:(nrow(geneMod) - 1)
    geneMod
}



## FIXME: make sure mutations within modules are ordered!!
## This next add to R code.
## FIXME: remember to pass num drivers!!
to.long.rt.original <- function(rt, verbosity = 0) {
    if(is.numeric(rt$parent))
        rt$parent <- as.character(rt$parent)
    if(!("0" %in% rt$parent))
        stop("0 must be one parent node")
    if(is.numeric(rt$child))
        rt$child <- as.character(rt$child)
    rt$parent <- unlist(lapply(rt$parent, nice.string))
    rt$child <- unlist(lapply(rt$child, nice.string))
   
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

    geneModule <- gene.to.module(srt)
    idm <- seq.int(length(names(long.rt)))
    names(idm) <- names(long.rt)
    idm <- c("0" = 0L, idm)
    geneModule$ModuleNumID <- idm[geneModule[, "Module"]]

    ## add integer IDs
    addIntID <- function(z, idm) {
        z$childID <- idm[z$child]
        z$parentsID <- idm[z$parents]
        if( any(is.na(z$parentsID)) ||
           any(is.na(z$childID)) ) {
            stop(paste("An ID is NA:",
                       "Is a gene part of two different modules?",
                       "(That includes being by itself and part",
                       "of a module.)"))
            
        }
        return(z)
    }
    long.rt <- lapply(long.rt, function(x) addIntID(x, idm = idm))
   
    if(verbosity >= 4) {
        message(paste("Number of drivers: ",
                      length(unique(geneModule[, "Gene"]))))
        message(paste("Number of modules: ",
                      length(unique(geneModule[, "Module"]))))
    }
    return(list(long.rt = long.rt, geneModule = geneModule))
}



## rt.to.cpp <- function(rt) {
##     lrt <- to.long.rt(rt)
##     rTable_to_Poset0(lrt)
## }

wrap.test.rt <- function(rt) {
    lrt <- to.long.rt(rt)
    ## wrap_test_rt(lrt$long.rt)
    wrap_test_rt(lrt$long.rt, lrt$geneModule)
}
## was called wrap.test.checkRestrictions






wrap.test.rt(rt12)
wrap.test.rt(rt6)
wrap.test.rt(rt7)
wrap.test.rt(rt8)
wrap.test.rt(rt11)
wrap.test.rt(rt11s)

## test the Inf
wrap.test.rt(rt4)
wrap.test.rt(rt5)



## These are not proper posets
wrap.test.rt(rt9)
wrap.test.rt(rt8.b)
wrap.test.rt(rt2) ## rt2 is not a proper poset
wrap.test.rt(rt3) ## rt3 is not a proper poset
## same gene as module and as isolated.



rt.nr <- data.frame(parent = c(1, 2),
                    child = c(
                        3,
                        4),
                  s = 0.1,
                  sh = 0.05,
                  typeDep = "MN",
                  stringsAsFactors = FALSE)
wrap.test.rt(rt.nr) ## fails, as it should



## FIXME show the semantically equivalent formulation without modules
## for example, this would be equivalent to rt9
rt9ok <- data.frame(
    parent = c(
        0,
        0,
        1,
        2,
        1,
        2
        ),
    child = c(
        1,
        2,
        3,
        3,
        4,
        4),
    s = c(0.1),
    sh = c(-1, -2, -1, -1, -2, -2),
    typeDep = c("MN", "MN", "SM", "SM", "SM", "SM"),
    stringsAsFactors = FALSE)

## FIXME have another way of specification, where we specify dependencies
## between modules AND module membership
## Anything that accepts an rt, has
## function(rt, geneModule = NULL)
## if geneModule is not null, it is the mapping.



## FIXME: if modules, verify they are in the union of rt, epist and order



m0 <- data.frame(parent = c(0, "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)
gM <- c("0" = 0, "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")

rtAndGeneModule(m0, gM)
rtAndGeneModule(m0)

wrap.test.rt(rtAndGeneModule(m0, gM))
wrap.test.rt(rtAndGeneModule(m0))



m1 <- data.frame(parent = c(0, "a", "b", 0),
                 child  = c("a", "b", "c", "d"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))

oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)
oeffects1b <- data.frame(ids = c("d>a", "c > d"),
                         s = c(0.4, -0.3))








to.long.epist.order(epistm1, ":")
to.long.epist.order(epistm1b, ":")

to.long.epist.order(oeffects1, ">")
to.long.epist.order(oeffects1b, ">")




## FIXME user specifies the rt, the epist and the order, and optionally
## the modules.




## FIXME printing of tables: if modules then both in terms of modules and genes.
## FIXME: plotting: similar, giving options of modules or genes.







## yes, at least one element needs to be a character to force all to be chars
m1 <- data.frame(parent = c(0, 1, "2"),
                 child  = c(1, 2, "3"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)
gM1 <- c("0" = 0, "1" = "1,2", "2" = "3,4,5", "3" = "6")

rtAndGeneModule(m1, gM1)
wrap.test.rt(rtAndGeneModule(m1, gM1))

m2 <- data.frame(parent = c("0", "A", "B"),
                 child  = c("A", "B", "C"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)
gM2 <- c("0" = 0, "A" = "pten,myc", "B" = "p53", "C" = "BRCA1, BRCA2")
rtAndGeneModule(m2, gM2)
wrap.test.rt(rtAndGeneModule(m2, gM2))


## FIXME: create an "evalAllGenotypes" that produces a table of the
## fitness of all possible genotypes?

evalGenotype(rt7, c(1L))

evalGenotype(rt7, c(1L))
evalGenotype(rt7, c(2L))

evalGenotype(rt7, c(3L))
evalGenotype(rt7, c(4L))

evalGenotype(rt7, c(1L, 2L))
evalGenotype(rt7, c(1L, 3L))
evalGenotype(rt7, c(1L, 4L))
evalGenotype(rt7, c(2L, 3L))
evalGenotype(rt7, c(2L, 4L))

evalGenotype(rt7, c(1L, 2L, 3L))
evalGenotype(rt7, c(1L, 2L, 4L))
evalGenotype(rt7, c(1L, 3L, 4L))

evalGenotype(rt7, c(2L, 3L, 4L))
evalGenotype(rt7, c(1L, 2L, 3L, 4L))

evalGenotype(rt8, c(1L, 5L))
evalGenotype(rt8, c(1L, 2L, 5L))

evalGenotype(rt8, c(1L, 3L))
evalGenotype(rt8, c(2L, 3L))
evalGenotype(rt8, c(1L, 2L, 3L))
evalGenotype(rt8, c(5L, 3L))
evalGenotype(rt8, c(4L, 3L))

evalGenotype(rt8, c(5L, 4L))
evalGenotype(rt8, c(1L, 4L))
evalGenotype(rt8, c(2L, 4L))
evalGenotype(rt8, c(1L, 2L, 4L))

evalGenotype(rt8, c(1L, 5L, 4L))
evalGenotype(rt8, c(2L, 5L, 4L))

evalGenotype(rt8, c(1L, 2L, 5L, 4L))


evalGenotype(rt8.sm, c(5L, 4L))
evalGenotype(rt8.sm, c(1L, 4L))
evalGenotype(rt8.sm, c(2L, 4L))
evalGenotype(rt8.sm, c(1L, 2L, 4L))


identical(evalGenotype(rt8, c(1L, 5L, 4L)),
          evalGenotype(rt8.sm, c(1L, 5L, 4L)))

identical(evalGenotype(rt8, c(2L, 5L, 4L)),
          evalGenotype(rt8.sm, c(2L, 5L, 4L)))

identical(evalGenotype(rt8, c(1L, 2L, 5L, 4L)),
          evalGenotype(rt8.sm, c(1L, 2L, 5L, 4L)))

         
## guardar todos los tests en un RData para future testing



## FIXME: store output of each wrap, and use as test cases later.





library(rbenchmark)

benchmark(to.long.rt(rt2), replications = 1000) ## 3.75
benchmark(to.long.rt(rt3), replications = 1000) ## 7.48

benchmark(rt.to.cpp(rt2), replications = 1000) ## 3.8
benchmark(rt.to.cpp(rt3), replications = 1000) ## 7.5
## so the R to CPP, in C++ part, is less than a milisecond
## and the R code is about 3 to 7 miliseconds. And this is
## only done once.






## For tests: this should fail it
rt1 <- data.frame(parent = c(
                      0, 0, 0,
                      1,
                      "2, 3",
                      4,
                      "5, 6, 7",
                      "5, 6, 7",
                      4
                      ),
           child = c(
               1,
               6,
               7,
               2,
               4,
               5,
               8,
               9,
               9),
                  s = 0.1,
                  sh = 0.05,
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

to.long.rt(rt1)





## FIXME  how to specify the genotype
## right now, it is doing it wrong! it is mapping letters to 1, 2, incorrectly
rt7.l <- data.frame(parent = c(
                      0, 0, "a", "b"
                      ),
           child = c(
               "a",
               "b",
               "c",
               "d"),
                  s = c(0.1, 0.2, 0.3, 0.4),
                  sh = c(99, 99, -0.05, -Inf),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(rt7.l, c(1L, 2L))
## FIXME  the next should work!!!
evalGenotype(rt7.l, c("a", "b"))

## synthetic viability
sv1 <- data.frame(parent = c(
                      0, "a"
                      ),
           child = c(
               "a",
               "b"),
                  s = c(-0.1, 0.2),
                  sh = c(99, -0.3),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(sv1, c(1L))
evalGenotype(sv1, c(2L))
evalGenotype(sv1, c(1L, 2L))



##

## A   B   fitness
## wt  wt  s
## wt  M   s1
## M   wt  s2
## M   M   s3

s1 <- -0.3
s2 <- -0.1
s3 <- 0.4

sv2a <- data.frame(parent = c(
                      0, "a"
                      ),
           child = c(
               "a",
               "b"),
                  s = c(s1, s3),
                  sh = c(99, s2),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(sv2a, c(1L))
evalGenotype(sv2a, c(2L))
evalGenotype(sv2a, c(1L, 2L))


sv2b <- data.frame(parent = c(
                      0, "b"
                      ),
           child = c(
               "b",
               "a"),
                  s = c(s2, s3),
                  sh = c(99, s1),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(sv2b, c(1L))
evalGenotype(sv2b, c(2L))
evalGenotype(sv2b, c(1L, 2L))


identical( evalGenotype(sv2a, c(2L)), evalGenotype(sv2b, c(1L)) )

evalGenotype(sv2a, c(1L))
evalGenotype(sv2b, c(2L)) ## FIXME this is wrong now. because of the letters?






## numbers not letters
sv2a <- data.frame(parent = c(
                      0, 1
                      ),
           child = c(
               1,
               2),
                  s = c(s1, s3),
                  sh = c(99, s2),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(sv2a, c(1L))
evalGenotype(sv2a, c(2L))
evalGenotype(sv2a, c(1L, 2L))


sv2b <- data.frame(parent = c(
                      0, 2
                      ),
           child = c(
               2,
               1),
                  s = c(s2, s3),
                  sh = c(99, s1),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)

evalGenotype(sv2b, c(1L))
evalGenotype(sv2b, c(2L))
evalGenotype(sv2b, c(1L, 2L))


identical( evalGenotype(sv2a, c(2L)), evalGenotype(sv2b, c(1L)) )

## these would be identical if we added, or similar, s and sh
evalGenotype(sv2a, c(1L))
evalGenotype(sv2b, c(1L)) 

## ditto
evalGenotype(sv2a, c(2L))
evalGenotype(sv2b, c(2L)) 


## works if we add the s.
evalGenotype(sv2a, c(1L, 2L))
evalGenotype(sv2b, c(1L, 2L)) 

## FIXME: evalGenotype should return the final fitness, final birth, and
## final death, if we would also specify the model


## Should we just add (or whatever) all the s and sh?
## Yes. Well, all the birth and death rate expressions are of the form
##  (1 +/- s)^number.of.events
## so just do \prod_{i} (1 + s_i), where s keeps the sign.



## FIXME allow s to be a symbol? Not, just define externally



## The general object has:
## - the restriction table
## - the epistasis
## - the order

## At least one of the above.



rt2 <- data.frame(parent = c(
                      0, 0, 0,
                      1,
                      "2, 3",
                      4,
                      "5, 6, 7",
                      "5, 6, 7",
                      4,
                      4
                      ),
           child = c(
               1,
               7,
               3,
               2,
               4,
               5,
               8,
               9,
               9,
               6),
                  s = 0.1,
                  sh = 0.05,
                  typeDep = "MN",
                  stringsAsFactors = FALSE)











## list.of.deps1 <- function(x) {
##     ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
##     ##                 "SM" = 2, "semimonotone" = 2)
##     lookupTypeDep <- c("MN" = "MN", "monotone" = "MN",
##                        "SM" = "SM", "semimonotone" = "SM",
##                        "XM" = "XM", "xmpn" = "XM")
##     ## FIXME: check values of typeDep

##     if(length(x) == 1)
##         return(list(
##             child = unique(as.integer(unlist(strsplit(x$child, ",")))),
##             s = x$s,
##             sh = x$sh,
##             typeDep = lookupTypeDep[x$typeDep],
##             parent = list(
##                 as.integer(unlist(strsplit(x$parent, ","))))))
##     else {
##         if(length(unique(x$s))!= 1)
##             stop("Not all s identical within a child")
##         if(length(unique(x$sh))!= 1)
##             stop("Not all sh identical within a child")
##         if(length(unique(x$typeDep))!= 1)
##             stop("Not all typeDep identical within a child")
##         return(list(
##             child = unique(as.integer(unlist(strsplit(x$child, ",")))),
##             s = x$s[1],
##             sh = x$sh[1],
##             typeDep = lookupTypeDep[x$typeDep[1]],
##             parent = lapply(strsplit(x$parent, ","), as.integer)))
##     }
## }









## list.of.deps0 <- function(x) {
##     ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
##     ##                 "SM" = 2, "semimonotone" = 2)
##     lookupTypeDep <- c("MN" = "MN", "monotone" = "MN",
##                        "SM" = "SM", "semimonotone" = "SM")
##     ## FIXME: check values of typeDep
##     if(length(x) == 1)
##         return(list(
##             child = x$child,
##             s = x$s,
##             sh = x$sh,
##             typeDep = lookupTypeDep[x$typeDep],
##             parent = list(
##                 as.integer(unlist(strsplit(x$parent, ","))))))
##     else {
##         if(length(unique(x$s))!= 1)
##             stop("Not all s identical within a child")
##         if(length(unique(x$sh))!= 1)
##             stop("Not all sh identical within a child")
##         if(length(unique(x$typeDep))!= 1)
##             stop("Not all typeDep identical within a child")
##         return(list(
##             child = x$child[1],
##             s = x$s[1],
##             sh = x$sh[1],
##             typeDep = lookupTypeDep[x$typeDep[1]],
##             parent = lapply(strsplit(x$parent, ","), as.integer)))
##     }
## }


## to.long.rt0 <- function(rt, verbosity = 0) {
##     if(is.numeric(rt$parent))
##         rt$parent <- as.character(rt$parent)
##     srt <- rt[order(rt$child), ]
##     ## check all childs
##     if(!identical(as.integer(sort(unique(rt$child))),
##                   seq.int(max(rt$child))))
##         stop("Not all children present")
##     if(verbosity >= 4)
##         message("Setting number of drivers to ",
##                 max(rt$child))
##     ## splitted <- split(srt, srt$child)
##     return(lapply(split(srt, srt$child), list.of.deps0))
## }





## list.of.deps.00 <- function(x) {
##     ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
##     ##                 "SM" = 2, "semimonotone" = 2)
##     lookupTypeDep <- c("MN" = "monotone",
##                        "monotone" = "monotone",
##                        "SM" = "semimonotone",
##                        "semimonotone" = "semimonotone")
##     ## FIXME: check values of typeDep
   
##     if(length(x) == 1)
##         return(list(
##             child = nice.string(x$child),
##             s = x$s,
##             sh = x$sh,
##             typeDep = lookupTypeDep[x$typeDep],
##             parent = nice.string(x$parent)))
##     else {
##         if(length(unique(x$s))!= 1)
##             stop("Not all s identical within a child")
##         if(length(unique(x$sh))!= 1)
##             stop("Not all sh identical within a child")
##         if(length(unique(x$typeDep))!= 1)
##             stop("Not all typeDep identical within a child")
##         return(list(
##             child = nice.string(x$child),
##             s = x$s[1],
##             sh = x$sh[1],
##             typeDep = lookupTypeDep[x$typeDep[1]],
##             parent = lapply(x$parent, nice.string)))
##     }
## }




## nice.vector.epist <- function(z) {
##     setdiff(sort(unique(unlist(lapply(strsplit(z, " "),
##                                     function(u) strsplit(u, ":"))))), "")
## }


## epist.element <- function(x, y) {
##     list(ids = nice.vector(x, sep = ":"), s = y)
## }
