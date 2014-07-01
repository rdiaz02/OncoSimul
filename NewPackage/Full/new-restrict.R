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
    typeDep = "MN",
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
    paste(sort(unique(unlist(lapply(strsplit(z, " "),
                                    function(u) strsplit(u, ","))))),
          collapse = ", ")
}

list.of.deps <- function(x) {
    ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
    ##                 "SM" = 2, "semimonotone" = 2)
    lookupTypeDep <- c("MN" = "MN", "monotone" = "MN",
                       "SM" = "SM", "semimonotone" = "SM")
    ## FIXME: check values of typeDep
   
    if(length(x) == 1)
        return(list(
            child = nice.string(x$child),
            s = x$s,
            sh = x$sh,
            typeDep = lookupTypeDep[x$typeDep],
            parent = nice.string(x$parent)))
    else {
        if(length(unique(x$s))!= 1)
            stop("Not all s identical within a child")
        if(length(unique(x$sh))!= 1)
            stop("Not all sh identical within a child")
        if(length(unique(x$typeDep))!= 1)
            stop("Not all typeDep identical within a child")
        return(list(
            child = nice.string(x$child),
            s = x$s[1],
            sh = x$sh[1],
            typeDep = lookupTypeDep[x$typeDep[1]],
            parent = lapply(x$parent, nice.string)))
    }
}

gene.to.module <- function(rt) {
    gtm <- function(x) {
        data.frame(cbind(unlist(strsplit(x, ", ")), x))
    }
    all.modules <- unique(unlist(lapply(c(rt$parent, rt$child), nice.string)))
    geneMod <- as.data.frame(rbindlist(lapply(all.modules, gtm)))
    colnames(geneMod) <- c("Gene", "Module")
    geneMod$Gene <- as.character(geneMod$Gene)
    geneMod$Module <- as.character(geneMod$Module)
    geneMod <- geneMod[order(geneMod$Gene), ]
    geneMod$NumericID <- 0:(nrow(geneMod) - 1)
    geneMod
}



## FIXME: make sure mutations within modules are ordered!!
## This next add to R code.
## FIXME: remember to pass num drivers!!
to.long.rt <- function(rt, verbosity = 0) {
    if(is.numeric(rt$parent))
        rt$parent <- as.character(rt$parent)
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
    wrap_test_rt(lrt$long.rt, lrt$geneModule)
}

wrap.test.checkRestrictions <- function(rt, genotype) {
    lrt <- to.long.rt(rt)
    wrap_test_checkRestriction(lrt, genotype)
}



library(Rcpp)
## setwd("../../")

sourceCpp("new-restrict.cpp",
          verbose = TRUE)


wrap.test.rt(rt3)
wrap.test.rt(rt2)
wrap.test.rt(rt6)
wrap.test.rt(rt7)


## test the Inf
wrap.test.rt(rt4)
wrap.test.rt(rt5)


wrap.test.checkRestrictions(rt7, c(1L, 2L))
wrap.test.checkRestrictions(rt7, c(1L, 3L))
wrap.test.checkRestrictions(rt7, c(1L, 4L))
wrap.test.checkRestrictions(rt7, c(2L, 3L))
wrap.test.checkRestrictions(rt7, c(2L, 4L))

wrap.test.checkRestrictions(rt7, c(1L, 2L, 3L))
wrap.test.checkRestrictions(rt7, c(1L, 2L, 4L))
wrap.test.checkRestrictions(rt7, c(1L, 3L, 4L))

wrap.test.checkRestrictions(rt7, c(2L, 3L, 4L))
wrap.test.checkRestrictions(rt7, c(1L, 2L, 3L, 4L))

wrap.test.checkRestrictions(rt8, c(1L, 5L))
wrap.test.checkRestrictions(rt8, c(1L, 2L, 5L))

wrap.test.checkRestrictions(rt8, c(1L, 3L))
wrap.test.checkRestrictions(rt8, c(2L, 3L))
wrap.test.checkRestrictions(rt8, c(1L, 2L, 3L))
wrap.test.checkRestrictions(rt8, c(5L, 3L))
wrap.test.checkRestrictions(rt8, c(4L, 3L))

wrap.test.checkRestrictions(rt8, c(5L, 4L))
wrap.test.checkRestrictions(rt8, c(1L, 4L))
wrap.test.checkRestrictions(rt8, c(2L, 4L))
wrap.test.checkRestrictions(rt8, c(1L, 2L, 4L))






## do test with modules with multiple


## rt.to.cpp(rt2)
## rt.to.cpp(rt3)






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














list.of.deps1 <- function(x) {
    ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
    ##                 "SM" = 2, "semimonotone" = 2)
    lookupTypeDep <- c("MN" = "MN", "monotone" = "MN",
                       "SM" = "SM", "semimonotone" = "SM")
    ## FIXME: check values of typeDep

    if(length(x) == 1)
        return(list(
            child = unique(as.integer(unlist(strsplit(x$child, ",")))),
            s = x$s,
            sh = x$sh,
            typeDep = lookupTypeDep[x$typeDep],
            parent = list(
                as.integer(unlist(strsplit(x$parent, ","))))))
    else {
        if(length(unique(x$s))!= 1)
            stop("Not all s identical within a child")
        if(length(unique(x$sh))!= 1)
            stop("Not all sh identical within a child")
        if(length(unique(x$typeDep))!= 1)
            stop("Not all typeDep identical within a child")
        return(list(
            child = unique(as.integer(unlist(strsplit(x$child, ",")))),
            s = x$s[1],
            sh = x$sh[1],
            typeDep = lookupTypeDep[x$typeDep[1]],
            parent = lapply(strsplit(x$parent, ","), as.integer)))
    }
}









list.of.deps0 <- function(x) {
    ## lookupTypeDep <- c("MN" = 1, "monotone" = 1,
    ##                 "SM" = 2, "semimonotone" = 2)
    lookupTypeDep <- c("MN" = "MN", "monotone" = "MN",
                       "SM" = "SM", "semimonotone" = "SM")
    ## FIXME: check values of typeDep
    if(length(x) == 1)
        return(list(
            child = x$child,
            s = x$s,
            sh = x$sh,
            typeDep = lookupTypeDep[x$typeDep],
            parent = list(
                as.integer(unlist(strsplit(x$parent, ","))))))
    else {
        if(length(unique(x$s))!= 1)
            stop("Not all s identical within a child")
        if(length(unique(x$sh))!= 1)
            stop("Not all sh identical within a child")
        if(length(unique(x$typeDep))!= 1)
            stop("Not all typeDep identical within a child")
        return(list(
            child = x$child[1],
            s = x$s[1],
            sh = x$sh[1],
            typeDep = lookupTypeDep[x$typeDep[1]],
            parent = lapply(strsplit(x$parent, ","), as.integer)))
    }
}


to.long.rt0 <- function(rt, verbosity = 0) {
    if(is.numeric(rt$parent))
        rt$parent <- as.character(rt$parent)
    srt <- rt[order(rt$child), ]
    ## check all childs
    if(!identical(as.integer(sort(unique(rt$child))),
                  seq.int(max(rt$child))))
        stop("Not all children present")
    if(verbosity >= 4)
        message("Setting number of drivers to ",
                max(rt$child))
    ## splitted <- split(srt, srt$child)
    return(lapply(split(srt, srt$child), list.of.deps0))
}
