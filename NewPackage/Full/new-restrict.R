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
                  s = c(0.1, 0.3, 0.2, 0.25),
                  sh = c(99, 99, -0.05, -Inf),
                  typeDep = "MN",
                  stringsAsFactors = FALSE)






list.of.deps <- function(x) {
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

## FIXME: make sure mutations within modules are ordered!!
## This next add to R code.
## FIXME: remember to pass num drivers!!
to.long.rt <- function(rt, verbosity = 0) {
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
    return(lapply(split(srt, srt$child), list.of.deps))
}


## rt.to.cpp <- function(rt) {
##     lrt <- to.long.rt(rt)
##     rTable_to_Poset0(lrt)
## }

wrap.test.rt <- function(rt) {
    lrt <- to.long.rt(rt)
    wrap_test_rt(lrt)
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


