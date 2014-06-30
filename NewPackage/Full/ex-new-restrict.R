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
                  type = "MN",
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
                  type = "MN",
                  stringsAsFactors = FALSE)


rt4 <- data.frame(parent = c(
                      0, 0
                      ),
           child = c(
               1,
               2),
                  s = 0.1,
                  sh = c(0.05, -Inf),
                  type = "MN",
                  stringsAsFactors = FALSE)

rt5 <- data.frame(parent = c(
                      0, 0
                      ),
           child = c(
               1,
               2),
                  s = -Inf,
                  sh = Inf,
                  type = "MN",
                  stringsAsFactors = FALSE)


list.of.deps <- function(x) {
    lookupType <- c("MN" = 1, "monotone" = 1,
                    "SM" = 2, "semimonotone" = 2)
    if(length(x) == 1)
        return(list(
            child = x$child,
            s = x$s,
            sh = x$sh,
            type = lookupType[x$type],
            parent = list(
                as.integer(unlist(strsplit(x$parent, ","))))))
    else {
        if(length(unique(x$s))!= 1)
            stop("Not all s identical within a child")
        if(length(unique(x$sh))!= 1)
            stop("Not all sh identical within a child")
        if(length(unique(x$type))!= 1)
            stop("Not all type identical within a child")
        return(list(
            child = x$child[1],
            s = x$s[1],
            sh = x$sh[1],
            type = lookupType[x$type[1]],
            parent = lapply(strsplit(x$parent, ","), as.integer)))
    }
}

to.long.rt <- function(rt) {
    if(is.numeric(rt$parent))
        rt$parent <- as.character(rt$parent)
    srt <- rt[order(rt$child), ]
    ## check all childs
    if(!identical(as.integer(sort(unique(rt$child))),
                  seq.int(max(rt$child))))
        stop("Not all children present")
    ## splitted <- split(srt, srt$child)
    return(lapply(split(srt, srt$child), list.of.deps))
}


rt.to.cpp <- function(rt) {
    lrt <- to.long.rt(rt)
    restrictTable_to_cpp0(lrt)
}

wrap.test.rt <- function(rt) {
    lrt <- to.long.rt(rt)
    wrap_test_rt(lrt)
}


library(Rcpp)
## setwd("../../")

sourceCpp("ex-multimap-for-restric.cpp",
          verbose = TRUE)

wrap.test.rt(rt3)
wrap.test.rt(rt2)

## test the Inf
wrap.test.rt(rt4)
wrap.test.rt(rt5)


rt.to.cpp(rt2)
rt.to.cpp(rt3)

wrap.test.rt(rt3)



## f4()
rt.to.cpp(rt2) ## will not work now, as the code for creating the object
               ## is outside restrictTable_to_cpp.

## turn into a list, each element of the list is the poset.
## similar to the C++ structure I'll have.

## How do I turn that into my structure inside C++?



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
                  type = "MN",
                  stringsAsFactors = FALSE)

to.long.rt(rt1)


