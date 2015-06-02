## data with genes in columns, subjects in rows
x1 <- replicate(100, sample(c(0, 1), 10, p = c(0.8, 0.2), replace = TRUE)); dim(x1)
x1 <- t(x1)
## play with column names
nn <- sample(1:100, ncol(x1))
colnames(x1) <- paste("g", nn, sep = "_")



##### Oncogenetic trees, Szabo & Boucher
library(Oncotree)

onco.fit <- oncotree.fit(x1)
## root gets number 1
onco.fit$parent$parent.num
## convert to adjacency matrix
## but remember 1 is actually the root, so we add a 1, and later rename
p.to.child <- cbind( onco.fit$parent$parent.num[-1], 2:(ncol(x1) + 1))
m1 <- matrix(0, nrow = ncol(x1)+1, ncol = ncol(x1)+1)
m1[p.to.child] <- 1
rownames(m1) <- colnames(m1) <- 0:ncol(x1) ## the adjacency matrix.

## check figure
par(mfrow = c(2, 2))
plot(onco.fit)
library(igraph)
gg1 <- graph.adjacency(m1, mode = "directed")
plot(gg1, vertex.label = colnames(m1))
library(gRbase)
g1 <- as(m1, "graphNEL")
plot(g1)


run.oncotree(x1)

######################### Rtreemix
library(Rtreemix)
rtm <- new("RtreemixData", Sample = x1, Events = as.character(0:ncol(x1)))
o1 <- fit(data = rtm, K = 1, noise = FALSE) ## same as onco.fit
o2 <- fit(data = rtm, K = 3, noise = TRUE, equal.edgeweights = TRUE)

## these already contain graphnel objects.

run.rtreemix(x1, 3, noise = TRUE)


####### CBN

run.cbn(x1, file = "testcbn")





## selection
x1 <- replicate(100, sample(c(0, 1), 10, p = c(0.8, 0.2), replace = TRUE))
selectSB(x1, threshold = 5)



##### Clonal ordering

run.clonalordering <- function(x, p.v.thresh = 0.05,
                               p.adjust.method = "none",
                               type.out = "graphNEL") {
  ## This just establishes an order, since higher freq. ought to happen
  ## earlier, but the actual counts alllow for reversion
  omuts <- sort(colSums(x), decreasing = TRUE)
  ngenes <- length(omuts)
  namesg <- names(omuts)
  lps <- matrix(0, nrow = choose(ngenes, 2), ncol = 3)
  lpn <- matrix("no_name", nrow = choose(ngenes, 2), ncol = 2)
  k <- 1
 
  for(i in 1:(ngenes - 1)) {
    for(j in (i + 1):ngenes) {
      tmpx <- table(x[, i], x[, j])[c(2, 3)]
      tmp.p <- binom.test(tmpx)$p.value
      if(tmpx[1] > tmpx[2]) {
        lpn[k, ] <- c(namesg[i], namesg[j])
        lps[k, ] <- c(i, j, tmp.p)
      } else {
        cat("\n NOTE: reversion w.r.t. to total counts")
        lpn[k, ] <- c(namesg[j], namesg[i])
        lps[k, ] <- c(j, i, tmp.p)
      } 
      k <- k + 1
    }
  } 
  if(p.adjust.method != "none")
    lps[, 3] <- p.adjust(lps[, 3], p.adjust = p.adjust.method)
  i.select <- which(lps[, 3] < p.v.thresh)
  lps <- lps[i.select, ]
  lpn <- lpn[i.select, ]
  browser()
  ## lps contains the posets
  poset.to.graph(lps[, c(1, 2)], names = c("Root", namesg),
                 type = type.out)
}

## Trying to implement logic in Barrett et al., 1999.
##  - for each subject, find order of mutations (mut.order)
##  - with data for all subjects, do a test.
##  - these data are not timed. I.e., we do not take some samples
##    before other samples.

## But we could do that: have mut.order.timed and see which comes first.


mut.order <- function(x, y) {
  ## return codes
  ## 4: conflicting evidence or no evidence. Like an NA
  ## 3: same time
  ## 1: x before y
  ## 2: y before x
  t1 <- table(x, y)
  if( (t1[2] > 0) && (t1[3] > 0) ) {
    return(4)
  } else if( (t1[2] > 0) && (t1[4] > 0) ) {
    return(1)
  } else if( (t1[3] > 0) && (t1[4] > 0) ) {
    return(2)
  } else if( (t1[4] > 0)) { ## and t1[3] == 0 and t1[2] == 0
    return(3)
  } else return(4)
}

indiv.clonalordering.notime <- function(x) {
  ## combs <- t(combn(ncol(x), 2))
  combs <- t(combn(sort(colnames(x)), 2))
  oo <- apply(combs, 1, function(z) mut.order(x[, z[1]], x[ , z[2]]))
  return(data.frame(combs, oo))
}


clonal.ordering <- function(x) {
  ## x is a list of observations (rows) by genes (columns)
  ## Each list is for a different subject
  ## But each list should have the same genes.

}

## Simulation return: observations by genes.
## Extra columns:
##    - subject
##    - time
##    - model code as string

## An extra column: subject (which might be )







## run.cbn <- function(data, file = "testcbn", eparam = 0.05,
##                     temp = 10, steps = 200) {
## ## I assume h-cbn and ct-cbn are available in local dir
##   write.linear.poset(data, file)


##   ## writeLines(as.character(c(ncol(data), 0)),
##   ##            con = paste(file, ".poset", sep = ""))


  
##   data2 <- cbind(1, data)
##   write(c(nrow(data2), ncol(data2)),
##         file = paste(file, ".pat", sep = ""),
##         sep = " ")
##   write(t(data2), file = paste(file, ".pat", sep = ""),
##         ncolumns = ncol(data2),
##         append = TRUE, sep = " ")
##   write(c("dummy", colnames(data)),
##         file = paste(file, ".prf", sep = ""),
##         sep = " ")
  
##   system(paste("export OMP_NUM_THREADS=", detectCores(), sep = ""))
##   dir.create(file)
##   ## first call is to create the lambda file
##   ## I think I did this two first calls to get an initial poset
##   ## from ct-cbn
##   system(paste("./h-cbn -f",  file, "-w"))
##   cat("\n\n")
##   system(paste("./h-cbn -f",  file, "-e", eparam,
##                "-w -m"))
##   cat("\n\n")
##   system(paste("./h-cbn -f",  file, "-s", 
##                "-T", temp,  "-N", steps,
##                "-m -w"))

##   ## the final poset in file/00000.poset
  
## }




## a direct translation of linear_poset and write_poset
linear.poset.RDU <- function(x) {
  nr <- nrow(x)
  nc <- ncol(x)
  sorted <- order(colMeans(x), decreasing = TRUE)
  poset <- matrix(0, ncol = nc, nrow = nc)
  s <- sorted[1]
  for (t in sorted[2:nc]) {
    poset[s, t] <- 1
    s <- t
  }
  
  ## now, translate write_poset.
  ## posetw <- matrix(0, ncol = 2, nrow = nr)
  ## for (i in 1:nc) {
  ##   for(j in 1:nc) {
  ##     if(poset[i, j])
  ##       posetw[i, ] <- c(i, j)
  ##   }
  ## }
  
  ## do the R way
  posetw <- which(poset == 1, arr.ind = TRUE)
  posetw <- posetw[order(posetw[, 1]), ]
  return(posetw)
}



linear.poset.Ricardo <- function(x) {
  cols <- ncol(x) + 1 ## unless we have added the 1
  m1 <- matrix(0, ncol = cols -1 , nrow = cols -1)
  o1 <- order(colMeans(x), decreasing = TRUE)
  s <- o1[1]
  v1 <- vector()
  v2 <- vector()
  ## change loop and "c" for something better
  for(t in o1[2:(cols - 1)]) {
    v1 <- c(v1, s)
    v2 <- c(v2, t)
    m1[s, t] <- 1
    s <- t
  }
  vsorted <- order(v1, decreasing = FALSE)
  m2 <- matrix(0, ncol = 2, nrow = cols)
  v1 <- v1[vsorted]
  v2 <- v2[vsorted]
  m2[1, 1] <- cols - 1

  for(x in 1:(cols - 2)) {
    m2[ x + 1, 1] <- v1[x]
    m2[ x + 1, 2] <- v2[x]
  }
  return(m2)
}







### in the Pathways example dir, the following might help understand
### if we do not use the linear posets created, it takes forever.

### (You can create the two-liner by hand; or you can run the command h-cbn -f jones_core -w
### but this will also be veeery slow, and I think gives you nothing extra)


### If a posets file with entries existes, then the command
### h-cbn -f jones-core is fast. If it is just the simple two-row file, it takes forever.

### Now, use two different poset files, and compare what happens
### posets.linear comes from the python run. posets.nada is the two-liner


### cp jones_core.poset.linear jones_core.poset; ../../Oncogenetic-trees/h-cbn -f jones_core -s -N 3 -m -v

### cp jones_core.poset.nada jones_core.poset; ../../Oncogenetic-trees/h-cbn -f jones_core -s -N 3 -m -v


### the "-v" makes it verbose, and you see in what it has to spend a lot of time



### I probably want to understand what linear_poset (in cbn.py) is doing, and probably
### call it directly from R using and appropriate pat file.








##################################################

### How many cases shared if sampling from a pool of simulations


p.e.shared <- function(N, n, m, e) {
  ## N: size of pool
  ## n: size of sample 1
  ## m: size of sample 2
  ## e: number shared
  cc <- choose(m, e)

  num <- prod(seq(from = n, to = n - e + 1)) *
    prod( seq(from = N-n, to = N - n - (m - e) + 1) )

  den <- prod(seq(from = N, to = N - m +1 ))

  return(cc * num/den)
}
## NaNs (from Infs) for large values


i.s0 <- function(N, n, m) {
  s1 <- sample(1:N, n, replace = FALSE)
  s2 <- sample(1:N, m, replace = FALSE)
  return(length(intersect(s1, s2)))
}



p.s.sim <- function(N, n, m, reps = 1e4) {
  ints <- replicate(reps, i.s0(N, n, m))
  tt <- table(ints)/reps
  ev <- mean(ints)
  return(list(tt, ev))
}


p.s.sim2 <- function(N, n, m, reps = 1e4) {
  s1 <- sample(1:N, n, replace = FALSE)
  
  i.s0b <- function(N, m) {
    s2 <- sample(1:N, m, replace = FALSE)
    return(length(intersect(s1, s2)))
  }

  ints <- replicate(reps, i.s0b(N, m))
  tt <- table(ints)/reps
  ev <- mean(ints)
  if(is.na(tt["0"])) {
     p.a.o <- sum(tt)
   } else {
     p.a.o <- 1 - tt["0"]
   }
  return(list(table.shared = tt,
              prob.at.least.one.shared = p.a.o,
              expected.number.shared = ev))
}


p.s.sim3 <- function(N, n, m, reps = 1e4) {
  s1 <- sample(1:N, n, replace = FALSE)
  
  i.s0b <- function(N, m) {
    s2 <- sample(1:N, m, replace = FALSE)
    return(length(intersect(s1, s2)))
  }

  ints <- unlist(mclapply(1:reps,
                          function(x) i.s0b(N, m),
                          mc.cores = detectCores()))

  tt <- table(ints)/reps
  ev <- mean(ints)
  if(is.na(tt["0"])) {
     p.a.o <- sum(tt)
   } else {
     p.a.o <- 1 - tt["0"]
   }
  return(list(table.shared = tt,
              prob.at.least.one.shared = p.a.o,
              expected.number.shared = ev))
}


## compare, e.g.,
p.e.shared(40, 7, 5, 2)
# with third entry in table
p.s.sim2(40, 7, 5, 1e5)

## now, examples
p.s.sim2(1000, 10, 10, 1e5)


p.s.sim2(500000, 1000, 1000, 1e4)
## > $table.shared
## ints
##      0      1      2      3      4      5      6      7      8      9 
## 0.1335 0.2658 0.2766 0.1827 0.0877 0.0378 0.0115 0.0039 0.0003 0.0002 

## $prob.at.least.one.shared
##      0 
## 0.8665 

## $expected.number.shared
## [1] 2.007

p.s.sim2(1000000, 1000, 1000, 1e4)
## > $table.shared
## ints
##      0      1      2      3      4      5      6      7 
## 0.3756 0.3619 0.1835 0.0600 0.0144 0.0039 0.0006 0.0001 

## $prob.at.least.one.shared
##      0 
## 0.6244 
library(parallel)
## $expected.number.shared
## [1] 0.9903
p.s.sim3(1000000, 1000, 1000, 1e4)
## ints
##      0      1      2      3      4      5      6      8 
## 0.3671 0.3698 0.1806 0.0647 0.0147 0.0022 0.0008 0.0001 

## $prob.at.least.one.shared
##      0 
## 0.6329 

## $expected.number.shared
## [1] 1
p.s.sim3(1000000, 1000, 1000, 1e5)
## $table.shared
## ints
##       0       1       2       3       4       5       6       7       8 
## 0.36742 0.36846 0.18307 0.06176 0.01557 0.00321 0.00042 0.00006 0.00003 

## $prob.at.least.one.shared
##       0 
## 0.63258 

## $expected.number.shared
## [1] 1.00139
p.s.sim3(1000000, 1000, 1000, 5e5)
## $table.shared
## ints
##        0        1        2        3        4        5        6        7 
## 0.367196 0.369252 0.184014 0.060470 0.015390 0.003084 0.000502 0.000082 
##        8        9 
## 0.000008 0.000002 

## $prob.at.least.one.shared
##        0 
## 0.632804 

## $expected.number.shared
## [1] 0.999338


p.s.sim3(500000, 1000, 1000, 5e5)
## $table.shared
## ints
##        0        1        2        3        4        5        6        7 
## 0.135080 0.270344 0.271254 0.180820 0.090178 0.036272 0.011574 0.003388 
##        8        9       10       11       12       13 
## 0.000858 0.000194 0.000032 0.000002 0.000002 0.000002 

## $prob.at.least.one.shared
##       0 
## 0.86492 

## $expected.number.shared
## [1] 1.999546


p.s.sim3(1000000, 1000, 200, 5e5)
## $table.shared
## ints
##        0        1        2        3        4        5 
## 0.818432 0.164238 0.016114 0.001168 0.000046 0.000002 

## $prob.at.least.one.shared
##        0 
## 0.181568 

## $expected.number.shared
## [1] 0.200164


p.s.sim3(1000000, 1000, 50, 5e5)

## $table.shared
## ints
##        0        1        2        3        4 
## 0.951400 0.047436 0.001144 0.000018 0.000002 

## $prob.at.least.one.shared
##      0 
## 0.0486 

## $expected.number.shared
## [1] 0.049786
