## See if I get right using a Moran process with Mather et al., 2012, approach.

## I try to differentiate the functino to compute from the computed value.
## I add a ".f"


library(rbenchmark)





W.f <- function(death, growth, mu) return(death + growth + mu)



R.f <- function(death, growth, mu)
  return(sqrt( (growth - death)^2 + (2*growth + 2*death + mu)*mu ) )



pM.f <- function(t, R, W, death, growth) {

## FIXME: large R * t will lead to overflows here.  
  Ct <- cosh(R * t /2)
  St <- sinh(R * t /2)

  if(is.infinite(Ct) || is.infinite(St)) {
    cat("\n\n Ct or St infinite: Ct = ", Ct, "  St", St, "\n")
    browser()
  }
##  return( (R*Ct + 2 * death * St - W * St)/(R*Ct - 2 * growth * St + W * St)     )
## factor out a few things

  ## I bet value.o is much more stable than value.
##  value.o <- (R*Ct + 2 * death * St - W * St)/(R*Ct - 2 * growth * St + W * St)     
  value <- (R * Ct + St * (2 * death - W ))/(R * Ct + St * (W - 2 * growth))
  if(is.nan(value)) {
    cat("\n\n pM is NAN \n")
    browser()
  }
  return(value)  
  
  ## return( (R * Ct + St * (2 * death - W ))/(R * Ct + St * (W - 2 * growth))     )
  
}




pE.f <- function(pM, W, death, growth) {
  return( (death * (1 - pM ) )/(W - death - growth * pM ) )
  
}



pB.f <- function(pE, death, growth) return((growth * pE)/death)



basic.values <- function(t, growth, death, mu) {
    W <- W.f(death, growth, mu)
    R <- R.f(death, growth, mu)
    pM <- pM.f(t, R, W, death, growth)
    pe <- pE.f(pM, W, death, growth)
    return(c(pe = pe, 
             pm = pM,
             pb = pB.f(pe, death, growth)))


}

algo2.ex <- function(t, num, growth, death, mutation) {
    W <- W.f(death, growth, mutation)
    R <- R.f(death, growth, mutation)
    pm <- pM.f(t, R, W, death, growth)
    pe <- pE.f(pm, W, death, growth)
    pb <- pB.f(pe, death, growth)

  if( ((1 - pe/pm) > 1) || ((1 - pe/pm) < 0) || (pb > 1) || (pb <= 0) ) {
    cat("\n Algo 2, pe, pm, pb in the limit\n")
    print(c(pe, pm, pb))
  }

    
    m <- rbinom(n = 1, size = num , prob = (1 - (pe/pm) ))
    if(m == 0) {
        n <- 0
    } else {
        n <- m + rnbinom(n = 1, size = m, prob = 1 - pb)
    }
    return(c(W = W, R = R, pm = pm, pe = pe , pb = pb,
             m  = m, n = n))
}


ti.ex <- function(num, growth, death, mutation) {
    W <- W.f(death, growth, mutation)
    R <- R.f(death, growth, mutation)
    return(ti.original(R, W, death, growth, num))
}

ti.log <- function(R, W, death, growth, n) {
  ## new version, using logs
  
  ## Can be made faster by generating a bunch of runifs in a call?
  ## Yes, FIXME. Use that later in C or even in R
  ## Not really: only when we initialize or sample. Otherwise,
  ## we draw one random number only
  r <- runif(1)

  ## rr is r ^ (1/n), computed in a more stable (?) way
  rr <- exp((1/n) * log(r))


  ## I stop here! Not clear this improves anything
  if( (((R - W + 2 * death)/(R + W - 2 * growth))^n)   < r ) {
    ## I add the max, to prevent negative numbers
    eq.11 <- (1/R) * log(
                       ((r^(1/n)) * (R - W + 2 * growth) - W - R + 2 * death) /
                       ((r^(1/n)) * (-R - W + 2 * growth) - W + R + 2 * death))
    ## FIXME: what do they do in their code???!!!  They do nothing
    if (eq.11 <= 0)
      browser()
    return(max(0, eq.11))
    
  } else {
    return(Inf)
  } 
}



ti.original <- function(R, W, death, growth, n) {
  ## a hack, to signal mutation is cero, and thus ti is Inf
  if( (R == -99) && (W == -99)) {
    ## the "cannot mutation" condition
    return(Inf)
  }


## C++
  ## double a = std::numeric_limits<double>::infinity();
  
## we could use simply W, which is necessarily possitive
  ## or just R, which is also necessarily positive
  
## ti is begin called in an mapply, and that signal is being passed in Algo5.
  
  
  ## original version, without using logs.
  
  ## Can be made faster by generating a bunch of runifs in a call?
  ## Yes, FIXME. Use that later in C or even in R
  ## Not really: only when we initialize or sample. Otherwise,
  ## we draw one random number only
  r <- runif(1)
  if( (((R - W + 2 * death)/(R + W - 2 * growth))^n)   < r ) {
    ## I add the max, to prevent negative numbers
    eq.11 <- (1/R) * log(
                       ((r^(1/n)) * (R - W + 2 * growth) - W - R + 2 * death) /
                       ((r^(1/n)) * (-R - W + 2 * growth) - W + R + 2 * death))
    ## FIXME: what do they do in their code???!!!  They do nothing
    if (eq.11 <= 0) {
      cat("\n\n eq.11 <= 0. eq.11 = ", eq.11, "\n")
      browser()
    }
    return(max(0, eq.11))
    
  } else {
    return(Inf)
  } 
}

## for now, leave as it was
ti <- ti.original



ti.divide <- function(R, W, death, growth, n) {
  ## a hack, to signal mutation is cero, and thus ti is Inf
  if( (R == -99) && (W == -99)) {
    ## the "cannot mutation" condition
    return(Inf)
  }


## C++
  ## double a = std::numeric_limits<double>::infinity();
  
## we could use simply W, which is necessarily possitive
  ## or just R, which is also necessarily positive
  
## ti is begin called in an mapply, and that signal is being passed in Algo5.
  
  
  ## original version, without using logs.
  
  ## Can be made faster by generating a bunch of runifs in a call?
  ## Yes, FIXME. Use that later in C or even in R
  ## Not really: only when we initialize or sample. Otherwise,
  ## we draw one random number only
  r <- runif(1)
  if( (((R - W + 2 * death)/(R + W - 2 * growth))^n)   < r ) {
    ## I add the max, to prevent negative numbers
    eq.11 <- (1/R) * (log(-1 * ((r^(1/n)) * (R - W + 2 * growth) - W - R + 2 * death )) -
                      log(-1 * ((r^(1/n)) * (-R - W + 2 * growth) - W + R + 2 * death)))
    ## FIXME: what do they do in their code???!!!  They do nothing
    if (eq.11 <= 0) {
      cat("\n\n eq.11 <= 0. eq.11 = ", eq.11, "\n")
      browser()
    }
    return(max(0, eq.11))
    
  } else {
    return(Inf)
  } 
}





Algo2 <- function(num, t, R, W, death, growth) {

  ## cat("\n Arguments", num, t, R, gm, W, g, "\n")
  pm <- pM.f(t, R, W, death, growth)
  pe <- pE.f(pm, W, death, growth)
  pb <- pB.f(pe, death, growth)

  if(any(is.na(c(pm, pe, pb)))) browser()

##  browser()
  ## if( (1 - pe/pm) > 1) browser()
  ## if( (1 - pe/pm) < 0) browser()
  ## if(pb > 1) browser()
  ## if(pb <= 0) browser()


  ## FIXME: ask authors about underflow, here. They make no provission for it in their code. 

  if(isTRUE(all.equal(pe, pm))) {
    cat("\n WARNING Algo2: pe == pm\n")
    return(0) ## n cannot differ from 0 
  }

  
  ## if( ((1 - pe/pm) > 1) || ((1 - pe/pm) < 0) || (pb > 1) || (pb <= 0) ) {
  ##   cat("\n Algo 2, pe, pm, pb in the limit\n")
  ##   print(c(pe, pm, pb))
  ##   browser()
  ## }


  if( ((1 - pe/pm) > 1) || ((1 - pe/pm) < 0) || (pb > 1) || (pb <= 0) ) {
    cat("\n Algo 2, pe, pm, pb in the limit\n")
    print(c(pe, pm, pb))
    browser()
  }

  
## FIXME: we probably want to use all.equal in R and set to 0 if
  ## all.equal(pe/pm), etc.
  
  
  if(num < 0) browser()
  m <- rbinom(n = 1, size = num , prob = (1 - (pe/pm) ))

  if(m == 0) {
    n <- 0
  } else {
    n <- m + rnbinom(n = 1, size = m, prob = 1 - pb)
  }

##  cat("\n            Algo. 2.   m = ", m, ".  n = ", n)
  return(n)
}



Algo3 <- function(num, t, R, W, death, growth) {
  pm <- pM.f(t, R, W, death, growth)
  pe <- pE.f(pm, W, death, growth)
  pb <- pB.f(pe, death, growth)

  if(any(is.na(c(pm, pe, pb)))) browser()

  if(isTRUE(all.equal(pe, pm))) {
    cat("\n WARNING Algo3: pe == pm\n")
    return(0) ## n cannot differ from 0 
  }

##  browser()
  ## checks
  if( (1 - pe/pm) > 1) browser()
  if( (1 - pe/pm) < 0) browser()
  if(pb > 1) browser()
  if(pb <= 0) browser()
  
  if(num < 1) browser()
  ## cat(" Inside algo 3", "num ", num, "prob ", (1 - (pe/pm) ),
  ##     "pe ", pe, "pm", pm)
  ## pex <<- pe
  ## pmx <<- pm
  ## sizex <<- num - 1
  m <- rbinom(n = 1, size = num - 1, prob = (1 - (pe/pm) ))

  
  n <- m + 1 + rnbinom(n = 1, size = m + 2, prob = 1 - pb)
##  cat("\n          Algo. 3.   m = ", m, ".  n = ", n)
  return(n)
}



## mutations start with 0. So 0 is a mutation
## NOPE! Makes no sense


fitness.CBN2 <- function(mutatedPos, genotype, restrict.table, num.drivers,
                         birth.rate, s, fitness.parent, CBN.multiple = "CBN",
                         fitnessNo = 0, fitnessYes = fitness.linear) {

  if(mutatedPos > num.drivers) { # the new mutation is a passenger
    return(fitness.parent)
    ##    return(fitnessYes(genotype, birth.rate, s, num.drivers))
  } else {
    sum.present <- 0
    num.dependencies <- restrict.table[mutatedPos, 2]
    if(!num.dependencies) {
      out <- "Yes"
    } else {
      out <- "No"
      for(i in 3:(2 + num.dependencies)) { ## ojito indices in C
        sum.present <- sum.present + genotype[restrict.table[mutatedPos, i]]
      }
      if(CBN.multiple == "Multiple") {
        if(sum.present) out <- "Yes" #
      } else{ ## if(CBN.multiple = "CBN")
        if(sum.present == restrict.table[mutatedPos, 2])
          out <- "Yes" # return(fitnessYes(genotype, birth.rate, s, num.drivers))
      }
    }
    if(out == "Yes") return(fitnessYes(genotype, birth.rate, s, num.drivers))
    else return(fitnessNo)
  }
##  return(fitnessNo)
}



fitness.CBN <- function(mutatedPos, genotype, restrict.table, num.drivers,
                         birth.rate, s, CBN.multiple = "CBN",
                         fitnessNo = 0, fitnessYes = fitness.linear) {
## withouth fitness.parent
  if(mutatedPos > num.drivers) { # the new mutation is a passenger
##    return(fitness.parent)
    return(fitnessYes(genotype, birth.rate, s, num.drivers))
  } else {
    sum.present <- 0
    num.dependencies <- restrict.table[mutatedPos, 2]
    if(!num.dependencies) {
      out <- "Yes"
    } else {
      out <- "No"
      for(i in 3:(2 + num.dependencies)) { ## ojito indices in C
        sum.present <- sum.present + genotype[restrict.table[mutatedPos, i]]
      }
      if(CBN.multiple == "Multiple") {
        if(sum.present) out <- "Yes" #
      } else{ ## if(CBN.multiple = "CBN")
        if(sum.present == restrict.table[mutatedPos, 2])
          out <- "Yes" # return(fitnessYes(genotype, birth.rate, s, num.drivers))
      }
    }
    if(out == "Yes") return(fitnessYes(genotype, birth.rate, s, num.drivers))
    else return(fitnessNo)
  }
##  return(fitnessNo)
}


fitness.linear <- function(newGenotype, basal.birth.rate, s = 0.005, num.drivers = 20) {
  ## later, we will use num.drivers. For now not
  ## (1 + s)^(sum of mutations) as in beerenw, but need to have 1 be different.
  ## and this is exponential. Use a simple linear model
  return(basal.birth.rate + s * sum(newGenotype[1:num.drivers]))
}







#### Comparisons and benchmarks


## comparing the ti functions. 
ss <- 100; set.seed(ss); d <- 0.3; g <- 0.300005; mu <- 0.005; n <- 300; ti.c(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g, n, 99); set.seed(ss); ti.original(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g, n)


ss <- 1; set.seed(ss); d <- 0.3; g <- 0.300005; mu <- 0.005; n <- 300; exwrap(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g, n, 99); set.seed(ss); ti.original(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g, n)



d <- 0.3; g <- 0.300005; mu <- 0.005; n <- 300;
## R <- R.f.c(d, g, mu, 0);
## W <- W.f.c2(d, g, mu, 0)
R <- R.f(d, g, mu)
W <- W.f(d, g, mu)

set.seed(1);
ti.original(R, W, d, g, n)
set.seed(1);
wrapTi(R, W, d, g, n, -99)

benchmark(ti.original(R, W, d, g, n), 
          wrapTi(R, W, d, g, n, -99),
          columns=c("test", "replications", "elapsed", "relative", "user.self", "sys.self"),
          order="relative",
          replications=50000)

d <- 0.3; g <- 0.405; mu <- 0.005; n <- 30;
## R <- R.f.c(d, g, mu, 0);
## W <- W.f.c2(d, g, mu, 0)

R <- R.f(d, g, mu)
W <- W.f(d, g, mu)

t <- 20

set.seed(10)
Algo2(n, t, R, W, d, g)
set.seed(10)
wrapAlgo2(n, t, R, W, d, g, 1, -99)
set.seed(10)
wrapAlgo2(n, t, R, W, d, g, 0, -99)

set.seed(3)
Algo3(n, t, R, W, d, g)
set.seed(3)
wrapAlgo3(n, t, R, W, d, g, 1, -99)
set.seed(3)
wrapAlgo3(n, t, R, W, d, g, 0, -99)



debug <- 1

benchmark(Algo2(n, t, R, W, d, g), 
          wrapAlgo2(n, t, R, W, d, g, debug, -99),
          columns=c("test", "replications", "elapsed", "relative"),
          order="relative",
          replications=10000)

benchmark(Algo3(n, t, R, W, d, g), 
          wrapAlgo3(n, t, R, W, d, g, debug, -99),
          columns=c("test", "replications", "elapsed", "relative"),
          order="relative",
          replications=10000)




## if not "as.integer", no in-place modifications!
## example
##
m1 <- matrix(rep(1, 15), ncol = 3); class(m1[, 1])
## and pointer location changes
wrapFitnessLinearVerbose(m1, 2, 0.3, 0.005, 4, 99); m1
wrapFitnessLinearVerbose(m1, 2, 0.3, 0.005, 4, 99); m1
gc()
wrapFitnessLinearVerbose(m1, 2, 0.3, 0.005, 4, 99); m1

m2 <- matrix(rep(1L, 15), ncol = 3); class(m2[, 1])
## and pointer location changes
wrapFitnessLinearVerbose(m2, 2, 0.3, 0.005, 4, 99); m2
wrapFitnessLinearVerbose(m2, 2, 0.3, 0.005, 4, 99); m2
gc()
wrapFitnessLinearVerbose(m2, 2, 0.3, 0.005, 4, 99); m2



allGenot <- matrix(as.integer(sample(c(0, 1), 15, 0.5)), ncol = 3)
allGenot
wrapFitnessLinear(allGenot, 0, 0.3, 0.005, 4, 99)
wrapFitnessLinearVerbose(allGenot, 0, 0.3, 0.005, 4, 99)
allGenot ## this ain't modified! it was not returned
## R creates a tmp when an SEXP is created when calling the function?
## note what happens with the position with and without gc()

wrapFitnessLinear(allGenot, 1, 0.3, 0.005, 4, 99)
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 4, 99)

wrapFitnessLinear(allGenot, 0, 0.3, 0.005, 3,  99)
wrapFitnessLinear(allGenot, 1, 0.3, 0.005, 3, 99)
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 3, 99)

wrapFitnessLinear(allGenot, 0, 0.3, 0.005, 2, 99)
wrapFitnessLinear(allGenot, 1, 0.3, 0.005, 2, 99)
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 2, 99)

wrapFitnessLinear(allGenot, 0, 0.3, 0.005, 1, 99)
wrapFitnessLinear(allGenot, 1, 0.3, 0.005, 1, 99)
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 1, 99)

wrapFitnessLinear(allGenot, 0, 0.3, 0.005, 0, 99)
wrapFitnessLinear(allGenot, 1, 0.3, 0.005, 0, 99)
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 0, 99)



## bound check? Yes, on column
wrapFitnessLinear(allGenot, 3, 0.3, 0.005, 2, 99)
wrapFitnessLinear(allGenot, -9, 0.3, 0.005, 2, 99)

## not on accumulate
wrapFitnessLinear(allGenot, 2, 0.3, 0.005, 77, 99)





##
restrictTable <- rbind(
  c(1, 0, -9, -9),
  c(2, 0, -9, -9),
  c(3, 1, 1, -9),
  c(4, 1, 3, -9),
  c(5, 2, 4, 2)
  )


convertRestrictTable <- function(x) {
  ## to convert the table to the format for C
  ## as there the mutations are numbered from 0

  ## In R the format for a row is:
  ##  - the mutation,
  ##  - the number of mutations on which it depends
  ##  - the actual mutations on which it depends
  ##  - the rest are "-9"
  
  t.restrictTable <- matrix(as.integer(x),
                            ncol = nrow(x), byrow = TRUE)

  t.restrictTable[-2, ] <- t.restrictTable[-2, ] - 1
  return(t.restrictTable)
}

restrictTableC <- convertRestrictTable(restrictTable)

## note genotypes has subject as column
genotypes <- t(
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 1, 1, 0, 0),
    c(0, 1, 1, 1, 0, 0),
    c(0, 1, 0, 1, 1, 0),
    c(0, 0, 1, 1, 1, 0),
    c(0, 0, 1, 1, 1, 1),
    c(1, 0, 1, 1, 1, 1))
  )


## OJO: in the C code, positions start at 0.
wrapFitnessCBNRcpp(1, genotypes, 0, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNArma(1, genotypes, 0, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNstd(1, genotypes, 0, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)

## the mutated position, the genotype, .., num.of drivers
## for the wrap, genotypes, the actual genotype col
fitness.CBN2(2, c(0, 1, 0, 0, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")
fitness.CBN2(2, genotypes[, 1], restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")



wrapFitnessCBNRcpp(1, genotypes, 0, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(2, c(0, 1, 0, 0, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")

wrapFitnessCBNRcpp(3, genotypes, 1, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(4, c(0, 0, 0, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")

wrapFitnessCBNRcpp(3, genotypes, 1, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(4, c(0, 0, 0, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")

wrapFitnessCBNRcpp(3, genotypes, 2, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(4, c(0, 0, 1, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")

wrapFitnessCBNRcpp(3, genotypes, 2, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(4, c(0, 0, 1, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")

wrapFitnessCBNRcpp(3, genotypes, 3, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(4, c(0, 1, 1, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")

wrapFitnessCBNRcpp(3, genotypes, 3, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(4, c(0, 1, 1, 1, 0, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")


wrapFitnessCBNRcpp(4, genotypes, 4, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(5, c(0, 1, 0, 1, 1, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")

wrapFitnessCBNRcpp(4, genotypes, 4, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(5, c(0, 1, 0, 1, 1, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")

wrapFitnessCBNRcpp(4, genotypes, 5, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(5, c(0, 0, 1, 1, 1, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")

wrapFitnessCBNRcpp(4, genotypes, 5, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(5, c(0, 0, 1, 1, 1, 0), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")


wrapFitnessCBNRcpp(4, genotypes, 6, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(5, c(0, 0, 1, 1, 1, 1), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")
fitness.CBN2(5, genotypes[, 7], restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")


wrapFitnessCBNRcpp(4, genotypes, 7, restrictTableC, 5,
               0.03, 0.25, 0.1, "Multiple", -99)
fitness.CBN2(5, c(1, 0, 1, 1, 1, 1), restrictTable, 5,
             0.03, 0.25, 0.1, "Multiple")

wrapFitnessCBNRcpp(4, genotypes, 7, restrictTableC, 5,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(5, c(1, 0, 1, 1, 1, 1), restrictTable, 5,
             0.03, 0.25, 0.1, "CBN")






## check it breaks
restrictTable <- rbind(
  c(1, 2, -9, -9),
  c(2, 0, -9, -9),
  c(3, 1, 1, -9),
  c(4, 1, 3, -9),
  c(5, 2, 4, 2)
  )

## problem: mutation 0 is said to depend on 2,
## and their indices are -10. Will lead to trouble with indexing.

(restrictTablepocha <- convertRestrictTable(restrictTable))
genotypes <- t(
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(0, 0, 0, 1, 0, 0),
    c(0, 0, 1, 1, 0, 0),
    c(0, 1, 1, 1, 0, 0),
    c(0, 1, 0, 1, 1, 0),
    c(0, 0, 1, 1, 1, 0),
    c(0, 0, 1, 1, 1, 1),
    c(1, 0, 1, 1, 1, 1))
  )

## OK, error
wrapFitnessCBNRcpp(3, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNRcpp(2, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNRcpp(0, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)

## OK, no error
wrapFitnessCBNRcpp(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNRcpp(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)



fitness.CBN2(2, genotypes[, 1], restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)

## does not break, but it would be nice if it did.
wrapFitnessCBNRcpp(0, genotypes, 7, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
## 0.03 + 4 * 0.25
fitness.CBN2(1, genotypes[, 8], restrictTable, 5,
               0.03, 0.25, 0.1, "CBN", -99)
## 0.03 + 5 * 0.25
fitness.CBN2(1, genotypes[, 8], restrictTable, 6,
               0.03, 0.25, 0.1, "CBN", -99)
fitness.CBN2(1, genotypes[, 8], restrictTable, 3,
               0.03, 0.25, 0.1, "CBN", -99)



wrapFitnessCBNstd(3, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNstd(2, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNstd(0, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)

## OK, no error
wrapFitnessCBNstd(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNstd(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
## does not break, but it would be nice if it did.
wrapFitnessCBNstd(0, genotypes, 7, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)




wrapFitnessCBNArma(3, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNArma(2, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNArma(0, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)

## OK, no error
wrapFitnessCBNArma(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
wrapFitnessCBNArma(1, genotypes, 0, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)
## does not break, but it would be nice if it did.
## it does seem to break with Arma
wrapFitnessCBNArma(0, genotypes, 7, restrictTablepocha, 5,
               0.03, 0.25, 0.1, "CBN", -99)











## to test thoroughly
generate.genotypes <- function(num = 50, genes = 30, p = 0.3) {
  matrix(sample(c(0L, 1L), size = num * genes,
                prob = c(1-p, p), replace = TRUE),
         ncol = num)
}

mutateThePos <- function(genotypes) {
  mut.pos <- sample(1:nrow(genotypes), size = ncol(genotypes),
                    replace = TRUE)
  for(i in 1:ncol(genotypes))
    genotypes[mut.pos[i], i] <- 1
  return(list(genot = genotypes, mutatedPos = mut.pos))
}


g1 <- generate.genotypes(num = 5000)
g.and.pos <- mutateThePos(g1)


call.both.fit1 <- function(restrictTable, lista, CBN.multiple = "CBN", 
                          birth.rate = 0.4, s = 0.017, fitness.parent = 0.003,
                          this.num = NULL) {

  rtC <- convertRestrictTable(restrictTable)
  

  drivers <- nrow(restrictTable)
  
  genotypes <- lista$genot
  poss <- lista$mutatedPos
  num <- length(poss)
  rr <- matrix(-999, nrow = num, ncol = 2)
#  if(is.null(this.num)) {
    for(i in 1:num) {
      genot <- genotypes[, i]
      rr[i, 1] <- wrapFitnessCBNRcpp(poss[i] - 1, genotypes, i - 1, rtC, drivers,
                                 birth.rate, s, fitness.parent, CBN.multiple, -88)
      rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
                               birth.rate, s, fitness.parent, CBN.multiple)
      if(rr[i, 1] != rr[i, 2])
        browser()
      
    }
    print(summary(rr[, 1] - rr[, 2]))
    return(rr)
  ## } else {
  ##   i <- this.num
  ##   genot <- genotypes[, i]
  ##   rr[i, 1] <- wrapFitnessCBN(poss[i] - 1, genotypes, i - 1, rtC, drivers,
  ##                              birth.rate, s, fitness.parent, CBN.multiple, -88)
  ##   rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
  ##                            birth.rate, s, fitness.parent, CBN.multiple)
  ##   browser()
  ##   browser()
  ## }
}


call.both.fit2 <- function(restrictTable, lista, CBN.multiple = "CBN", 
                          birth.rate = 0.4, s = 0.017, fitness.parent = 0.003,
                          this.num = NULL) {

  rtC <- convertRestrictTable(restrictTable)
  

  drivers <- nrow(restrictTable)
  
  genotypes <- lista$genot
  poss <- lista$mutatedPos
  num <- length(poss)
  rr <- matrix(-999, nrow = num, ncol = 2)
#  if(is.null(this.num)) {
    for(i in 1:num) {
      genot <- genotypes[, i]
      rr[i, 1] <- wrapFitnessCBNArma(poss[i] - 1, genotypes, i - 1, rtC, drivers,
                                 birth.rate, s, fitness.parent, CBN.multiple, -88)
      rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
                               birth.rate, s, fitness.parent, CBN.multiple)
      if(rr[i, 1] != rr[i, 2])
        browser()
      
    }
    print(summary(rr[, 1] - rr[, 2]))
    return(rr)
  ## } else {
  ##   i <- this.num
  ##   genot <- genotypes[, i]
  ##   rr[i, 1] <- wrapFitnessCBN(poss[i] - 1, genotypes, i - 1, rtC, drivers,
  ##                              birth.rate, s, fitness.parent, CBN.multiple, -88)
  ##   rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
  ##                            birth.rate, s, fitness.parent, CBN.multiple)
  ##   browser()
  ##   browser()
  ## }
}




call.both.fit3 <- function(restrictTable, lista, CBN.multiple = "CBN", 
                          birth.rate = 0.4, s = 0.017, fitness.parent = 0.003,
                          this.num = NULL) {

  rtC <- convertRestrictTable(restrictTable)
  

  drivers <- nrow(restrictTable)
  
  genotypes <- lista$genot
  poss <- lista$mutatedPos
  num <- length(poss)
  rr <- matrix(-999, nrow = num, ncol = 2)
#  if(is.null(this.num)) {
    for(i in 1:num) {
      genot <- genotypes[, i]
      rr[i, 1] <- wrapFitnessCBNstd(poss[i] - 1, genotypes, i - 1, rtC, drivers,
                                 birth.rate, s, fitness.parent, CBN.multiple, -88)
      rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
                               birth.rate, s, fitness.parent, CBN.multiple)
      if(rr[i, 1] != rr[i, 2])
        browser()
      
    }
    print(summary(rr[, 1] - rr[, 2]))
    return(rr)
  ## } else {
  ##   i <- this.num
  ##   genot <- genotypes[, i]
  ##   rr[i, 1] <- wrapFitnessCBN(poss[i] - 1, genotypes, i - 1, rtC, drivers,
  ##                              birth.rate, s, fitness.parent, CBN.multiple, -88)
  ##   rr[i, 2] <- fitness.CBN2(poss[i], genot, restrictTable, drivers,
  ##                            birth.rate, s, fitness.parent, CBN.multiple)
  ##   browser()
  ##   browser()
  ## }
}







restrictTable2 <- matrix(as.integer(c(
    1,    0,   -9,   -9,   -9,
    2,    0,   -9,   -9,   -9,
    3,    0,    1,   -9,   -9,
    4,    1,    3,   -9,   -9,
    5,    2,    4,    2,   -9,
    6,    1,    4,   -9,   -9,
    7,    2,    5,    6,   -9,
    8,    2,    3,    7,   -9,
    9,    2,    1,    2,   -9,
   10,    3,    1,    2,    4,
   11,    3,    2,    3,    5,
   12,    3,    1,    4,    5)), ncol = 5, byrow = TRUE)



rtC <- convertRestrictTable(restrictTable2)

tmp <- call.both.fit1(restrictTable2, g.and.pos)
tmpA <- call.both.fit2(restrictTable2, g.and.pos)
tmpB <- call.both.fit2(restrictTable2, g.and.pos)


tmp2 <- call.both.fit1(restrictTable2, g.and.pos, "Multiple")
tmp2A <- call.both.fit2(restrictTable2, g.and.pos, "Multiple")
tmp2B <- call.both.fit3(restrictTable2, g.and.pos, "Multiple")



gg <- g.and.pos[[1]][, 1]
pp <- 11
drivers <- 10
br <- 0.03
s <- 0.02
fp <- 0.1
CBN.multiple <- "Multiple"

ggg <- g.and.pos[[1]][,1:20]
wrapFitnessCBN(pp - 1, ggg, 0, rtC, drivers,
               br, s, fp, CBN.multiple, -88)
fitness.CBN2(pp, ggg[, 1], restrictTable2, drivers,
               br, s, fp, CBN.multiple)


## it is all the large matrix being passed
## but, for real, the table will already be there.

benchmark(wrapFitnessCBN(pp - 1, ggg, 0, rtC, drivers,
               br, s, fp, CBN.multiple, -88),
          fitness.CBN2(pp, ggg[, 1], restrictTable2, drivers,
                       br, s, fp, CBN.multiple),
          columns=c("test", "replications", "elapsed", "relative"),
          order="relative",
          replications=500)


##########################################################


### Algorithm 5,


### First version, without passengers.


### Number of different species: at most, number of events.


### What follows is all within-individual.
### Mixtures of CBNs can be thought of as among-individual
### as they induce a different covariation between drivers and passengers.


### "Classic CBN":
###      - lack of any dependence drives fitness to fitnessNo;
###      - thus: pressence of all deps.: fitnessYes.
### "Multiple routes":
###      - pressence of at least one dependence: fitnessYes
###      - lack of all deps: fitnessNo

### And a model with both? More complicated; later (if ever ;-)


### How to include Hjelm approach?
### Classic CBN: fitnessNo is not 0 (but less than fitnessYes)
### Multiple routes: fitnessNo is not 0 (but less than fitnessYes)


### For now, a single function, that counts. No early exit, but less
### branching and thus options for loop unrolling.


### First idea: two tables: Restrictions and indexRestrictions
### Restrictions:
###   - first column mutation in question
###   - second column

### IndexRestrictions: like a hash. Where restrictions are for each mut.
## restrict.table <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(3, 4),
##   c(1, 5),
##   c(5, 6),
##   c(6, 7),
##   c(6, 8)
##   )
## restrict.table <- restrict.table[, c(2, 1)]


### Second idea: single table.
### Columns:
###   - Mutation
###   - Number of deps
###   - All deps
###   - Rest of columns are -9
###
### Should be a ragged array, but simpler if fixed size, determined
###  at run-time.
###

### Will need to think how to go from general tree to this table,
### via adjacency matrix, formed with "isParent"?



fitness.CBN <- function(mutatedPos, genotype, restrict.table, num.drivers,
                        birth.rate, s, 
                        fitnessNo = 0, fitnessYes = fitness.linear) {
  ## I am no longer using this function
  ## this is loop format, for C

  if(mutatedPos > num.drivers) { # a passenger
    return(fitnessYes(genotype, birth.rate, s, num.drivers))
  } else {
    notMet <- FALSE
    
    for(i in 1:nrow(restrict.table)) {
      if(restrict.table[i] == mutatedPos) {
        if(genotype[restrict.table[i]] != 1) {
          notMet <- TRUE
          break
        }
      }
    }
    
    if(notMet) {
      return(fitnessNo)
    } else {
      return(fitnessYes(genotype, birth.rate, s, num.drivers))
    }
  }
}



## FIXME: can we use a log-shaped?
fitness.linear0 <- function(newGenotype, basal.birth.rate, s = 0.005, num.drivers = 20) {
  ## later, we will use num.drivers. For now not
  ## (1 + s)^(sum of mutations) as in beerenw, but need to have 1 be different.
  ## and this is exponential. Use a simple linear model
  return(basal.birth.rate + s * sum(newGenotype))
}







fitness.log <- function(newGenotype, basal.birth.rate, s = 0.005, num.drivers = 20) {
  return( basal.birth.rate + log(1 + sum(newGenotype)) )
}
  
## fitness <- fitness.linear


## Bozic et al. use stagnation prob:
## dk = 0.5 * (1 - s) ^k   (k num drivers). division = 1 - dk; division + stag = 1
## so at each step, either divide or stagnate
## plot(function(x, s = 0.04) {1 - (0.5 * (1 - s)^x)}, from = 0, to = 20)
## mutation: 10-5



### Wait! this is too convoluted? Just sample based upong number of cells?
## This is what Sprouffske et al. do, but they do have just mutations in drivers.
## However, for use, unless there are driver mutations, pop growth will be slow or cero.
## Do that: use a fitness that does not increase unless drivers.

areWeDone0 <- function(popsizes, genotypes, num.drivers = 20) {
  ## for now, keep it simple. Later, differentiate between drivers and passengers
}

## Later, turn size.for.detection into a random variable
areWeDone <- function(popsizes, size.for.detection) {
  ifelse(sum(popsizes) > size.for.detection, 1, 0)
}






## set.seed(3); oo <- Algo5(death = 0.2, birth.rate = 0.201, InitSize = 100, mu = 1e-7)

## MaxSpecies probably has to be much larger
## Algo5(death = 0.5, birth.rate = 0.5, InitSize = 100, mu = 1e-4, detectionSize = 1e7, sample.every = 1, s = 0.4)$muts.by.time


## The following seem reasonable, but we do not get to many mutations
## Algo5(death = 0.5, birth.rate = 0.5, InitSize = 100, mu = 1e-4, detectionSize = 1e6, sample.every = 1, s = 0.4, NumGenes = 20)$muts.by.time

## Algo5(death = 0.5, birth.rate = 0.5, InitSize = 10, mu = 1e-4, detectionSize = 1e7, sample.every = 1, s = 0.4, NumGenes = 30)$muts.by.time

## Algo5(death = 0.5, birth.rate = 0.5, InitSize = 10, mu = 5e-5, detectionSize = 5e6, sample.every = .5, s = 0.4, NumGenes = 30)$muts.by.time

## Algo5(death = 0.5, birth.rate = 0.5, InitSize = 10, mu = 5e-5, detectionSize = 5e6, sample.every = 1.0, s = 0.4, NumGenes = 30)$muts.by.time

## tmp <- Algo5(death = 0.5, birth.rate = 0.5, InitSize = 10, mu = 5e-5, detectionSize = 5e6, sample.every = 1.0, s = 0.4, NumGenes = 30)


## FIXME: to do:
## - write in C
## - check their code
## - define a range of param. values: really now?
## - trigger checks based on PopSize, so I avoid getting stuck if small mut. times: really?
## - spell out what is it we are doing.
## - think models to use


## simple model 1-> 2 -> 3 -> 4; 1 -> 5 -> 6 -> 7; 6 -> 8;
##
restrict.table <- rbind(
  c(1, 2),
  c(2, 3),
  c(3, 4),
  c(1, 5),
  c(5, 6),
  c(6, 7),
  c(6, 8)
  )
restrict.table <- restrict.table[, c(2, 1)]

set.seed(2)
tmp <- Algo5(death = 0.5, birth.rate = 0.501,
             InitSize = 50, mu = 2e-5, detectionSize = 5e6,
             fitness.function = fitness.CBN,
             restrict.table = restrict.table,
             num.drivers = 8,
             sample.every = 2.0, s = 0.4, NumGenes = 30)


## For thinking what should be close in memory:

## W: with birth and mutation
## R: with birth and mutation
## ti: assigned to ti and called with R, W, birth, PopSizes (tis[Flag])
## NextMutationTime: ???
## Algo3: R, W, birth, PosSizes
## fitness: genotypes
## Algo2: R, W, birth, PosSizes

## So: one matrix for genotypes and another for rest of stats.
## Arrays as provided by R and RCpp?

## What do I want to return to R?
## Where/when/how do I sample the individuals?

Algo5 <- function(MaxSpecies = 5000, NumGenes = 10, InitSize = 50,
                  death = 0.1, birth.rate = 0.1, mu = 1e-6,
                  ## maxTime = 100000,
                  sample.every = 500,
                  detectionSize = 1e6,
                  fitness.function = fitness.CBN,
                  restrict.table = restrict.table,
                  s = 0.005,
                  num.drivers = 20) {
##  on.exit(browser())
  
  ## In C: use a struct for all data for a species? Or arrays?
  ## For Genotype
  ## In C or C++ think this through. Could use an array of bits or similar.
  ## http://stackoverflow.com/questions/2525310/how-to-define-and-work-with-an-array-of-bits-in-c
  ## And, most crucially, I can quickly identify if identical genotype by doing a == b.
  ## these are like numbers in binary base


  ## sanity check
  if(max(restrict.table[, 1]) != num.drivers)
    stop("max(restrict.table[, 1]) != num.drivers")

  
  Genotype <- matrix(0, nrow = MaxSpecies, ncol = NumGenes)
  ## Too general. We set death to same values, birth increases.
  ## Like linear model, 3.1.1 in Mather et al.
  ## and similar to Durrett et al., 2010
  
  ## Rates <- matrix(nrow = MaxSpecies, ncol = 3) ## growth, death, mu
  ## too general for now. Will probably use a constant death
  ## and will only change mutation with mutator phenotype, later.

  ## selective.advantage is the w of Beerenw. and the a_m of Mather
  
  birth <- vector("numeric", length = MaxSpecies)
  mutation <- vector("numeric", length = MaxSpecies)
  Flag <- vector(mode = "logical", length = MaxSpecies)
  Flag <- rep(FALSE, length = MaxSpecies)
  PopSizes <- vector(mode = "integer", length = MaxSpecies)
  TimeLastUpdate <- vector(mode = "numeric", length = MaxSpecies)
  Exists <- vector(mode = "logical", length = MaxSpecies)
  Exists <- c(TRUE, rep(FALSE, MaxSpecies-1))


  
  
  Ws <- vector(mode = "numeric", length = MaxSpecies)
  Rs <- vector(mode = "numeric", length = MaxSpecies)

  NextMutationTime <- tis <- rep(NA, length = MaxSpecies)
  ## this probably ain't needed but makes life simpler
  ## indices <- seq.int(MaxSpecies)
  
  ## Might want, certainly in C, to have the W, R, etc, preallocated
  
  ## 5.1. Initialize
  birth[1] <- birth.rate
  Flag[1] <- TRUE
  mutation[1] <- mu * NumGenes
  PopSizes[1] <- InitSize
  num.species <- 1

  Ws[1] <- W.f(death, birth[1], mutation[1])
  Rs[1] <- R.f(death, birth[1], mutation[1])

  CurrentTime <- 0
  iter <- 0
  TimeAllPopSample <- CurrentTime + sample.every
  
  TimeLastUpdate <- rep(0, MaxSpecies)

  ## I crucially assume that new mutations are placed in the next empty slot.
  ## And if a species becomes extinct, its hole is left in all the data structures.

  out.ns <- matrix(ncol = MaxSpecies + 2, 
                   nrow = MaxSpecies * 2)

  ## FIXME: dimensions of out.ns!!!
  ## First column is iter, second is time
  
  out.ns[1, ] <- c(0, 0, PopSizes)
  out.i <- 1

  ## to keep track; not used
  ## mut.time <- vector(mode = "numeric", length = MaxSpecies)

  ## FIXME: how are we bailing out? do we need to sample? check their C++ code.

  ## stopping:
  ## a) num.species can never get larger than number of iterations.
  ## b) if(sum(PopSizes) < 0) this got screwed up, so stop it. It is an error.
  ## set a) and b) as limits. If we gt over it, abort. Nothing to be saved,
  ## except to figure out errors.


  ## I need to impose my conditions: fraction of cells with more than a given number of mutations.
  ## This is computed at sampling.

##  while((CurrentTime < maxTime) && (sum(PopSizes) > 0) && (num.species < MaxSpecies)){
  simulsDone <- 0
  while(!simulsDone) {
    
    iter <- iter + 1
    num.species <- sum(Exists)


    ########### 5.2 

    #### Ws and Rs do NOT need to be recomputed!! Only computed for new individuals

    if(any(is.na(c(Rs[Flag], Ws[Flag], birth[Flag], PopSizes[Flag])))) {
      cat("\n NAs in Rs, Ws, ...")
      browser()
    }

    tis[Flag] <- mapply(ti, Rs[Flag], Ws[Flag], death, birth[Flag],
                  PopSizes[Flag]) 
    
    ##    indices1 <- indices[Flag]

    ## debugging and monitor. code.
    if(any(tis[Exists] < 0)) {
      cat("\n tis < 0 \n")
      browser()
    }
    ## tis identically 0 should never happen. but see ti,
    ## which sets it, to prevent negative values from unerflow and rounding
    if(any(tis[Exists] == 0)) cat("\n BEWARE some tis == 0 \n")
    ## if(iter > 100) { ## FIXME: this is just for debugging, I guess??
    ##   cat("\n")
    ##   stop()
    ## }

    ## a vector
    NextMutationTime[Flag] <- tis[Flag] + CurrentTime

    ## cat("\n iter ", iter, "\n tis\n", tis[Exists], "\nNextMutationTime\n", NextMutationTime[Exists], "\n")
    ## cat("\n iter ", iter)

    
    ## in C: beware of infs
        
    ## if I used the next, I'd screw up the sampling
    ## TimeLastUpdate[Flag] <- CurrentTime
    Flag[Flag] <- FALSE

    ## This is really part of 5.3 but we need it to decide if we need to
    ## sample the whole pop
    NextMutant <- which.min(NextMutationTime[Exists])
    ## min.ti <- NextMutationTime[NextMut]
    ## MinimumNextMutationTime <- min.ti + CurrentTime
    MinimumNextMutationTime <- NextMutationTime[NextMutant]


  
    ## Do we need to sample the population?
    if(TimeAllPopSample >= MinimumNextMutationTime) { ## we are not sampling
      ## thus, we do 5.3 to 5.7 or 5.8

      
      
      ## 5.3
      ## See before, where we have computed MinimumNextMutationTime
      ## and found out the corresponding class
      CurrentTime <- MinimumNextMutationTime
      ##CurrentTime <- min.ti + CurrentTime 

      ## mut.time[iter] <- CurrentTime
      ## this is very wrong! we want NextMut
      ## indexMut <- indices1[NextMut]

      ## 5.4
      MutantTimeSinceLastUpdate <- CurrentTime - TimeLastUpdate[NextMutant]

      ## debug
      if(MutantTimeSinceLastUpdate > sample.every) {
        cat("\n\n MutantTimeSinceLastUpdate > sample.every \n")
        browser()
      }
      
      ## Isn't this always the minimum tis?
      ## Nope! Only first iterations
      ## if(MutantTimeSinceLastUpdate != MinimumNextMutationTime) {
      ##   cat("\ MutantTImeSinceLastUpdate, MinimymNextMutationTime\n")
      ##   cat(MutantTimeSinceLastUpdate, MinimumNextMutationTime, "\n")
      ##   browser()
      ## }

      ## cat("\n MutantTimeSinceLastUpdate ", MutantTimeSinceLastUpdate, "\n")

      PopSizes[NextMutant] <- Algo3(PopSizes[NextMutant],
                                MutantTimeSinceLastUpdate,
                                Rs[NextMutant],
                                Ws[NextMutant], death,
                                birth[NextMutant])

      TimeLastUpdate[NextMutant] <- CurrentTime

      ## 5.5
      ## watch out! sample(x, n) where length(x) = 1 does not return that element.
      mutable.pos <- which(!Genotype[NextMutant, ])
      if(length(mutable.pos) > 1)
        mutatedPos <- sample(mutable.pos, 1)
      else if(length(mutable.pos) == 1)
        mutatedPos <- mutable.pos
      else
        stop("\n no places left for mutations. What should we do? How could this happen?")
      
      ## 5.6 By decree, impossible to have a second recorded mutation in
      ## the same gene.
      
      newGenotype <- Genotype[NextMutant, ]
      ## FIXME: keep fitness of parent here; it should be fitness[NextMutant]
      ## or actually birth.rate

      newGenotype[mutatedPos] <- 1
      
      j <- 1
      
      ## FIXME: triple-check this condition always works!!
      while((j <= num.species) && !(identical(Genotype[j, ], newGenotype))) j <- j + 1
      if(j > num.species) {
        ## A new species has appeared
        Genotype[j, ] <- newGenotype
        num.species <- num.species + 1
        if(num.species > MaxSpecies) {
          stop("\n !!!!!!!!! num.species > MaxSpecies \n")
        }

        PopSizes[j] <- 1
        birth[j] <- fitness.function(mutatedPos, newGenotype,
                                     restrict.table,
                                     num.drivers= num.drivers,
                                     birth.rate,
                                     s = s)

        ## birth[j] <- fitness.function(newGenotype, birth.rate,
        ##                              s = s, num.drivers= num.drivers)

        ## Use a more general mutation function. Ditto above for fitness.
        ## Mut.rate is mut.rate per gene, so mutation rate goes down, as we
        ## only interested in mutation rate in different genes and this is
        ## Poisson
        mutation[j] <- mu * (sum(!newGenotype))
        if(mutation[j] == 0) { ## signal the "will not mutate" condition
         Ws[j] <- Rs[j] <- -99 
        } else {
          Ws[j] <- W.f(death, birth[j], mutation[j])
          Rs[j] <- R.f(death, birth[j], mutation[j])
        }
        Exists[j] <- TRUE
        TimeLastUpdate[j] <- CurrentTime
      } else {
        ## A mutation to a pre-existing species
        TimeSinceLastUpdate <- CurrentTime - TimeLastUpdate[j]
        if (TimeSinceLastUpdate <= 0 ) {
          cat("\n TimeSinceLastUpdate <= 0     \n")
          browser()
        }
        ## eh, we missing the "Exists!!"" set it to TRUE!!!
        tmpNum <- Algo2(PopSizes[j], TimeSinceLastUpdate,
                        Rs[j],
                        Ws[j], death,
                        birth[j])
        PopSizes[j] <- tmpNum + 1
        TimeLastUpdate[j] <- CurrentTime
      }
      
      if(NextMutant == j) stop ("SEVERE PROBLEM: NextMutant == j")
      ## 5.7
      Flag[NextMutant] <- TRUE
      Flag[j] <- TRUE  

    } else { ### We are sampling
      cat("\n WE ARE SAMPLING at time", TimeAllPopSample, "\n")
      ## FIXME: OJO here!!!

      ## I add the "Pop" to distinguish from the scalar;
      ## in C, predefine the matrix and have same dims as TimeLastSam

      ## but it is not needed.
      PopTimeSinceLastUpdate <- TimeAllPopSample - TimeLastUpdate[Exists]
      if (any(PopTimeSinceLastUpdate <= 0) ) {
        cat("\n ******* some PopTime < 0\n")
        browser()}
      ## do only with existing pops
      ## which, I guess, could include at least one with current 0 pop size
      ## check OK FIXME: try to update a zero size population

      ## THIS IS WRONG!! SAMPLE ONLY UNFLAGGED!!! Step 9 in algorithm!!!
      ## PopSizes[Exists] <- mapply(Algo2, PopSizes[Exists],
      ##                            PopTimeSinceLastUpdate,
      ##                            Rs[Exists],
      ##                            Ws[Exists],
      ##                            death,
      ##                            birth[Exists])
      PopSizes[!Flag] <- mapply(Algo2, PopSizes[!Flag],
                                 PopTimeSinceLastUpdate,
                                 Rs[!Flag],
                                 Ws[!Flag],
                                 death,
                                 birth[!Flag])

      
      CurrentTime <- TimeAllPopSample
      TimeLastUpdate[!Flag] <- CurrentTime
      ## fix later size of out.ns FIXME
      out.i <- out.i + 1
      out.ns[out.i, ] <- c(iter, CurrentTime, PopSizes)
      ## so we go to 5.2 with all of them
      Flag[!Flag] <- TRUE
      
      TimeAllPopSample <- TimeAllPopSample + sample.every

      ## check the stopping condition
      ## simulsDone <- areWeDone(PopSizes[Exists], Genotype[Exists, ])

      ## we might not enter here for a looong number of iter
      simulsDone <- areWeDone(PopSizes[Exists], detectionSize)

      if(sum(PopSizes) <= 1) {
        simulsDone <- 1
        cat("\n WATCH OUT: PopSizes is 0!!!!")
      }

    }
  }
  ## before exit, to get idea of what happened
  ## browser()

  ## FIXME: OJO, need to separate Mutations in drivers
  ## from those in passengers


  ## NOTES for C++
  ## From here on can be done in R itself.
  ## Return two matrices: Genotypes and the rest of stuff

  ## following is wrong, as a species could have gone extinct,
  ## and we'd want to see it
  ## or does Exists just flags that something came into existence?
  ## Yeap. Exists is never turned to false.
  NumMutations <- apply(Genotype[Exists, ], 1, sum)
  ll <- list(NumSpecies = sum(Exists),
             TotalPopSize = sum(PopSizes[Exists]),
             PopSizes = PopSizes[Exists],
             FinalTime = CurrentTime,
             Fitness = birth[Exists],
             NumMutations = NumMutations
             )
  print(ll)

  ## following is wrong, as a species could have gone extinct,
  ## and we'd want to see it
  pops.by.time <- out.ns[seq.int(out.i), 1:(sum(Exists) + 2) ]

  muts.by.time <- cbind(pops.by.time[, c(1, 2)],
                    t(apply(pops.by.time[, -c(1, 2)], 1,
                            function(x) tapply(x, NumMutations, sum))))
  colnames(muts.by.time)[c(1, 2)] <- c("Iteration", "Time")
  return(list(pops.by.time = pops.by.time,
              muts.by.time = muts.by.time, 
              ll))
}




## An attempt (not finished) to use relative fitness (i.e., growth depends
## on others too).
## Why not?
## a) it is an approximation and I don't want to have to mess with it.
## b) Bozic et al., Durrett et al., etc, do not feel compelled to doing this.
##    Why should I?
## c) more parameters (e.g., a scale factor for growth rate of pop.)

Algo6 <- function(MaxSpecies = 1000, NumGenes = 10, InitSize = 50,
                  death = 0.1, birth.rate = 0.1, mu = 1e-6,
                  ## maxTime = 100000,
                  sample.every = 500,
                  detectionSize = 1e6,
                  fitness.function = fitness.linear,
                  s = 0.005,
                  num.drivers = 20,
                  rel.fitness = TRUE) {
  on.exit(browser())
  ## In C: use a struct for all data for a species? Or arrays?
  ## For Genotype
  ## In C or C++ think this through. Could use an array of bits or similar.
  ## http://stackoverflow.com/questions/2525310/how-to-define-and-work-with-an-array-of-bits-in-c
  ## And, most crucially, I can quickly identify if identical genotype by doing a == b.
  ## these are like numbers in binary base
  
  Genotype <- matrix(0, nrow = MaxSpecies, ncol = NumGenes)
  ## Too general. We set death to same values, birth increases.
  ## Like linear model, 3.1.1 in Mather et al.
  ## and similar to Durrett et al., 2010
  
  ## Rates <- matrix(nrow = MaxSpecies, ncol = 3) ## growth, death, mu
  ## too general for now. Will probably use a constant death
  ## and will only change mutation with mutator phenotype, later.

  ## selective.advantage is the w of Beerenw. and the a_m of Mather
  if(rel.fitness) selective.advantage <- vector("numeric", length = MaxSpecies)
  
  birth <- vector("numeric", length = MaxSpecies)
  mutation <- vector("numeric", length = MaxSpecies)
  Flag <- vector(mode = "logical", length = MaxSpecies)
  Flag <- rep(FALSE, length = MaxSpecies)
  PopSizes <- vector(mode = "integer", length = MaxSpecies)
  TimeLastUpdate <- vector(mode = "numeric", length = MaxSpecies)
  Exists <- vector(mode = "logical", length = MaxSpecies)
  Exists <- c(TRUE, rep(FALSE, MaxSpecies-1))

  
  
  Ws <- vector(mode = "numeric", length = MaxSpecies)
  Rs <- vector(mode = "numeric", length = MaxSpecies)

  NextMutationTime <- tis <- rep(NA, length = MaxSpecies)
  ## this probably ain't needed but makes life simpler
  indices <- seq.int(MaxSpecies)
  
  ## Might want, certainly in C, to have the W, R, etc, preallocated
  
  ## 5.1. Initialize
  birth[1] <- birth.rate
  if(rel.fitness) {
    selective.advantage[1] <- 1
    birth[1] <- 1 ## which makes birth.rate useless here?
  }
  browser()
  Flag[1] <- TRUE
  mutation[1] <- mu * NumGenes
  PopSizes[1] <- InitSize
  num.species <- 1

  Ws[1] <- W.f(death, birth[1], mutation[1])
  Rs[1] <- R.f(death, birth[1], mutation[1])

  CurrentTime <- 0
  iter <- 0
  TimeAllPopSample <- CurrentTime + sample.every
  
  TimeLastUpdate <- rep(0, MaxSpecies)

  ## I crucially assume that new mutations are placed in the next empty slot.
  ## And if a species becomes extinct, its hole is left in all the data structures.

  out.ns <- matrix(ncol = MaxSpecies + 2, 
                   nrow = MaxSpecies * 2)

  ## FIXME: dimensions of out.ns!!!
  ## First column is iter, second is time
  
  out.ns[1, ] <- c(0, 0, PopSizes)
  out.i <- 1

  ## to keep track
  mut.time <- vector(mode = "numeric", length = MaxSpecies)

  ## FIXME: how are we bailing out? do we need to sample? check their C++ code.

  ## stopping:
  ## a) num.species can never get larger than number of iterations.
  ## b) if(sum(PopSizes) < 0) this got screwed up, so stop it. It is an error.
  ## set a) and b) as limits. If we gt over it, abort. Nothing to be saved,
  ## except to figure out errors.


  ## I need to impose my conditions: fraction of cells with more than a given number of mutations.
  ## This is computed at sampling.

##  while((CurrentTime < maxTime) && (sum(PopSizes) > 0) && (num.species < MaxSpecies)){
  simulsDone <- 0
  while(!simulsDone) {
    
    iter <- iter + 1
    num.species <- sum(Exists)


    ########### 5.2 

    #### Ws and Rs do NOT need to be recomputed!! Only computed for new individuals

    if(any(is.na(c(Rs[Flag], Ws[Flag], birth[Flag], PopSizes[Flag])))) {
      cat("\n NAs in Rs, Ws, ...")
      browser()
    }

    tis[Flag] <- mapply(ti, Rs[Flag], Ws[Flag], death, birth[Flag],
                  PopSizes[Flag]) 
    
    ##    indices1 <- indices[Flag]

    ## debugging and monitor. code.
    if(any(tis[Exists] < 0)) {
      cat("\n tis < 0 \n")
      browser()
    }
    ## tis identically 0 should never happen. but see ti,
    ## which sets it, to prevent negative values from unerflow and rounding
    if(any(tis[Exists] == 0)) cat("\n BEWARE some tis == 0 \n")
    ## if(iter > 100) { ## FIXME: this is just for debugging, I guess??
    ##   cat("\n")
    ##   stop()
    ## }

    ## a vector
    NextMutationTime[Flag] <- tis[Flag] + CurrentTime

    ## cat("\n iter ", iter, "\n tis\n", tis[Exists], "\nNextMutationTime\n", NextMutationTime[Exists], "\n")
    ## cat("\n iter ", iter)

    
    ## in C: beware of infs
        
    ## if I used the next, I'd screw up the sampling
    ## TimeLastUpdate[Flag] <- CurrentTime
    Flag[Flag] <- FALSE

    ## This is really part of 5.3 but we need it to decide if we need to
    ## sample the whole pop
    NextMutant <- which.min(NextMutationTime[Exists])
    ## min.ti <- NextMutationTime[NextMut]
    ## MinimumNextMutationTime <- min.ti + CurrentTime
    MinimumNextMutationTime <- NextMutationTime[NextMutant]


  
    ## Do we need to sample the population?
    if(TimeAllPopSample >= MinimumNextMutationTime) { ## we are not sampling
      ## thus, we do 5.3 to 5.7 or 5.8

      
      
      ## 5.3
      ## See before, where we have computed MinimumNextMutationTime
      ## and found out the corresponding class
      CurrentTime <- MinimumNextMutationTime
      ##CurrentTime <- min.ti + CurrentTime 

      mut.time[iter] <- CurrentTime
      ## this is very wrong! we want NextMut
      ## indexMut <- indices1[NextMut]

      ## 5.4
      MutantTimeSinceLastUpdate <- CurrentTime - TimeLastUpdate[NextMutant]

      ## debug
      if(MutantTimeSinceLastUpdate > sample.every) {
        cat("\n\n MutantTimeSinceLastUpdate > sample.every \n")
        browser()
      }
      
      ## Isn't this always the minimum tis?
      ## Nope! Only first iterations
      ## if(MutantTimeSinceLastUpdate != MinimumNextMutationTime) {
      ##   cat("\ MutantTImeSinceLastUpdate, MinimymNextMutationTime\n")
      ##   cat(MutantTimeSinceLastUpdate, MinimumNextMutationTime, "\n")
      ##   browser()
      ## }

      ## cat("\n MutantTimeSinceLastUpdate ", MutantTimeSinceLastUpdate, "\n")

      browser()
      PopSizes[NextMutant] <- Algo3(PopSizes[NextMutant],
                                MutantTimeSinceLastUpdate,
                                Rs[NextMutant],
                                Ws[NextMutant], death,
                                birth[NextMutant])

      TimeLastUpdate[NextMutant] <- CurrentTime

      ## 5.5
      ## watch out! sample(x, n) where length(x) = 1 does not return that element.
      mutable.pos <- which(!Genotype[NextMutant, ])
      if(length(mutable.pos) > 1)
        mutatedPos <- sample(mutable.pos, 1)
      else if(length(mutable.pos) == 1)
        mutatedPos <- mutable.pos
      else
        stop("\n no places left for mutations. What should we do?")
      
      ## 5.6
      ## By decree, impossible to have a mutation in the same gene.
      
      newGenotype <- Genotype[NextMutant, ]
      newGenotype[mutatedPos] <- 1
      
      j <- 1
      
      ## FIXME: triple-check this condition always works!!
      while((j <= num.species) && !(identical(Genotype[j, ], newGenotype))) j <- j + 1
      if(j > num.species) {
        ## A new species has appeared
        Genotype[j, ] <- newGenotype
        num.species <- num.species + 1
        if(num.species > MaxSpecies) {
          stop("\n !!!!!!!!! num.species > MaxSpecies \n")
        }

        PopSizes[j] <- 1
        if(rel.fitness) {
          selective.advantage[j] <- fitness.function(newGenotype,
                                                     birth.rate, s = s,
                                                     num.drivers= num.drivers)
          average.fitness <- sum(birth[Exists] * PopSizes[Exists])
          birth[j] <- selective.advantage[j]/average.fitness
        } else {
          birth[j] <- fitness.function(newGenotype, birth.rate,
                                       s = s, num.drivers= num.drivers)
        }
        ## Use a more general mutation function. Ditto above for fitness.
        ## Mut.rate is mut.rate per gene, so mutation rate goes down, as we
        ## only interested in mutation rate in different genes and this is
        ## Poisson
        mutation[j] <- mu * (sum(!newGenotype))
        if(mutation[j] == 0) { ## signal the "will not mutate" condition
         Ws[j] <- Rs[j] <- -99 
        } else {
          Ws[j] <- W.f(death, birth[j], mutation[j])
          Rs[j] <- R.f(death, birth[j], mutation[j])
        }
        Exists[j] <- TRUE
        TimeLastUpdate[j] <- CurrentTime
      } else {
        ## A mutation to a pre-existing species
        TimeSinceLastUpdate <- CurrentTime - TimeLastUpdate[j]
        if (TimeSinceLastUpdate <= 0 ) {
          cat("\n TimeSinceLastUpdate <= 0     \n")
          browser()
        }
        
        tmpNum <- Algo2(PopSizes[j], TimeSinceLastUpdate,
                        Rs[j],
                        Ws[j], death,
                        birth[j])
        PopSizes[j] <- tmpNum + 1
        TimeLastUpdate[j] <- CurrentTime
      }
      
      if(NextMutant == j) stop ("SEVERE PROBLEM: NextMutant == j")
      ## 5.7
      Flag[NextMutant] <- TRUE
      Flag[j] <- TRUE  

    } else { ### We are sampling
      cat("\n WE ARE SAMPLING at time", TimeAllPopSample, "\n")
      ## FIXME: OJO here!!!

      ## I add the "Pop" to distinguish from the scalar;
      ## in C, predefine the matrix and have same dims as TimeLastSam
      PopTimeSinceLastUpdate <- TimeAllPopSample - TimeLastUpdate[Exists]
      if (any(PopTimeSinceLastUpdate <= 0) ) {
        cat("\n ******* some PopTime < 0\n")
        browser()}
      ## do only with existing pops
      ## which, I guess, could include at least one with current 0 pop size
      ## check OK FIXME: try to update a zero size population
      PopSizes[Exists] <- mapply(Algo2, PopSizes[Exists],
                                 PopTimeSinceLastUpdate,
                                 Rs[Exists],
                                 Ws[Exists],
                                 death,
                                 birth[Exists])
      if(rel.fitness) {
        average.fitness <- sum(birth[Exists] * PopSizes[Exists])
        birth[Exists] <- selective.advantage[Exists]/average.fitness
      }
      CurrentTime <- TimeAllPopSample
      TimeLastUpdate[Exists] <- CurrentTime
      ## fix later size of out.ns FIXME
      out.i <- out.i + 1
      out.ns[out.i, ] <- c(iter, CurrentTime, PopSizes)
      ## so we go to 5.2 with all of them
      Flag[Exists] <- TRUE
      
      TimeAllPopSample <- TimeAllPopSample + sample.every

      ## check the stopping condition
      ## simulsDone <- areWeDone(PopSizes[Exists], Genotype[Exists, ])

      ## we might not enter here for a looong number of iter
      simulsDone <- areWeDone(PopSizes[Exists], detectionSize)

      if(sum(PopSizes) <= 1) {
        simulsDone <- 1
        cat("\n WATCH OUT: PopSizes is 0!!!!")
      }

    }
  }
  ## before exist, to get idea of what happened
  ## browser()
  NumMutations <- apply(Genotype[Exists, ], 1, sum)
  ll <- list(NumSpecies = sum(Exists),
             TotalPopSize = sum(PopSizes[Exists]),
             PopSizes = PopSizes[Exists],
             FinalTime = CurrentTime,
             Fitness = birth[Exists],
             NumMutations = NumMutations
             )
  print(ll)

  pops.by.time <- out.ns[seq.int(out.i), 1:(sum(Exists) + 2) ]

  muts.by.time <- cbind(pops.by.time[, c(1, 2)],
                    t(apply(pops.by.time[, -c(1, 2)], 1,
                            function(x) tapply(x, NumMutations, sum))))
  return(list(pops.by.time = pops.by.time,
              muts.by.time = muts.by.time, 
              ll))
}






## the next one will not work. left as a partial step to fixing the sampling.
Algo5.v0 <- function(MaxSpecies = 1000, NumGenes = 10, InitSize = 50,
                   death = 0.01, birth.rate = 0.012, mu = 1e-6,
                   maxTime = 10000,
                   sample.every = 100) {
  ## In C: use a struct for all data for a species? Or arrays?
  ## For Genotype
  ## In C or C++ think this through. Could use an array of bits or similar.
  ## http://stackoverflow.com/questions/2525310/how-to-define-and-work-with-an-array-of-bits-in-c
  ## And, most crucially, I can quickly identify if identical genotype by doing a == b.
  ## these are like numbers in binary base
  
  Genotype <- matrix(0, nrow = MaxSpecies, ncol = NumGenes)
  ## Too general. We set death to same values, birth increases.
  ## Like linear model, 3.1.1 in Mather et al.
  ## and similar to Durrett et al., 2010
  
  ## Rates <- matrix(nrow = MaxSpecies, ncol = 3) ## growth, death, mu
  ## too general for now. Will probably use a constant death
  ## and will only change mutation with mutator phenotype, later.
  
  birth <- vector("numeric", length = MaxSpecies)
  mutation <- vector("numeric", length = MaxSpecies)
  Flag <- vector(mode = "logical", length = MaxSpecies)
  Flag <- rep(FALSE, length = MaxSpecies)
  PopSizes <- vector(mode = "integer", length = MaxSpecies)
  TimeLastSampling <- vector(mode = "numeric", length = MaxSpecies)
  Exists <- vector(mode = "logical", length = MaxSpecies)
  Exists <- c(TRUE, rep(FALSE, MaxSpecies-1))

  Ws <- vector(mode = "numeric", length = MaxSpecies)
  Rs <- vector(mode = "numeric", length = MaxSpecies)
  
  
  ## this probably ain't needed but makes life simpler
  ## not used, really
  indices <- seq.int(MaxSpecies)
  
  ## Might want, certainly in C, to have the W, R, etc, preallocated
  
  ## 5.1. Initialize
 
  birth[1] <- birth.rate
  Flag[1] <- TRUE
  mutation[1] <- mu * NumGenes
  PopSizes[1] <- InitSize
  num.species <- 1

  Ws[1] <- W.f(death, birth[1], mutation[1])
  Rs[1] <- R.f(death, birth[1], mutation[1])

  
  CurrentTime <- 0
  iter <- 0
  TimeAllPopSample <- CurrentTime
  
  TimeLastSampling <- rep(0, MaxSpecies)

  ## I crucially assume that new mutations are placed in the next empty slot.
  ## And if a species becomes extinct, its hole is left in all the data structures.

  out.ns <- matrix(ncol = MaxSpecies + 2, 
                   nrow = MaxSpecies * 2)

  ## FIXME: dimensions of out.ns!!!
  ## First column is iter, second is time
  
  out.ns[1, ] <- c(0, 0, PopSizes)
  out.i <- 1

  ## to keep track
  mut.time <- vector(mode = "numeric", length = MaxSpecies)
  
  
  while((CurrentTime < maxTime) && (sum(PopSizes) > 0) && (num.species < MaxSpecies)){
    iter <- iter + 1
    num.species <- sum(Exists)
    TimeAllPopSample <- TimeAllPopSample + sample.every

    ########### 5.2 and 5.3 together

    #### Ws and Rs do NOT need to be recomputed!! Only computed for new individuals

    if(any(is.na(c(Rs[Flag], Ws[Flag], birth[Flag], PopSizes[Flag])))) {
      cat("\n NAs in Rs, Ws, ...")
      browser()
    }
      
    tis <- mapply(ti, Rs[Flag], Ws[Flag], death, birth[Flag],
                  PopSizes[Flag])

    indices1 <- indices[Flag]
    
    if(any(tis < 0)) {
      cat("\n tis < 0 \n")
      browser()
    }
    if(any(tis == 0)) cat("\n BEWARE some tis = 0 \n")


    cat("\n iter ", iter, "\t")
    cat("tis ", tis, "\n")

    if(iter > 100) { ## FIXME: this is just for debugging, I guees??
      cat("\n")
      stop()
    }

    ##?? or is it with + tis?
    TimeLastSampling[Flag] <- CurrentTime
    Flag[Flag] <- FALSE


    ## Note: either one or two. Never three or more. 
    ## this is really part of 5.3 but we need it for below
    NextMut <- which.min(tis)
    min.ti <- tis[NextMut]

    
    ## Do we need to sample the population?
##    if(TimeAllPopSample < CurrentTime) {
    if(TimeAllPopSample < min.ti) {
      
      ## FIXME: OJO here!!!

      ## I add the "Pop" to distinguish from the scalar;
      ## in C, predefine the matrix and have same dims as TimeLastSam
      PopTimeSinceLastUpdate <- TimeAllPopSample - TimeLastSampling[Exists]
      if (any(PopTimeSinceLastUpdate <= 0) ) {cat("\n ******* some PopTime < 0\n"); browser()}
      browser()
      ## do only with existing pops
      ## which, I guess, could include at least one with current 0 pop size
      ## check OK FIXME: try to update a zero size population
      PopSizes[Exists] <- mapply(Algo2, PopSizes[Exists],
                                 PopTimeSinceLastUpdate,
                                 Rs[Exists],
                                 Ws[Exists],
                                 death,
                                 birth[Exists])

      
      ## fix later size of out.ns FIXME
      out.i <- out.i + 1
      out.ns[out.i, ] <- c(iter, CurrentTime, PopSizes)
      ## note that non-existing popsizes should remain at 0
      ##      MutantTimeSinceLastUpdate <- CurrentTime - TimeAllPopSample
      ## The idea when updating the pop, for the mutant, is:
      ## Doing Algo3 with t1 is the same as
      ## Algo2 with t and Algo3 with t1 - t (and new pop.size).
      ## But that is not correct.
      ## See split.algo below, for numerical examples.
      ## Problematic if the population goes extinct.

    } else {

      ## I think I need to place all that follows here, all steps 5.3 to
      ## 5.9 NO! the other way around: if no need to sample, do normal
      ## stuff. Otherwise, go straight to 5.9

      
     MutantTimeSinceLastUpdate <- CurrentTime - TimeLastSampling[indexMut] 
    }




    


    
    ## Note: either one or two. Never three or more. 
    ## this is moved before the check for whether we need pop sample
    ## NextMut <- which.min(tis)
    ## min.ti <- tis[NextMut]

    ## 5.3 (though see also above)
    ## CurrentTime <- tis[NextMut] + CurrentTime 
    CurrentTime <- min.ti + CurrentTime 
    ## MinimumNextMutationTime <- tis[NextMut] + CurrentTime
    ## CurrentTime <- MinimumNextMutationTime
    mut.time[iter] <- CurrentTime
  
    indexMut <- indices1[NextMut]
     

       
    

    
    ########### 5.4
    

    MutantTimeSinceLastUpdate <- CurrentTime - TimeLastSampling[indexMut]
    ## Isn't this always the minimum tis?
    ## Yes it is if no sampling before next mut. time:
    ## uncomment below if you want if no sampling.

    ## if(!all.equal(TimeSinceLastUpdate, tis[NextMut])) {
    ##   cat("\n ejje not identical\n")
    ##   browser()
    ## }

    PopSizes[indexMut] <- Algo3(PopSizes[indexMut],
                                MutantTimeSinceLastUpdate,
                                Rs[indexMut],
                                Ws[indexMut], death,
                                birth[indexMut])

    ########### 5.5
    
    mutatedPos <- sample(which(!Genotype[indexMut, ]), 1)
    ## mutatedPos <- sample(NumGenes, size = 1)

    ########### 5.6
    ## by decree, impossible to have a mutation in the same gene.
    ## Mut.rate is mut.rate per gene, so mutation rate goes down, as we
    ## only interested in mutation rate in different genes
    
  
    newGenotype <- Genotype[indexMut, ]
    newGenotype[mutatedPos] <- 1

    j <- 1


    ## FIXME: triple-check this condition always works!!
    while((j <= num.species) && !(identical(Genotype[j, ], newGenotype))) j <- j + 1
    if(j > num.species) {
      ## A new species has appeared
      Genotype[j, ] <- newGenotype
      num.species <- num.species + 1
      PopSizes[j] <- 1
      birth[j] <- fitness(newGenotype, birth.rate)
      ## Use a more general mutation function. Ditto above for fitness.
      mutation[j] <- mu * (sum(!newGenotype))
      Ws[j] <- W.f(death, birth[j], mutation[j])
      Rs[j] <- R.f(death, birth[j], mutation[j])
      Exists[j] <- TRUE
      
    } else {
      TimeSinceLastUpdate <- CurrentTime - TimeLastSampling[j]
      if (TimeSinceLastUpdate <= 0 ) {browser()}

      tmpNum <- Algo2(PopSizes[j], TimeSinceLastUpdate,
                      Rs[j],
                      Ws[j], death,
                      birth[j])
      PopSizes[j] <- tmpNum + 1
    }

    if(indexMut == j) stop ("indexMut == j")
    Flag[indexMut] <- TRUE
    Flag[j] <- TRUE  
        
      ## newSpecies <- TRUE
      ## ## not very R'ish
      ## for(j in 1:num.species) {
      ##   if(Genotype[j] == newGenotype) {
      ##     newSpecies <- FALSE
      ##     break()
      ##   }
      ## }
      ## if(newSpecies)



  }
  return(out.ns[seq.int(out.i), ])
}





## In Algo.5:
## If we want to sample by iteration, not time,
## we must be careful, because we cannot update the pop. of the "to mutate"
## with algorithm 2. Update with algo2 all except the one to mutate.
## So we would need this code below? After 5.3.

## if(!(iter %% sample.every) ) {
##   TimeSinceLastUpdate <- CurrentTime - TimeLastSampling
##   browser()
##   ## do only with existing pops
##   ## which, I guess, could include at least one with current 0 pop size
##   ## check OK FIXME
##   tmpNum <- mapply(Algo2, PopSizes, TimeSinceLastUpdate,
##                    Rs,
##                    Ws, death,
##                    birth)

##   ## fix later size of out.ns FIXME
##   out.i <- out.i + 1
##   out.ns[out.i, ] <- c(iter, CurrentTime,
##                        PopSizes)
## }







## R bug?
## save(file = "this.crashes.rbinom.RData", sizex, pex, pmx)
## rbinom(n = 1, size = sizex , prob = (1 - (pex/pmx) ))


## gs, e1
## checking eq. S7 to S12

## Cf <- function(R, t) cosh(R * t/2)
## Sf <- function(R, t) sinh(R * t/2)

## g.A <- function(W, t, R, gamma, growth, n, s) {
##   (((W * Sf(R, t) - R * Cf(R, t))* exp(s) - 2 * gamma * Sf(R, t))/
##     ( (2 * growth * Sf(R, t) * exp(s)) - W * Sf(R, t) - R * Cf(R, t)) )^n

## }



## g.B <- function(W, t, R, gamma, growth, n, s) {
##   pm <- pM(t, R, W, gamma, growth)
##   pe <- pE(pm, W, gamma, growth)
##   pb <- pB(pe, gamma, growth)
##   g1 <- (1 - pb)/(1 - pb * exp(s))

##   return(
##     ((pm - pe)*exp(s)*g1+pe)^n
    
##     )
## }

## ## yes, they look the same
## g.A(3, 1, 2, 0.1, 0.1, 10, 2) - g.B(3, 1, 2, 0.1, 0.1, 10, 2)

## g.A(3, 10, 2, 0.1, 0.1, 10, 2) - g.B(3, 10, 2, 0.1, 0.1, 10, 2)








### the difference between the above and these versiones of algo 2 and 3
### is whether we use the definition of the nb as in wikipedia (below)
### or in R. 


## Algo2_0 <- function(num, t, R, W, death, growth) {

##   ## cat("\n Arguments", num, t, R, gm, W, g, "\n")
##   pm <- pM(t, R, W, death, growth)
##   pe <- pE(pm, W, death, growth)
##   pb <- pB(pe, death, growth)

##   if(any(is.na(c(pm, pe, pb)))) browser()

## ##  browser()
##   if( (1 - pe/pm) > 1) browser()
##   if( (1 - pe/pm) < 0) browser()
##   if(pb > 1) browser()
##   if(pb <= 0) browser()
  
  
##   if(num < 0) browser()
##   m <- rbinom(n = 1, size = num , prob = (1 - (pe/pm) ))

##   if(m == 0) {
##     n <- 0
##   } else {
##     n <- m + rnbinom(n = 1, size = m, prob = pb)
##   }

## ##  cat("\n            Algo. 2.   m = ", m, ".  n = ", n)
##   return(n)
## }



## Algo3_0 <- function(num, t, R, W, death, growth) {
##   pm <- pM(t, R, W, death, growth)
##   pe <- pE(pm, W, death, growth)
##   pb <- pB(pe, death, growth)

##   if(any(is.na(c(pm, pe, pb)))) browser()
## ##  browser()
##   ## checks
##   if( (1 - pe/pm) > 1) browser()
##   if( (1 - pe/pm) < 0) browser()
##   if(pb > 1) browser()
##   if(pb <= 0) browser()
  
##   if(num < 1) browser()
##   ## cat(" Inside algo 3", "num ", num, "prob ", (1 - (pe/pm) ),
##   ##     "pe ", pe, "pm", pm)
##   ## pex <<- pe
##   ## pmx <<- pm
##   ## sizex <<- num - 1
##   m <- rbinom(n = 1, size = num - 1, prob = (1 - (pe/pm) ))

  
##   n <- m + 1 + rnbinom(n = 1, size = m + 2, prob = pb)
## ##  cat("\n          Algo. 3.   m = ", m, ".  n = ", n)
##   return(n)
## }





split.algo <- function(num, t1, t2, death = 0.01, growth = 0.011,
                       mutation = 1e-6) {
  ## Trying to understand the early stopping and what to do about updating. 
  ## Nope, this ain't fully correct.
  W <- W.f(death, growth, mutation)
  R <- R.f(death, growth, mutation)
  a <- Algo3(num, t1 + t2, R, W, death, growth)

  b1 <- Algo2(num, t1, R, W, death, growth)
  b2 <- ifelse(b1 == 0, 0, Algo3(b1, t2, R, W, death, growth))
  return(c(a, b1, b2, b2 - a, (b2 - a)/a))
}


## split.algo(5000, 1, 2000)

## summary(replicate(10000, split.algo(500, 1, 2000)[4]))
## summary(replicate(10000, split.algo(50, 1000, 20)[4]))







Moran.1 <- function(ns, mus, deaths, growths, max.events = 1000) {
  ## two classes
  ## no selection

  ## ds: what Mather call lambda. The death rate.
  ## gs: what Mather call g. Birth rate.
  ## mus: mutation rates.
  
  out.ns <- matrix(nrow = max.events, ncol = 4)
  
  t <- 0
  N <- sum(ns)
  num.classes <- length(ns)
  
  event <- 0
  while(event < max.events) {

#    browser()
    ## 4.2
    Ws <- mapply(W, deaths, growths, mus)
    Rs <- mapply(R, deaths, growths, mus)
    tis <- mapply(ti, Rs, Ws, deaths, growths, ns)

    if(any(tis < 0)) browser()
    if(any(tis == 0)) cat("\n BEWARE some tis = 0 \n")
    
    ## 4.3
    mutated <- which.min(tis)
    tm <- tis[mutated]
    if(is.infinite(tm)) {
      event <- event + 1
      out.ns[event, ] <- c(Inf, NA, NA, NA)
      break()
    }
    t <- t + tm

##     browser()
    
    ## For better R code, use mapply and sapply, and evaluate 4.5
    ## before 4.4?
    
    ## 4.4
    ns[mutated] <- Algo3(ns[mutated], tm, Rs[mutated],
                         Ws[mutated], deaths[mutated],
                         growths[mutated])
    ## 4.5
    ## as written here, this is too general for just two classes
    non.mutated <- seq_len(num.classes)[-mutated]
    for(i in non.mutated) {
      ns[i] <- Algo2(ns[i], tm, Rs[i], Ws[i], deaths[i], growths[i])
    }
##    browser()

    ## 4.6
    ## more work to do here for more complex scenarios
    if(mutated == 1) {
      ns[2] <- ns[2] + 1
    } else {
      ns[1] <- ns[1] + 1
    }

    ## cat(" t = ", t, ". Mutated was ", mutated,
    ##     ".  ns", paste(ns, collapse = ","), "\n")

    event <- event + 1
    out.ns[event, ] <- c(t, mutated, ns[1], ns[2])

    ## need to reset intensities
    ## gm <- 1/sum(ns)
    ## gs <- 1/sum(ns)

    
##    if( sum(ns) != N) stop(" Not equal N!!!!")
    ## if(any(ns == 0)) {
    ##   cat("\n bailing out at ns == 0\n ")
    ##   break()
    ## }
    if(sum(ns) > 10e15) break()
  }
  
  return(out.ns[1:event, ])
}


## debug(Moran.1)
## set.seed(1);
## oo1 <- Moran.1(ns = c(500, 500),
##                mus = c(.00001, .0001),
##                deaths = c(0.1, 0.1),
##                growths = c(0.2, 0.2),
##                max.events = 1000)


plot.traject.Moran <- function(matrix, n = nrow(matrix)) {
  matplot(x = matrix[1:n , 1],
          y = matrix[1:n, c(3, 4)],
          col = c("red", "blue"),
          type = "b")
}

## plot.traject(Moran.1(ns = c(500, 500),
##                mus = c(.001, .001),
##                deaths = c(0.1, 0.1),
##                growths = c(0.2, 0.2),
##                max.events = 1000))








require(inline)
testfun <- cxxfunction(
  signature(x = "numeric"),
  body = '
       NumericVector xx(x);
       return( xx );',
  plugin = "Rcpp"
  )

## testfun(1:4)




##############################

### Using inline, etc

## return as a new object
## W.f.c <- cxxfunction(signature(death = "numeric",
##                                growth = "numeric",
##                                mu = "numeric"),
##                      body = '
##                        NumericVector xd(death);
##                        NumericVector xg(growth);
##                        NumericVector xmu(mu);
##                        NumericVector W(1);
##                        W = xd + xg + xmu;
##                        return (W);
##                          ',
##                      plugin = "Rcpp"
##                      )

## return as a passed object
## W.f.c2 <- cxxfunction(signature(death = "numeric",
##                                 growth = "numeric",
##                                 mu = "numeric",
##                                 w = "numeric"),
##                       body = '
##                        NumericVector xd(death);
##                        NumericVector xg(growth);
##                        NumericVector xmu(mu);
##                        NumericVector W(w);
##                        W = xd + xg + xmu;
##                        return (W); // yes, return needed
##                          ',
##                       plugin = "Rcpp"
##                       )

## turning them into doubles, and using wrap
W.f.c2 <- cxxfunction(signature(death = "numeric",
                                growth = "numeric",
                                mu = "numeric",
                                w = "numeric"),
                      body = '
                       double xd = as<double>(death);
                       double xmu = as<double>(mu);
                       double xg = as<double>(growth);
                       double W = as<double>(w);
                       W = xd + xg + xmu;
                       return (wrap(W)); // yes, return needed
                         ',
                      plugin = "Rcpp"
                      )


R.f.c <- cxxfunction(signature(death = "numeric",
                               growth = "numeric",
                               mu = "numeric",
                               r = "numeric"),
                     body = '
                       double xd = as<double>(death);
                       double xmu = as<double>(mu);
                       double xg = as<double>(growth);
                       double R = as<double>(r);
                       R = sqrt( pow( xg - xd, 2) + ( 2 * (xg + xd) + xmu) * xmu  );
                       return (wrap(R));
                         ',
                     plugin = "Rcpp"
                     )


pM.f.c <- cxxfunction(
  signature( t = "numeric",
            R = "numeric",
            W = "numeric",
            death = "numeric",
            growth = "numeric",
            pM = "numeric"),
  body =
  '
           double xdeath = as<double>(death);
           double xgrowth = as<double>(growth);
           double xR = as<double>(R);
           double xW = as<double>(W);
           double xT = as<double>(t);
            // Maybe use long doubles, and then cast?
           double xpM = as<double>(pM);
           double Ct = cosh(xR * xT/2);
           double St = sinh(xR * xT/2);

           if( (!std::isfinite(Ct) ) || (!std::isfinite(St)) ) {
                 throw std::range_error("pM.f: Ct or St too big");
           }
           
           xpM = (xR * Ct + St * (2 * xdeath - xW ))/(xR * Ct + St * (xW - 2 * xgrowth));

           if( !std::isfinite(xpM) ) {
               throw std::range_error("pM.f: xpM not finite");
           }
           return(wrap(xpM));
',
  plugin = "Rcpp" , verbose = TRUE)


pE.f.c <- cxxfunction(
  signature(pM = "numeric",
            W = "numeric",
            death = "numeric",
            growth = "numeric",
            pE = "numeric"),
  body = '
double xpM = as<double>(pM);
double xW = as<double>(W);
double xdeath = as<double>(death);
double xgrowth = as<double>(growth);
double xpE = as<double>(pE);

xpE = (xdeath * (1.0 - xpM ) )/(xW - xdeath - xgrowth * xpM );
return(wrap(xpE));

',
  plugin = "Rcpp", verbose = TRUE)


pB.f.c <- cxxfunction(
  signature(pE = "numeric",
            death = "numeric",
            growth = "numeric",
            pB = "numeric"),
  body = '
double xpB = as<double>(pB);
double xdeath = as<double>(death);
double xgrowth = as<double>(growth);
double xpE = as<double>(pE);

xpB =  (xgrowth * xpE)/xdeath; 
return(wrap(xpB));
',
  plugin = "Rcpp", verbose = TRUE
  )

ti.c <- cxxfunction(
  signature(R = "numeric",
            W = "numeric",
            death = "numeric",
            growth = "numeric",
            n = "numeric",
            ti = "numeric"),
          include = '#include <limits>
                     #include <iostream>',
  body = '

double xdeath = as<double>(death);
double xgrowth = as<double>(growth);
double xn = as<double>(n);
double xR = as<double>(R);
double xW = as<double>(W);
double xti = as<double>(ti);
double eq11;
double r;
double rr;

// W < 0 is a signal that mutation is zero, and thus ti is Inf
if(xW <= -99.0) {
 xti = std::numeric_limits<double>::infinity();
 } else {

// using version with log
RNGScope scope;
r = ::Rf_runif(0.0, 1.0);
rr = exp((1/xn) * log(r));

std::cout << "r = " << r << " rr = " << rr << std::endl;

// eq.12
if( ((xR - xW + 2 * xdeath)/(xR + xW - 2 * xgrowth)) < rr ) {


std::cout << "numerator = " << (xR - xW + 2 * xgrowth) - xW - xR + 2 * xdeath << std::endl;
std::cout << "denominator = " << (-xR -xW + 2 * xgrowth) - xW + xR + 2 * xdeath << std::endl;


   //eq. 11
   xti = (1/xR) * log( (rr * (xR - xW + 2 * xgrowth) - xW - xR + 2 * xdeath) /
                     (rr * (-xR -xW + 2 * xgrowth) - xW + xR + 2 * xdeath));
   if(xti < 0.0) {
          throw std::range_error("ti: eq.11 < 0");
          }
   // return(wrap(xti));
   } else {
    xti = std::numeric_limits<double>::infinity();
 //return(wrap(xti));
  }
}
return(wrap(xti));
',
          plugin = "Rcpp" , verbose = TRUE)




Algo2.c <- cxxfunction(
             signature(
               num = "numeric",
               t = "numeric",
               R = "numeric",
               W = "numeric",
               death = "numeric",
               growth = "numeric"),
             plugin = "Rcpp",
             verbose = TRUE,
             body = '
// note we are modifying the num parameter

           double xdeath = as<double>(death);
           double xgrowth = as<double>(growth);
           double xR = as<double>(R);
           double xW = as<double>(W);
           double xT = as<double>(t);
           double xnum = as<double>(num);                      

           double pm;
double pe;
double pb;

// this seems not to work.It is the "." thing
pm = pM.f.c(xT, xR, xW, xdeath, xgrowth, pm);
'
       
             )


## for now, just CBN and Hjelm-like

## restrict.table
## mutatedPos, neededPrevious


fitness.CBN.c <- cxxfunction(
  signature(mutatedPos = "integer",
            genotype = "integer",
            restrictTable = "integer",
            numDrivers = "integer",
            birthRate = "numeric",
            s = "numeric",
            fitnessNo = "numeric",
            fitnessYes = "numeric",
            rvi = "integer"),
  body = '

int rv = as<int>(rvi);
IntegerMatrix restrictions(restrictTable);
int nr = restrictions.nrow();
int nc = restrictions.ncol();
// later row or col major order
rv = nr;

for(int i = 0; i < nr; i++) {
std::cout << "row = " << i << std::endl;
for(int j = 0; j < nc; j++) {
std::cout << "   " << restrictions(i, j) << std::endl;

}
}




// later with access to seq. pos 
return (wrap(rv));
',
  plugin = "Rcpp",
  verbose = TRUE)

m1 <- matrix(seq.int(6), ncol = 2, byrow = TRUE)

fitness.CBN.c(2, 3, m1, 2, 2, 2, 2, 2, 99)
