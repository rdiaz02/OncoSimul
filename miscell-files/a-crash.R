##
o1 <- allFitnessEffects(
    noIntGenes = c("a" = .1, "b" = 0.2, "c" = 0.3))
set.seed(1)

oncoSimulIndiv(o1, initMutant = "b, a, c", verbosity = 1)
oncoSimulIndiv(o1, initMutant = "b, a, c", verbosity = 3, mutationPropGrowth = FALSE)

oncoSimulIndiv(o1, initMutant = "b, a, c, a", verbosity = 3)

oncoSimulIndiv(o1, initMutant = "b, a, c, d", verbosity = 3)




## so test correctness of the following and no crashes as lead to stop in R

evalAllGenotypes(o1, order = TRUE, addwt = TRUE)
## this crashes R
oncoSimulIndiv(o1, initMutant = "a, b")

oncoSimulIndiv(o1, initMutant = "a > b")

oncoSimulIndiv(o1, initMutant = "b, a")

## segfault
oncoSimulIndiv(o1, initMutant = "b, a, c")



  ## c > a > b: (1 + 0.9) * (1 - 0.2)
  ## a > c > b: (1 + 0.9) * (1 - 0.2)
  ## c > b > a: (1 + 0.9) * (1 + 0.3)
  ## a > c = c > a:   (1 + .9)





## Some things to include in tests and things that crash R

## This shows same gene in orderEffects and epistasis
o99 <- allFitnessEffects(
    orderEffects = c("a>b" = -0.2, "b > a" = 0.3),
    epistasis = c("a:c" = 0.9))
    
  evalAllGenotypes(o99, order = TRUE, addwt = TRUE)

  ## c > a > b: (1 + 0.9) * (1 - 0.2)
  ## a > c > b: (1 + 0.9) * (1 - 0.2)
  ## c > b > a: (1 + 0.9) * (1 + 0.3)
  ## a > c = c > a:   (1 + .9)

  ## this crashes R
  oncoSimulIndiv(o99, initMutant = "a>b>c")

o999 <- allFitnessEffects(
    orderEffects = c("a>b" = -0.2, "b > a" = 0.3),
    epistasis = c("a:c" = 0.9, "b:e" = 0.1))
    
  evalAllGenotypes(o999, order = TRUE, addwt = TRUE)

## this crashes R
oncoSimulIndiv(o999, initMutant = "a>c>b>e")



## but so does this
o999 <- allFitnessEffects(
    ## orderEffects = c("a>b" = -0.2, "b > a" = 0.3),
    epistasis = c("a:c" = 0.9, "b:e" = 0.1))
    
  evalAllGenotypes(o999, order = TRUE, addwt = TRUE)

oncoSimulIndiv(o99, initMutant = "a>b>c>e")
oncoSimulIndiv(o99, initMutant = "a,b,c,e")


