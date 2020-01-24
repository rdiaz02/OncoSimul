\name{meanCompositionPop}
\alias{meanCompositionPop}
\alias{stripChartPop}
\encoding{UTF-8}

\title{
  Plot stripChart to see the distribution of the mean results from several 
  iterations of oncoSimulPop for the last fraction (fract) of the simulations.
}

\description{

  Visualize the distribution of results form several simulations from
  oncosimulPop. It generates a strip chart with the average of the
  results from the last fraction of the simulations (i.e., the average
  from all samples from total number of samples * fract until the end of
  the simulation) for each genotype. Lines connect data from the same
  simulation.
  
}

\usage{
  meanCompositionPop(objPop, ylab = "N", fract = 0.5, ...)
}

\arguments{
  \item{objPop}{
Data frame from the class oncosimulpop or oncosimul which contains the results
of the simulations in \code{pop.by.time}.
}
  \item{ylab}{
    Text for y axis label.
  }

  \item{fract}{ What fraction of the simulations to use. 
    The simulations used to compute the average are those from total
    number of simulations * fract until the last simulation.
  }
  
  \item{...}{
    Other arguments passed to plot.
  }
}

\author{
  Sara Dorado Alfaro,
  Miguel Hernández del Valle,
  Alvaro Huertas García,
  Diego Mañanes Cayero,
  Alejandro Martín Muñoz.
}

\seealso{
  \code{\link{oncoSimulIndiv}}, \code{\link{compositionPop2}},
  \code{\link{plot.oncosimul}}
}

\examples{


sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                             "Root" = "Root",
                             "A" = "a1, a2",
                             "B" = "b",
                             "C" = "c"))
evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)

e1 <- oncoSimulPop(5, sv2, model = "McFL",
                     mu = 5e-6,
                     sampleEvery = 0.02,
                     keepEvery = 1,
                     initSize = 5000,
                   finalTime = 1000,
                   mc.cores = 2,
                   onlyCancer = FALSE)


meanCompositionPop(e1)

## The next example will only work if
## we are using frequency-dependent fitness
## See https://github.com/rdiaz02/OncoSimul/tree/freq-dep-fitness

if("frequencyDependentFitness" %in% names(formals(allFitnessEffects))) {

## Map fitness
crs <- function (a, b, c){
  data.frame(Genotype = c("WT", "C", "R"),
             Fitness = c(paste0("1 + ", a, " * f_2 - ", b, " * f_1"),
                         paste0("1 + ", b, " * f_ - ", c, " * f_2"),
                         paste0("1 + ", c, " * f_1 - ", a, " * f_")
             ))
}

## Set frequency-dependent fitness
afcrs2 <- allFitnessEffects(genotFitness = crs(10, 1, 1), 
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

## Run simulations
resultscrs2 <- oncoSimulPop(10,
                           afcrs2,
                             model = "McFL", 
                             onlyCancer = FALSE,
                             finalTime = 100, 
                             mu = 1e-2,
                             initSize = 4000, 
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)

## Create strip chart plot
meanCompositionPop(resultscrs2, ylab = "Number of cells", cex.lab=1.1,
     las = 1)
}
}
\keyword{hplot}
