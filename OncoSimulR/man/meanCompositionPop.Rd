\name{meanCompositionPop}
\alias{meanCompositionPop}
\alias{stripChartPop}
\encoding{UTF-8}

\title{
  Plot stripChart to see the distribution of the mean results from several 
  iterations of oncoSimulPop since the half to the end of simulation.
}
\description{
  It allows to see the distribution of results form several simulations from
  oncosimulPop. It generates a strip chart with the average of the results from 
  the middle to the end of a simulation for each genotype. Lines connect
  data from the same simulation.
}
\usage{
  meanCompositionPop(objPop, ylab = "N", ...)
}

\arguments{
  \item{objPop}{
Data frame from the class oncosimulpop or oncosimul which contains the results
of the simulations in \code{pop.by.time}.
}
  \item{ylab}{
    Text for y axis label.
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
\keyword{hplot}
