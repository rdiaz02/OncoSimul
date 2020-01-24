\name{compositionPop2}
\alias{compositionPop2}
\alias{simul_boxplot2}
\encoding{UTF-8}

\title{
  Plot a boxplot to see the distribution of the results from several 
  iterations of oncoSimulPop.
}
\description{
  Extracts the last result for each genotype from different simulations of 
  oncoSimulPop to create a data frame which is used to show a boxplot
  (using ggplot2) to see the distribution of the results from several simulations. 
  The mean is represented by a red point on boxplots.
}
\usage{
compositionPop2(objPop, cols = NULL,  xlab = "Genotype",
                ylab = "N", main = FALSE, ...)
}
\arguments{
  \item{objPop}{
    Data frame from the class \code{oncosimulpop} or \code{oncosimul} which
    contains the results of the simulations in \code{pop.by.time}.
}
  \item{cols}{
    Colours used to represent the boxplot. If \code{cols == NULL} (by default),
    colours from stacked plots of OncoSimulR are used. 
  }
  \item{main}{
    Optional title for the boxplot. If \code{main == FALSE} (by default), 
    no title is shown.If main is different from FALSE the boxplot title will 
    be the string passed in main. 
  }
  \item{xlab}{
    Text for x axis label.
  }
  \item{ylab}{
    Text for y axis label.
  }
}

\author{
  Sara Dorado Alfaro,
  Miguel Hernández del Valle,
  Alvaro Huertas García,
  Diego Mañanes Cayero,
  Alejandro Martín Muñoz
}

\seealso{
  \code{\link{oncoSimulIndiv}}, \code{\link{plot.oncosimul}},
  \code{\link{meanCompositionPop}}
}

\examples{
S1S2_LV <- function(r1, r2, K1, K2, a_12, a_21, awt = 1,
                    gt = c("WT", "S1", "S2")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("max(0.1, 1 - ", awt, " * (n_2 + n_1))"),
               paste0("1 + ", r1,
                      " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                      ")"),
               paste0("1 + ", r2,
                      " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                      ")")
             ))
}

## Show expressions for birth rates
S1S2_LV("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")

S1S2_competition <-allFitnessEffects(
    genotFitness = S1S2_LV(1.2, 1.7, 350, 400, 0, 0, ## a_12, a_21 are 0
              gt = c("WT", "S1", "S2")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs", ## frequency expressed in absolute terms
    spPopSizes = c(50, 50, 50)) ## Same population, lower than K

## Plot fitness landscape of genotype "S1, S2" evaluation
data.frame("S1_fitness" = evalGenotype(genotype = "S1", 
                                       fitnessEffects = S1S2_competition), 
           "S2_fitness" = evalGenotype(genotype = "S2", 
                                       fitnessEffects = S1S2_competition))

S1S2_simulation <- oncoSimulPop(20, S1S2_competition, ## 10 simulations
                              model = "Exp",
                              mc.cores = n_cores,
                              onlyCancer = FALSE,
                              finalTime = 500,
                              mu = 1e-1, ## Allows appearance of S1 and S2
                              initSize = 4000,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
                              

## Plot boxplots                          
compositionPop2(S1S2_simulation, main = "No competition effect", 
                ylab = "Number of individuals")

}
\keyword{hplot}

