\name{compositionPop2}
\alias{compositionPop2}
\alias{simul_boxplot2}
\encoding{UTF-8}

\title{
  Plot a boxplot to see the distribution of the results from several 
  runs of oncoSimulPop.
}
\description{
  Extracts the last result for each genotype from different simulations of 
  oncoSimulPop to create a data frame which is used to show a boxplot
  (using ggplot2) to see the distribution of the results from several simulations. 
  The mean is represented by a red point on boxplots.
}
\usage{
compositionPop2(objPop, cols = NULL,  xlab = "Genotype",
                ylab = "N", main = FALSE)
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

compositionPop2(e1)
                
## The next example will only work if
## we are using frequency-dependent fitness
## See https://github.com/rdiaz02/OncoSimul/tree/freq-dep-fitness

if("frequencyDependentFitness" %in% names(formals(allFitnessEffects))) {

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
                   ),
                   stringsAsFactors = FALSE)
    }

    ## Show expressions for birth rates
    S1S2_LV("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")
    
    S1S2_competition <-allFitnessEffects(
        genotFitness = S1S2_LV(1.2, 1.7, 350, 400, 0, 0, ## a_12, a_21 are 0
                               gt = c("WT", "S1", "S2")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs", ## frequency expressed in absolute terms
        spPopSizes = c(50, 50, 50)) ## Same population, lower than K
    
    S1S2_simulation <- oncoSimulPop(10, S1S2_competition,
                                    model = "Exp",
                                    mc.cores = 2,
                                    onlyCancer = FALSE,
                                    finalTime = 100,
                                    mu = 1e-1, ## Allows appearance of S1 and S2
                                    initSize = 4000,
                                    seed = NULL,
                                    errorHitMaxTries = FALSE,
                                    errorHitWallTime = FALSE)
    
    
    ## Plot boxplots                          
    compositionPop2(S1S2_simulation, main = "No competition effect", 
                    ylab = "Number of individuals")
    ## Compare with
    meanCompositionPop(S1S2_simulation)
}

}
\keyword{hplot}

