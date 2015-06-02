\name{samplePop}
\alias{samplePop}

\title{
  Obtain a sample from a population of simulations.
  
}
\description{
  
  Obtain a sample (a matrix of individuals by genes or, equivalently, a
  vector of genotypes) from an oncosimulpop object (i.e., a simulation
  of multiple individuals) or a single oncosimul object. Sampling
  schemes include whole tumor and single cell sampling, and sampling at
  the end of the tumor progression or during the progression of the
  disease.
}

\usage{
samplePop(x, timeSample = "last", typeSample = "whole",
          thresholdWhole = 0.5)
}

\arguments{
  \item{x}{
    An object of class \code{oncosimulpop}.
  }
  
  \item{timeSample}{
    "last" means to sample each individual in the very last time period
    of the simulation. "unif" (or "uniform") means sampling each
    individual at a time choosen uniformly from all the times recorded
    in the simulation between the time when the first driver appeared
    and the final time period. "unif" means that it is almost sure that
    different individuals will be sampled at different times. "last"
    does not guarantee that different individuals will be sample at the
    same time unit, only that all will be sampled in the last time unit
    of their simulation.
  }
  
  \item{typeSample}{
    "singleCell" (or "single") for single cell sampling, where the
    probability of sampling a cell (a clone) is directly proportional to
    its population size.  "wholeTumor" (or "whole") for whole tumor
    sampling (i.e., this is similar to a biopsy being the entire tumor).
  }
  
  \item{thresholdWhole}{
    In whole tumor sampling, whether a gene is detected as mutated depends
    on thresholdWhole: a gene is considered mutated if it is altered in at
    least thresholdWhole proportion of the cells in that individual.
  }
}

\details{
  samplePop simply repeats the sampling process in each individual of
  the oncosimulpop object.
}

\value{
  A matrix. Each row is a genotype, where 0 denotes no alteration and 1
  alteration. To make it more clear that genes/laterations are in
  columns, columns are named starting by "G."  (for "gene").
}

\references{
  Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
             accumulation of mutations during tumor progression: effects
             of passengers, evolutionary models, and sampling, \emph{BMC
             Bioinformatics\/},
             \bold{16}(41). \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}.
}

\author{
  Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{oncoSimulPop}}
}

\examples{
data(examplePosets)
p705 <- examplePosets[["p705"]]

## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).

p1 <- oncoSimulPop(4, p705, mc.cores = 2)
samplePop(p1)


## Now single cell sampling

r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")
}

\keyword{manip}









