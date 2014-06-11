\name{samplePop}
\alias{samplePop}

\title{
  Obtain a sample from a population of simulations.
  %%  ~~function to do ... ~~
}
\description{
  
  Obtain a sample (a matrix of individuals by genes or, equivalently, a
  vector of genotypes) from an oncosimulpop object (i.e., a simulation
  of multiple individuals) or a single oncosimul object. Sampling
  schemes include whole tumor and single cell sampling, and sampling at
  the end of the tumor progression or during the progression of the
  disease.
  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
samplePop(x, timeSample = "last", typeSample = "whole", thresholdWhole = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    An object of class \code{oncosimulpop}.
%%     ~~Describe \code{x} here~~
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

  
%%     ~~Describe \code{typeSample} here~~
}
\item{thresholdWhole}{
  In whole tumor sampling, whether a gene is detected as mutated depends
  on thresholdWhole: a gene is considered mutated if it is altered in at
  least thresholdWhole proportion of the cells in that individual.
%%     ~~Describe \code{thresholdWhole} here~~
}
}
\details{

  samplePop simply repeats the sampling process in each individual of
  the oncosimulpop object.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A matrix. Each row is a genotype, where 0 denotes no alteration and 1
  alteration.

}
\references{
    Diaz-Uriarte, R. (2014). Inferring restrictions in the temporal order
of mutations during tumor progression: effects of passenger mutations,
evolutionary models, and sampling. \url{http://dx.doi.org/10.1101/005587}

%% ~put references to the literature/web site here ~
}
\author{
  Ramon Diaz-Uriarte
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{oncoSimulPop}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(examplePosets)
p705 <- examplePosets[["p705"]]
p1 <- oncoSimulPop(20, p705, detectionSize = 1e6,
                   sampleEvery = 10, mc.cores = 4)
samplePop(p1)

r1 <- oncoSimulIndiv(p705, detectionSize = 1e6)
samplePop(r1, typeSample = "single")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

