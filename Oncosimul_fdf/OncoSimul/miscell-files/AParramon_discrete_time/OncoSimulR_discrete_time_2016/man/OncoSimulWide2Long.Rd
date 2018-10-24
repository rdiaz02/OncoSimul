\name{OncoSimulWide2Long}
\alias{OncoSimulWide2Long}

\title{
  Convert the \code{pops.by.time} component of an \code{oncosimul}
  object into "long" format.
  
}
\description{
  
  Conver the \code{pops.by.time} component from its "wide" format (with
  one column for time, and as many columns as clones/genotypes) into
  "long" format, so that it can be used with other functions, for
  instance for plots.
}

\usage{
OncoSimulWide2Long(x)
}

\arguments{
  \item{x}{
    An object of class \code{oncosimul} or \code{oncosimul2}.
  }
  
}


\value{

  A data frame with four columns: Time; Y, the number of cells (the
  population size); Drivers, a factor with the number of drivers of the
  given genotype; Genotype, the genotyp.
}

\author{
  Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{oncoSimulIndiv}}
  
}

\examples{
data(examplePosets)
## An object of class oncosimul
p705 <- examplePosets[["p705"]]
p1 <- oncoSimulIndiv(p705)
class(p1)
lp1 <- OncoSimulWide2Long(p1)
head(lp1)
summary(lp1)

## An object of class oncosimul2
data(examplesFitnessEffects)

sm <-  oncoSimulIndiv(examplesFitnessEffects$cbn1,
                       model = "McFL", 
                       mu = 5e-7,
                       detectionSize = 1e8, 
                       detectionDrivers = 2,
                       sampleEvery = 0.025,
                       keepEvery = 5,
                       initSize = 2000,
                       onlyCancer = FALSE)
class(sm)
lsm <- OncoSimulWide2Long(sm)
head(lsm)
summary(lsm)
}

\keyword{manip}









