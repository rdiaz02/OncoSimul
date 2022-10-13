\name{intervention-examples}
\alias{osi}
\alias{osi_with_ints}
\alias{atex4}
\alias{atex5}
\alias{atex2b}
\alias{uvex3}


\docType{data}
\title{
Runs from simulations of interventions examples shown in the
vignette. 
}

\description{
Simulations shown in the vignette. Since running them can take a few
seconds, we have pre-run them, and stored the results.

 They are here mainly to facilitate creation of table from the vignette
 itself. The script is available under "inst/miscell". 

}

\usage{
data(osi)
data(osi_with_ints)
data(atex4)
data(atex5)
data(atex2b)
data(uvex3)
}

\format{
  Output from runs of oncoSimulIndiv, with some components removed to
  minimize size.
}


\examples{
data(atex2b)
plot(atex2b)
}

\keyword{datasets}


