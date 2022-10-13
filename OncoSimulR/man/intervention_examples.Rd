\name{intervention-examples}
\alias{benchmark_1}
\alias{benchmark_1_0.05}
\alias{benchmark_2}
\alias{benchmark_3}

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
data(benchmark_1)
data(benchmark_1_0.05)
data(benchmark_2)
data(benchmark_3)
}

\format{
  Output from runs of oncoSimulIndiv, with some components removed to
  minimize size.
}


\examples{
data(AdaptiveTherapyComplexExample)
plot(atex2b)
}

\keyword{datasets}


