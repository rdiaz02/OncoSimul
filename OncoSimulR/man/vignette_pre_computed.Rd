\name{vignette_pre_computed}
\alias{osi}
\alias{osi_with_ints}
\alias{atex4}
\alias{atex5}
\alias{atex2b}
\alias{uvex3}
\alias{smyelo3v57}
\alias{s_3_b}
\alias{uvex2}
\alias{simT2}
\alias{simul_period_1}
\alias{simT3}
\alias{s_3_a}


\docType{data}
\title{
Runs from simulations of interventions examples shown in the
vignette. Most, but not all, are from intervention examples.
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
data(smyelo3v57)
data(s_3_b)
data(uvex2)
data(simT2)
data(simul_period_1)
data(simT3)
data(s_3_a)
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


