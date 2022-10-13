\name{intervention-examples}
\alias{osi_intex2}
\alias{osi_with_ints}
\alias{HansenExample3}
\alias{HansenExampleAT4}
\alias{AdaptiveTherapyComplexExample}
\alias{userVarsBasicExample2}

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
data(osi_intex2)
data(osi_with_ints)
data(HansenExample3)
data(HansenExampleAT4)
data(AdaptiveTherapyComplexExample)
data(userVarsBasicExample2)
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


