\name{freq-dep-simul-examples}
\alias{woAntibS}
%% \alias{benchmark_1_0.05}
%% \alias{benchmark_2}
%% \alias{benchmark_3}

\docType{data}
\title{
Runs from simulations of frequency-dependent examples shown in the vignette.
}

\description{

  Simulations shown in the vignette. Since running them can take a few
  seconds, we have pre-run them, and stored the results.
}
\usage{
data(woAntibS)
}
\format{
  For output from runs of \code{oncoSimulIndiv} a list of classes
  \code{oncosimul} and \code{oncosimul2}.
}


\seealso{
  \code{\link{oncoSimulIndiv}}
}

\examples{
data(woAntibS)
plot(woAntibS, show = "genotypes", type = "line",
     col = c("black", "green", "red"))
}


\keyword{datasets}


