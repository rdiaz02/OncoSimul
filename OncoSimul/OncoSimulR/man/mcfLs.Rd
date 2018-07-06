\name{mcfLs}
\alias{mcfLs}
\docType{data}
\title{
mcfLs simulation from the vignette
}

\description{

  Trimmed output from the simulation mcfLs in the vignette. This is a
  somewhat long run, and we have stored here the object (after trimming
  the Genotype matrix) to allow for plotting it.  

}
\usage{data("mcfLs")}
\format{
  An object of class "oncosimul2". A list.
}
\seealso{
  \code{\link{plot.oncosimul}}
}


\examples{
data(mcfLs)

plot(mcfLs, addtot = TRUE, lwdClone = 0.9, log = "")
summary(mcfLs)

}
\keyword{datasets}


