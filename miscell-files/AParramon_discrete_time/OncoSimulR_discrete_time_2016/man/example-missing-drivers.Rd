\name{example-missing-drivers}
\alias{ex_missing_drivers_b11}
\alias{ex_missing_drivers_b12}
\docType{data}
\title{
An example where there are intermediate missing drivers.
}

\description{

  An example where there are intermediate missing drivers. This is
fictitious and I've never seen it. But it is here to check plots work
even if there are no cases of some intermediate value of drivers (2 in
this case). b11 contains the full, original data, whereas b12 contains
the same data where there are no cases with exactly 2 drivers.


}
\usage{data("ex_missing_drivers_b11"); data("ex_missing_drivers_b12")}
\format{
  Two objects of class "oncosimul". 
}
\seealso{
  \code{\link{plot.oncosimul}}
}


\examples{
data(ex_missing_drivers_b11)
plot(ex_missing_drivers_b11, type = "line")
dev.new()
data(ex_missing_drivers_b12)
plot(ex_missing_drivers_b12, type = "line")
}
\keyword{datasets}


