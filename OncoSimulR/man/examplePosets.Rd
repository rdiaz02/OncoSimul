\name{examplePosets}
\alias{examplePosets}
\docType{data}
\title{
Example posets
}

\description{

  Some example posets. For simplicity, all the posets are in a single
  list. You can access each poset by accessing each element of the
  list. The first digit or pair of digits denotes the number of nodes.
  
  
  Poset 1101 is the same as the one in Gerstung et al., 2009 (figure
  2A, poset 2). Poset 701 is the same as the one in Gerstung et al.,
  2011 (figure 2B, left, the pancreatic cancer poset). Those posets
  were entered manually at the command line: see \code{\link{poset}}.
  
  

}
\usage{data("examplePosets")}
\format{
The format is:
List of 13
$ p1101: num [1:10, 1:2] 1 1 3 3 3 7 7 8 9 10 ...
$ p1102: num [1:9, 1:2] 1 1 3 3 3 7 7 9 10 2 ...
$ p1103: num [1:9, 1:2] 1 1 3 3 3 7 7 8 10 2 ...
$ p1104: num [1:9, 1:2] 1 1 3 3 7 7 9 2 10 2 ...
$ p901 : num [1:8, 1:2] 1 2 4 5 7 8 5 1 2 3 ...
$ p902 : num [1:6, 1:2] 1 2 4 5 7 5 2 3 5 6 ...
$ p903 : num [1:6, 1:2] 1 2 5 7 8 1 2 3 6 8 ...
$ p904 : num [1:6, 1:2] 1 4 5 5 1 7 2 5 8 6 ...
$ p701 : num [1:9, 1:2] 1 1 1 1 2 3 4 4 5 2 ...
$ p702 : num [1:6, 1:2] 1 1 1 1 2 4 2 3 4 5 ...
$ p703 : num [1:6, 1:2] 1 1 1 1 3 5 2 3 4 5 ...
$ p704 : num [1:6, 1:2] 1 1 1 1 4 5 2 3 4 5 ...
$ p705 : num [1:6, 1:2] 1 2 1 1 1 2 2 5 4 6 ...
}
\source{
  Gerstung et al., 2009. Quantifying cancer progression with conjunctive
  Bayesian networks. \emph{Bioinformatics}, 21: 2809--2815.
  
  Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
  Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
  
}

\seealso{
  \code{\link{poset}}
}


\examples{
data(examplePosets)

## Plot all of them
par(mfrow = c(3, 5))

invisible(sapply(names(examplePosets),
                 function(x) {plotPoset(examplePosets[[x]],
                     main = x,
                     box = TRUE)}))

}
\keyword{datasets}


