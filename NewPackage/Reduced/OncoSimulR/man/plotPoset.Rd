\name{plotPoset}
\alias{plotPoset}
\title{
  Plot a poset.
}
\description{
  Plot a poset. Optionally add a root and change names of nodes.
}
\usage{
plotPoset(x, names = NULL, addroot = FALSE, box = FALSE, ...)
}

\arguments{
  \item{x}{
    A poset. A matrix with two columns where, in each row, the first
    column is the ancestor and the second the descendant. Note that
    there might be multiple rows with the same ancestor, and multiple
    rows with the same descendant. See \code{\link{poset}}.
  }
  \item{names}{
    If not NULL, a vector of names for the nodes, with the same length
    as the total number of nodes in a poset (which need not be the same
    as the number of rows; see \code{\link{poset}}). If
    \code{addroot = TRUE}, then 1 + the number of nodes in the poset.
  }
  \item{addroot}{
    Add a "Root" node to the graph?
  }
  \item{box}{
    Should the graph be placed inside a box?
  }
  
  \item{\dots}{ Additional arguments to \code{plot} (actually,
    \code{plot.graphNEL} in the \code{Rgraphviz} package}).
}
\details{
  The poset is converted to a \code{graphNEL} object.
  
}
\value{
  A plot is produced.
}
\author{
  Ramon Diaz-Uriarte

}

\seealso{
  \code{\link{examplePosets}}, \code{\link{poset}}
}
\examples{
data(examplePosets)
plotPoset(examplePosets[["p1101"]])

## If you will be using that poset a lot, maybe simpler if

poset701 <- examplePosets[["p701"]]
plotPoset(poset701, addroot = TRUE)

## Compare to Pancreatic cancer figure in Gerstung et al., 2011

plotPoset(poset701,
          names = c("KRAS", "SMAD4", "CDNK2A", "TP53",
                     "MLL3","PXDN", "TGFBR2"))

## If you want to show Root explicitly do

plotPoset(poset701, addroot = TRUE,
          names = c("Root", "KRAS", "SMAD4", "CDNK2A", "TP53",
                     "MLL3","PXDN", "TGFBR2"))


## Of course, names are in the order of nodes, so KRAS is for node 1,
##   etc, but the order of entries in the poset does not matter:

poset701b <- poset701[nrow(poset701):1, ]

plotPoset(poset701b,
          names = c("KRAS", "SMAD4", "CDNK2A", "TP53",
                     "MLL3","PXDN", "TGFBR2"))
}

\keyword{ hplot }

