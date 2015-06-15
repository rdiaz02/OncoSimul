\name{poset}
\alias{poset}

\title{
  Poset
}
\description{
  Poset: explanation.
}

\arguments{
  \item{x}{
    The poset. See details.
  }
}
\details{
  A poset is a two column matrix. In each row, the first column is the
  ancestor (or the restriction) and the second column the descendant (or
  the node that depends on the restriction). Each node is identified by
  a positive integer. The graph includes all nodes with integers between
  1 and the largest integer in the poset.
  
  Each node can be necessary for several nodes: in this case, the same
  node would appear in the first column in several rows.
  
  A node can depend on two or more nodes (conjunctions): in this case,
  the same node would appear in the second column in several rows.
  
  There can be nodes that do not depend on anything (except the Root
  node) and on which no other nodes depend. The simplest and safest way
  to deal with all possible cases, including these cases, is to have all
  nodes with at least one entry in the poset, and nodes that depend on
  no one, and on which no one depends should be placed on the second
  column (with a 0 on the first column).
  
  Alternatively, any node not named explicitly in the poset, but with a
  number smaller than the largest number in the poset, is taken to be a
  node that depends on no one and on which no one depends. See examples
  below.

  This specification of restrictions is for version 1. See
  \code{\link{allFitnessEffects}} for a much more flexible one for
  version 2. Both can be used with  \code{\link{oncoSimulIndiv}}.
}

\references{
  Posets and similar structures appear in several places. The following
  two papers use them extensively.
  
  Gerstung et al., 2009. Quantifying cancer progression with conjunctive
  Bayesian networks. \emph{Bioinformatics}, 21: 2809--2815.
  
  Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
  Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
}
\author{
  Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{examplePosets}}, \code{\link{plotPoset}},
  \code{\link{oncoSimulIndiv}}
}

\examples{

## Node 2 and 3 depend on 1, and 4 depends on no one
p1 <- cbind(c(1L, 1L, 0L), c(2L, 3L, 4L))
plotPoset(p1, addroot = TRUE)

## Node 2 and 3 depend on 1, and 4 to 7 depend on no one.
## We do not have nodes 4 to 6 explicitly in the poset.
p2 <- cbind(c(1L, 1L, 0L), c(2L, 3L, 7L))
plotPoset(p2, addroot = TRUE)

## But this is arguably cleaner
p3 <- cbind(c(1L, 1L, rep(0L, 4)), c(2L, 3L, 4:7 ))
plotPoset(p3, addroot = TRUE)

## A simple way to create a poset where no gene (in a set of 15) depends
## on any other.

p4 <- cbind(0L, 15L)
plotPoset(p4, addroot = TRUE)


## Specifying the pancreatic cancer poset in Gerstung et al., 2011
##   (their figure 2B, left). We use numbers, but for nicer plotting we
##   will use names: KRAS is 1, SMAD4 is 2, etc.

pancreaticCancerPoset <- cbind(c(1, 1, 1, 1, 2, 3, 4, 4, 5),
                               c(2, 3, 4, 5, 6, 6, 6, 7, 7))
storage.mode(pancreaticCancerPoset) <- "integer"

plotPoset(pancreaticCancerPoset,
          names = c("KRAS", "SMAD4", "CDNK2A", "TP53",
                    "MLL3","PXDN", "TGFBR2"))


## Specifying poset 2 in Figure 2A of Gerstung et al., 2009:

poset2 <- cbind(c(1, 1, 3, 3, 3, 7, 7, 8, 9, 10),
                c(2, 3, 4, 5, 6, 8, 9, 10, 10, 11))

storage.mode(poset2) <- "integer"
plotPoset(poset2)

}

\keyword{manip}

