\name{poset2AdjMat}
\alias{poset2AdjMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtain an adjacency matrix from a poset.
}
\description{
  Posets are probably simpler to prepare by hand than adjacency matrices, but
  our simulation functions use, as input, adjacency matrices. 
}
\usage{
poset2AdjMat(x)
}
%- maybe also 'usage' for other objects documented here.
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

  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An adjacency matrix (for a directed graph).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Posets and similar structures appear in several places. The following
  two papers use them extensively.

   Gerstung et al., 2009. Quantifying cancer progression with conjunctive
  Bayesian networks. \emph{Bioinformatics}, 21: 2809--2815.

  Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
  Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
%% ~put references to the literature/web site here ~
}
\author{
  Ramon Diaz-Uriarte
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{examplePosets}}, \code{\link{plotPoset}}, \code{\link{oncoSimulIndiv}}
  %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Node 2 and 3 depend on 1, and 4 depends on no one

p1 <- cbind(c(1, 1, 0), c(2, 3, 4))
poset2AdjMat(p1)
plotPoset(p1, addroot = TRUE)

## Node 2 and 3 depend on 1, and 4 to 7 depend on no one.
## We do not have nodes 4 to 6 explicitly in the poset.
p2 <- cbind(c(1, 1, 0), c(2, 3, 7))
poset2AdjMat(p2)
plotPoset(p2, addroot = TRUE)

## But this is arguably cleaner
p3 <- cbind(c(1, 1, rep(0, 4)), c(2, 3, 4:7 ))
poset2AdjMat(p3)
plotPoset(p3, addroot = TRUE)

identical(poset2AdjMat(p3), poset2AdjMat(p2))

## A simple way to create a poset where no gene (in a set of 15) depends
## on any other.

p4 <- cbind(0, 15)
poset2AdjMat(p4)
plotPoset(p4, addroot = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

