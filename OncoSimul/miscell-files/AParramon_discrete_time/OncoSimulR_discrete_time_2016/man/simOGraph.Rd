\name{simOGraph}
\alias{simOGraph}
\title{
  Simulate oncogenetic/CBN/XMPN DAGs.
}
\description{
  Simulate DAGs that represent restrictions in the accumulation of mutations.
}
\usage{
simOGraph(n, h = 4, conjunction = TRUE, nparents = 3,
multilevelParent = TRUE, removeDirectIndirect = TRUE, rootName = "Root")
}

\arguments{
  \item{n}{
    Number of nodes, or edges, in the graph. Like the number of genes.
}
  \item{h}{
    Approximate height of the graph. See details.
  }
  \item{conjunction}{
    If TRUE, conjunctions (i.e., multiple parents for a node) are allowed.
}
\item{nparents}{
  Maximum number of parents of a node, when \code{conjunction} is TRUE.
}
\item{multilevelParent}{
  Can a node have parents at different heights (i.e., parents that are
  at different distance from the root node)?
}
  \item{removeDirectIndirect}{
    Ensure that no two nodes are connected both directly (i.e., with an
  edge between them) and indirectly, through intermediate nodes. If
  TRUE, the direct connections are removed from the graph starting from
  the bottom.
}
\item{rootName}{
  The name you want to give the "Root" node.

}
}
\details{
  This is a simple, heuristic procedure for generating graphs of
  restrictions that seem compatible with published trees in the
  oncogenetic literature.

  The basic procedure is as follows: nodes (argument \code{n}) are split
  into approximately equally sized \code{h} groups, and then each node
  from a level is connected to nodes chosen randomly from nodes of the
  remaing superior (i.e., closer to the Root) levels. The number of
  edges comes from a uniform distribution between 1 and \code{nparents}.

  The actual depth of the graph can be smaller than \code{h} because
  nodes from a level might be connected to superior levels skipping
  intermediate ones.

  See the vignette for further discussion about arguments. 

}
\value{
  An adjacency matrix for a directed graph.
}
\author{Ramon Diaz-Uriarte}

\examples{
(a1 <- simOGraph(10))
library(graph) ## for simple plotting
plot(as(a1, "graphNEL"))
}

\keyword{ graphs }
\keyword{ datagen }
