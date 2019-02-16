\name{plotClonePhylog}
\alias{plotClonePhylog}
\title{
  Plot a parent-child relationship of the clones.
}

\description{

  Plot a parent-child relationship of the clones, controlling which
  clones are displayed, and whether to shown number of times of
  appearance, and time of first appearance of a clone. 


}


\usage{
plotClonePhylog(x, N = 1, t ="last", timeEvents = FALSE,
                 keepEvents = FALSE, fixOverlap = TRUE,
                 returnGraph = FALSE, ...)
}

\arguments{

  \item{x}{ The output from a simulation, as obtained from
  \code{oncoSimulIndiv}, \code{oncoSimulPop}, or \code{oncoSimulSample}
  (see \code{\link{oncoSimulIndiv}}). This must be from v.2 and forward
  (no phylogenetic information is stored for earlier objects).
}

\item{N}{ Show in the plot all clones that have a population size of at
  least N at time \code{time} and the parents of those clones (parents
  are shown regardless of population size ---i.e., you can see extinct
  parents).  If you want to show everything that ever appeared, set N = 0.}

\item{t}{ The time at which \code{N} should be satisfied. This can
  either be the string "last", meaning the last time of the simulation,
  or a range of two values. In the second case, all clones with
  population size of at least \code{N} in at least one time point
  between time[1] and time[2] will be shown (togheter with their parents).}

\item{timeEvents}{If TRUE, the vertical position of the nodes in the
  plot will be proportional to their time of first appearance. }

\item{keepEvents}{If TRUE, the graph will show all the birth
  events. Thus, the number of arrows shows the number of times a clone
  give rise to another. For large graphs with many events, this slows
  the graph considerably.}

\item{fixOverlap}{When using \code{timeEvents = TRUE} nodes can overlap
  (as we modify their vertical location after igraph has done the
  initial layout). This attempts to fix that problem by randomly
  relocating, along the X axis, the nodes that have the same X value.}

\item{returnGraph}{If TRUE, the igraph object is returned. You can use
  this to plot the object however you want or obtain the adjacency matrix.}
  
\item{\dots}{ Additional arguments. Currently not used.}.
  
}

\note{

  These are not, technically, proper phylogenetic trees and we use
  "phylogeny" here in an abuse of terminology. The plots we use, where
  we show parent child relationships are arguably more helpful in this
  context. But you could draw proper phylogenies with the information
  provided.

  
  
  If you want to obtain the adjacency matrix, this is trivial: just
  set \code{returnGraph = TRUE} and use \code{\link{get.adjacency}}. See
  an example below.  
}
\value{
  A plot is produced. If \code{returnGraph} the igraph object is returned.
}
\author{
  Ramon Diaz-Uriarte

}

\seealso{
  \code{\link{oncoSimulIndiv}}
}
\examples{

data(examplesFitnessEffects)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE)

## Show only those with N > 10 at end
plotClonePhylog(tmp, N = 10)

## Show only those with N > 1 between times 5 and 1000
plotClonePhylog(tmp, N = 1, t = c(5, 1000))

## Show everything, even if teminal nodes are extinct
plotClonePhylog(tmp, N = 0)


## Show time when first appeared
plotClonePhylog(tmp, N = 10, timeEvents = TRUE)

\dontrun{
## Show each event
## This can take a few seconds
plotClonePhylog(tmp, N = 10, keepEvents = TRUE)
}

## Adjacency matrix
require(igraph)
get.adjacency(plotClonePhylog(tmp, N = 10, returnGraph = TRUE))

}

\keyword{ hplot }

