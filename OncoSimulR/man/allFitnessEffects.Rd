\name{allFitnessEffects}
\alias{allFitnessEffects}
\title{Create fitness effects specification from restrictions,
  epistasis, and order effects.
}
\description{
  Given one or more of a set of poset restrictions, epistatic
  interactions, order effects, and genes without interactions, as well
  as, optionally, a mapping of genes to modules, return the complete
  fitness specification.

  The output of this function is not intended for user consumption, but
  as a way of preparing data to be sent to the C++ code.  }


\usage{

allFitnessEffects(rT = NULL, epistasis = NULL, orderEffects = NULL,
  noIntGenes = NULL, geneToModule = NULL, drvNames = NULL, keepInput =
  TRUE) }

\arguments{
  \item{rT}{A restriction table that is an extended version of a poset 
    (see \code{\link{poset}} ).
    A restriction table is a data frame where each row shows one edge
    between a parent and a child. A restriction table contains exactly these
    columns, in this order:
    \describe{
    \item{parent}{The identifiers of the parent nodes, in a
      parent-child relationship. There must be at least on entry with the
      name "Root".}
    \item{child}{The identifiers of the child nodes.}
    \item{s}{A numeric vector with the fitness effect that applies
      if the relationship is satisfied.}
    \item{sh}{A numeric vector with the fitness effect that applies if
      the relationship is not satisfied. This provides a way of
      explicitly modeling deviatons from the restrictions in the graph,
      and is discussed in Diaz-Uriarte, 2015. }
    \item{typeDep}{The type of dependency. Three possible types of
      relationship exist:
      \describe{
	\item{AND, monotonic, or CMPN}{Like in the CBN model, all parent nodes
	  must be present for a relationship to be satisfied. Specify it
	  as "AND" or "MN" or "monotone".}
	\item{OR, semimonotonic, or DMPN}{A single parent node is enough
	  for a relationship to be satisfied. Specify it as "OR" or
	  "SM" or "semimonotone".}
	\item{XOR or XMPN}{Exactly one parent node must be mutated for a
	  relationship to be satisfied. Specify it as "XOR" or "xmpn" or
	  "XMPN".}
      }

      In addition, for the nodes that depend only on the root node, you
      can use "--" or "-" if you want (though using any of the other
      three would have the same effects if a node that connects to root
      only connects to root).}
  }
  }
  \item{epistasis}{
    A named numeric vector. The names identify the relationship, and the
    numeric value is the fitness effect. For the names, each of the
    genes or modules involved is separated by a ":". A negative sign
    denotes the absence of that term.
  }
  \item{orderEffects}{
    A named numeric vector, as for \code{epistasis}. A ">" separates the
  names of the genes of modules of a relationship, so that "U > Z" means
  that the relationship is satisfied when mutation U has happened before
  mutation Z.
}
\item{noIntGenes}{
  A numeric vector (optionally named) with the fitness coefficients of genes
  (only genes, not modules) that show no interactions. These genes
  cannot be part of modules. But you can specify modules that have
  no epistatic interactions. See examples and vignette.

  Of course, avoid using potentially confusing characters in the
  names. In particular, "," and ">" are not allowed as gene names.
}

\item{geneToModule}{
      
  A named character vector that allows to match genes and modules. The
  names are the modules, and each of the values is a character vector
  with the gene names, separated by a comma, that correspond to a
  module. Note that modules cannot share genes. There is no need for
  modules to contain more than one gene.  If you specify a geneToModule
  argument, and you used a restriction table, the \code{geneToModule} 
  must necessarily contain, in the first position, "Root" (since the
  restriction table contains a node named "Root"). See examples below.
}
    

\item{drvNames}{The names of genes that are considered drivers. This is
  only used for: a) deciding when to stop the simulations, in case you
  use number of drivers as a simulation stopping criterion (see
  \code{\link{oncoSimulIndiv}}); b) for summarization purposes (e.g.,
  how many drivers are mutated); c) in figures. But you need not
  specifiy anything if you do not want to, and you can pass an empty
  vector (as \code{character(0)}). The default is to assume that all
  genes that are not in the \code{noIntGenes} are drivers.}

\item{keepInput}{
  If TRUE, whether to keep the original input. This is only useful for
  human consumption of the output. It is useful because it is easier to
  decode, say, the restriction table from the data frame than from the
  internal representation. But if you want, you can set it to FALSE and
  the object will be a little bit smaller.}
}
\details{
  This function is used for extremely flexible specification of fitness
  effects, including posets, XOR relationships, synthetic mortality and
  synthetic viability, arbitrary forms of epistatis, arbitrary forms of
  order effects, etc. Please, see the vignette for detailed and
  commented examples.
}

\value{
  An object of class "fitnessEffects". This is just a list, but it is not
  intended for human consumption.  The components are:

  \item{long.rt}{The restriction table in "long format", so as to be
    easy to parse by the C++ code.}

  \item{long.epistasis}{Ditto, but for the epistasis specification.}

  \item{long.orderEffects}{Ditto for the order effects.}

  \item{long.geneNoInt}{Ditto for the non-interaction genes.}

  \item{geneModule}{Similar, for the gene-module correspondence.}

  \item{graph}{An \code{igraph} object that shows the restrictions,
    epistasis and order effects, and is useful for plotting.}
  
  \item{drv}{The numeric identifiers of the drivers. The numbers
    correspond to the internal numeric coding of the genes.}

  \item{rT}{If \code{keepInput} is TRUE, the original restriction
    table.}

  \item{epistasis}{If \code{keepInput} is TRUE, the original epistasis
  vector.}

  \item{orderEffects}{If \code{keepInput} is TRUE, the original order
  effects vector.}

  \item{noIntGenes}{If \code{keepInput} is TRUE, the original 
    noIntGenes.}
}
\references{
    Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
  accumulation of mutations during tumor progression: effects of
  passengers, evolutionary models, and sampling
  \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}

    McFarland, C.~D. et al. (2013). Impact of deleterious passenger
  mutations on cancer progression.  \emph{Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{110}(8),
  2910--5.


}

\note{
  Please, note that the meaning of the fitness effects in the
  McFarland model is not the same as in the original paper; the fitness
  coefficients are transformed to allow for a simpler fitness function
  as a product of terms. This differs with respect to v.1. See the
  vignette for details.

  The names of the genes and modules can be fairly arbitrary. But if you
  try hard you can confuse the parser. For instance, using gene or
  module names that contain "," or ":", or ">" is likely to get you into
  trouble. Of course, you know you should not try to use those
  characters because you know those characters have special meanings to
  separate names or indicate epistasis or order relationships.  Right
  now, using those characters as names is caught (and result in
  stopping) if passed as names for noIntGenes.  }

\author{ Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{evalGenotype}}, \code{\link{oncoSimulIndiv}}, \code{\link{plot.fitnessEffects}}
}
\examples{
## A simple poset or CBN-like example

cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")

cbn1 <- allFitnessEffects(cs)

plot(cbn1)


## A more complex example, that includes a restriction table
## order effects, epistasis, genes without interactions, and moduels
p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                 child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                 sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                 typeDep = c(rep("--", 4), 
                     "XMPN", "XMPN", "MN", "MN", "SM", "SM"))

oe <- c("C > F" = -0.1, "H > I" = 0.12)
sm <- c("I:J"  = -1)
sv <- c("-K:M" = -.5, "K:-M" = -.5)
epist <- c(sm, sv)

modules <- c("Root" = "Root", "A" = "a1",
             "B" = "b1, b2", "C" = "c1",
             "D" = "d1, d2", "E" = "e1",
             "F" = "f1, f2", "G" = "g1",
             "H" = "h1, h2", "I" = "i1",
             "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")

set.seed(1) ## for repeatability
noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)

fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                         noIntGenes = noint, geneToModule = modules)

plot(fea)


## Modules that show, between them,
## no epistasis (so multiplicative effects).
## We specify the individual terms, but no value for the ":".

fnme <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2),
                          geneToModule = c("A" = "a1, a2",
                                           "B" = "b1"))

evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)

}

\keyword{ manip }
\keyword{ list }
