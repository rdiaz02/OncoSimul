\name{plot.fitnessEffects}
\alias{plot.fitnessEffects}
%%\alias{plotFitnessEffects}
\title{
  Plot fitnessEffects objects.
}
\description{
  Plot the restriction table/graph of restrictions, the epistasis, and
  the order effects in a fitnessEffects object. This is not a plot of
  the fitness landscape; for that, see \code{\link{plotFitnessLandscape}}.
}
\usage{
\method{plot}{fitnessEffects}(x, type = "graphNEL", layout = NULL,
expandModules = FALSE, autofit = FALSE,
scale_char = ifelse(type == "graphNEL", 1/10, 5),
return_g = FALSE, lwdf = 1, ...)

%%plotFitnessEffects(fe, type = "graphNEL", layout = NULL, expandModules = FALSE)
}

\arguments{
  \item{x}{A \code{fitnessEffects} object, as produced by \code{\link{allFitnessEffects}}.
}
\item{type}{
  Whether you want a "graphNEL" or an "igraph" graph.
}
\item{layout}{
  For "igraph", the layout. For example, if you know you really have
  only a tree you might want to use \code{layout.reingold.tilford}. Note
  that there is very limited support for passing options, etc. In most
  cases, it is either the default or the \code{layout.reingold.tilford}.
}
  \item{expandModules}{
   If there are modules with multiple genes, if you set this to TRUE
  modules will be replaced by their genes.
}

\item{autofit}{If TRUE, we try to fit the edges to the labels. This is a
  very experimental feature, likely to be not very robust.}

\item{scale_char}{If using \code{autofit = TRUE}, the scaling factor for
the size of the rectangles as a function of the number of
characters. You have to play with this because the best value can
depend on a number of things.}

\item{return_g}{It TRUE, the graph object (graphNEL or igrap) is
  returned.}

\item{lwdf}{The multiplier factor for \code{lwd} when using "graphNEL".}


 \item{\dots}{
    Other arguments passed to \code{plot}. Not used for now.
  }

}
\value{A plot.

  Order and epistatic relationships have orange edges. OR (semimonotone)
  relationships blue, and XOR red. All others have black edges (so AND
  and unique edges from root). Epistatic relationships, being
  symmetrical, have no arrows between nodes and have a dotted line
  type. Order relationships have an arrow from the earlier to the later
  event and have a different dotted line (lty 3).

  If \code{return_g} is TRUE, you are returned also the graph object
  (igraph or graphNEL) so that you can manipulate it further.
}
\author{
Ramon Diaz-Uriarte
}
\note{

  The purpose of the plot is to get a quick idea of the
relationships. Note that three-way (or higher order) epistatic
relationships cannot be shown as such (we would show all possible pairs,
but that is not quite the same thing). Likewise, there is no reasonable
way to convey the pressence of a "-" in the epistatic relationship.

Genes without interactions are not shown.

}

\seealso{
  \code{\link{allFitnessEffects}},
  \code{\link{plotFitnessLandscape}}
  
}
\examples{

cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")

cbn1 <- allFitnessEffects(cs)
plot(cbn1, "igraph")

library(igraph) ## to make layouts available
plot(cbn1, "igraph", layout = layout.reingold.tilford)


### A DAG with the three types of relationships
p3 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c", "f"),
                  child = c("a", "b", "d", "e", "c", "c", "f", "f", "g", "g"),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                  sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                  typeDep = c(rep("--", 4), 
                      "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
fp3 <- allFitnessEffects(p3)

plot(fp3)

plot(fp3, "igraph", layout = layout.reingold.tilford)


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

noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)

fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                         noIntGenes = noint, geneToModule = modules)


plot(fea)
plot(fea, expandModules = TRUE)
plot(fea, type = "igraph")

}

\keyword{ hplot }

