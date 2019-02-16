\name{plotFitnessLandscape}
\alias{plotFitnessLandscape}
\alias{plot.genotype_fitness_matrix}
\alias{plot.evalAllGenotypes}
\alias{plot.evalAllGenotypesMut}



\title{ Plot a fitness landscape.  }

\description{ Show a plot of a fitness landscape. The plot is modeled
  after (actually, mostly a blatant copy of) that of MAGELLAN,
  \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}.

  Note: this is not a plot of the fitnessEffects object; for that, see
  \code{\link{plot.fitnessEffects}}.

}

\usage{
plotFitnessLandscape(x, show_labels = TRUE,
                     col = c("green4", "red", "yellow"),
                     lty = c(1, 2, 3),
                     use_ggrepel = FALSE,
                     log = FALSE, max_num_genotypes = 2000,
                     only_accessible = FALSE,
                     accessible_th = 0,
                     ...)

\method{plot}{genotype_fitness_matrix}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max_num_genotypes = 2000,
                                   only_accessible = FALSE,
                                   accessible_th = 0,
                                   ...)

\method{plot}{evalAllGenotypes}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max_num_genotypes = 2000,
                                   only_accessible = FALSE,
                                   accessible_th = 0,
                                   ...)

\method{plot}{evalAllGenotypesMut}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max_num_genotypes = 2000,
                                   only_accessible = FALSE,
                                   accessible_th = 0,
                                   ...)
}

\arguments{
  
  \item{x}{ One of the following:
    \itemize{

      \item A matrix (or data frame) with g + 1 columns. Each of the
      first g columns contains a 1 or a 0 indicating that the gene of
      that column is mutated or not. Column g+ 1 contains the fitness
      values. This is, for instance, the output you will get from
      \code{\link{rfitness}}.

      \item A two column data frame. The second column is fitness, and
      the first column are genotypes, given as a character vector. For
      instance, a row "A, B" would mean the genotype with both A and B
      mutated.
      
      \item The output from a call to
      \code{\link{evalAllGenotypes}}. Make sure you use \code{order =
	FALSE} in that call.
      
      \item The output from a call to
      \code{\link{evalAllGenotypesMut}}. Make sure you use \code{order =
	FALSE}.
      
      \item The output from a call to \code{\link{allFitnessEffects}}.
      
    }
    
    The first two are the same as the format for the \code{genotFitness}
    component in \code{\link{allFitnessEffects}}.  }
  
  \item{show_labels}{If TRUE, show the genotype labels.}

  \item{col}{A three-element vector that  gives the colors to use for
    increase, decreases and no changes in fitness, respectively. The first
    two colours are also used for peaks and sinks.}
  
  \item{lty}{A three-element vector that  gives the line types to use for
    increase, decreases and no changes in fitness, respectively.}

  \item{use_ggrepel}{If TRUE, use the \link{ggrepel} package to avoid
    overlap of labels.}

  \item{log}{Log-scale the y axis.}

  \item{max_num_genotypes}{Maximum allowed number of genotypes. For some
  types of input, we make a call to \code{\link{evalAllGenotypes}}, and
  use this as the maximum.}

\item{only_accessible}{If TRUE, show only accessible paths. A path is
  considered accesible if, at each mutational step (i.e., with the
  addition of each mutation) fitness increases by at least
  \code{accessible_th}. If you set \code{only_accessible = TRUE}, the
  number of genotypes displayed can be much smaller than the number of
  existing genotypes if many of those genotypes are not accessible via
  any path.}

\item{accessible_th}{The threshold for the minimal change in fitness at
  each mutation step (i.e., between successive genotypes) to be used if
  \code{only_accessible = TRUE}.}

  \item{\dots}{
      Other arguments passed to \code{plot}. Not used for now.
  }
}

\value{

  A fitness landscape plot: a plot showing paths between
  genotypes and peaks and sinks (local maxima and minima).

}


\author{
Ramon Diaz-Uriarte
}

\references{
  MAGELLAN web site: \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}

  Brouillet, S. et al. (2015). MAGELLAN: a tool to explore small fitness landscapes. \emph{bioRxiv}, \bold{31583}. \url{http://doi.org/10.1101/031583}

}


\note{

  I have copied most of the ideas (and colors, and labels) of this plot
  from MAGELLAN (\url{http://wwwabi.snv.jussieu.fr/public/Magellan/})
  but MAGELLAN has other functionality that is not provided here such as
  epistasis stats for the landscape, and several visual manipulation
  options.

  One feature of this function that is not available in MAGELLAN is
  showing genotype labels (i.e., annotated by gene names), which can be
  helpful if the different genotypes mean something to you.

  
  In addition to the above differences, another difference between this
  plot and those of MAGELLAN is \bold{how sinks/peaks of more than one
  genotype are dealt with}. This plot will show as sinks or peaks sets
  of one or more genotypes that are of identical fitness (and separated
  by a Haming distance of one). So a sink or a peak might actually be
  made of more than one genotype. In MAGELLAN, as far as I can tell,
  peaks and sinks are always made of a single isolated genotype.

  Does this matter? In most realistic cases where not two genotypes can
  have exactly the same fittnes it does not. In some cases, though, it
  might matter. Are multi-genotype sinks/peaks really sinks/peaks?
  Arguably yes: suppose genotypes "AB" and "ABC" both have fitness 0,
  which is minimal among the fitness in the set of genotypes, and
  genotypes "A" and "ABCD" have fitness 0.1. To go from "A" to "ABCD",
  if you want to travel through "AB", you have to go through the valley
  of "AB" and "ABC"; once in "ABC" you can climb up to "ABCD"; and once
  in "AB" you can move to "ABC" since it has identical fitness to
  "AB". Mutatis mutandis for multi-genotype peaks. Ignoring the
  possibility of peaks/sinks made of more than one genotype actually
  makes code much simpler.

  
  Sometimes not showing the any links that involve a decrease in fitness
  can help see non-accessible pathways (in strong selection, no multiple
  mutations, etc); do this by passing, for instance, an NA for the
  second element of \code{col}.

  Finally, use common sense: for instance, if you pass a
  \code{\link{allFitnessEffects}} that specifies for, say, the fitness
  of a total of 5000 genotypes you'll have to wait a while for the plot
  to finish.

  

}

\seealso{
  \code{\link{allFitnessEffects}},
  \code{\link{evalAllGenotypes}},
  \code{\link{allFitnessEffects}},
  \code{\link{rfitness}},
  \code{\link{plot.fitnessEffects}}
}
\examples{

## Generate random fitness for four genes-genotypes
## and plot landscape.

r1 <- rfitness(4)
plot(r1)


## Specify fitness in a matrix, and plot it

m5 <- cbind(A = c(0, 1, 0, 1), B = c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
plotFitnessLandscape(m5)


## Specify fitness with allFitnessEffects, and plot it

fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))

plot(evalAllGenotypes(fe, order = FALSE))

## same as
plotFitnessLandscape(evalAllGenotypes(fe, order = FALSE))

}

\keyword{ hplot }

