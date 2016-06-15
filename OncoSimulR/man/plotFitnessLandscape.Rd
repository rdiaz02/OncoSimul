\name{plotFitnessLandscape}
\alias{plotFitnessLandscape}
\alias{plot.genotype_fitness_matrix}
\alias{plot.evalAllGenotypes}
\alias{plot.evalAllGenotypesMut}



\title{ Plot a fitness landscape.  }

\description{ Show a plot of a fitness landscape. The plot is modeled
  after (actually, a blatant copy of) that of MAGELLAN,
  \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}. 

}

\usage{
plotFitnessLandscape(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max.num.genotypes = 2000, ...)

\method{plot}{genotype_fitness_matrix}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max.num.genotypes = 2000, ...)

\method{plot}{evalAllGenotypes}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max.num.genotypes = 2000, ...)

\method{plot}{evalAllGenotypesMut}(x, show_labels = TRUE,
                                   col = c("green4", "red", "yellow"),
                                   lty = c(1, 2, 3),
                                   use_ggrepel = FALSE,
                                   log = FALSE, max.num.genotypes = 2000, ...)
}

\arguments{
  
  \item{x}{ One of the following:
    
    \describe{
      
      \item A matrix (or data frame) with g + 1 columns, where
      each of the first g columns contains a 1 or a 0 indicating that the
      gene of that column is mutated or not. Column g+ 1 contains the
      fitness values. This is, for instance, the output you will get from
      \code{\link{rfitness}}. 
      
      \item A two column data frame. The second column is fitness, and
      the first column are genotypes, given as a character vector. For
      instance, a row "A, B" would mean the genotype with both A and B
      mutated.
      
      \item The output from a call to \code{\link{evalAllGenotypes}}.
      
      \item The output from a call to \code{\link{evalAllGenotypesMut}}.

      \item The output from a call to \code{\link{allFitnessEffects}}.
      
    }
    The first two are the same as the format for the \code{genotFitness}
    component in \code{\link{allFitnessEffects}}.
  }
  
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

  Brouillet, S. et al. (2015). MAGELLAN: a tool to explore small fitness landscapes. \emph{bioRxiv\/}, \bold{31583}. \url{http://doi.org/10.1101/031583}

}


\note{

  I have copied most of the ideas (and colors, and labels) of this plot
  from MAGELLAN (\url{http://wwwabi.snv.jussieu.fr/public/Magellan/})
  but MAGELLAN has other functionality that is not provided here (e.g.,
  epistasis stats for the landscape, and several visual manipulation
  options).

  Use common sense: for instance, if you pass a
  \code{\link{allFitnessEffects}} that specifies for, say, the fitness
  of a total of 5000 genotypes you'll have to wait a while for the plot
  to finish.

 

}

\seealso{
  \code{\link{allFitnessEffects}},
  \code{\link{evalAllGenotypes}},
  \code{\link{allFitnessEffects}},
  \code{\link{rfitness}}
}
\examples{

## Generate random fitness for four genes-genotypes
## and plot landscape.

r1 <- rfitness(4)
plot(r1)


## Specify fitness in a matrix, and plot it

m5 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 2, 3, 5.5))
plot(m5)


## Specify fitness with allFitnessEffects, and plot it

fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
plot(evalAllGenotypes(fe))

## same as
plotFitnessLandscape(evalAllGenotypes(fe))

}

\keyword{ hplot }

