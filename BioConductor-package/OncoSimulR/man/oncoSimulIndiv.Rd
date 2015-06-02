\name{oncoSimulIndiv}
\alias{oncoSimulIndiv}
\alias{oncoSimulPop}
\alias{print.oncosimul}
\alias{print.oncosimulpop}
\alias{summary.oncosimul}
\alias{summary.oncosimulpop}



\title{
  Simulate tumor progression for one or more individuals.
}
\description{
  Simulate tumor progression including possible restrictions in the
  order of driver mutations. Optionally add passenger
  mutations. Simulation is done using the BNB algorithm of Mather et
  al., 2012.
}
\usage{
 oncoSimulIndiv(poset, model = "Bozic", numPassengers = 30, mu = 1e-6,
                detectionSize = 1e8, detectionDrivers = 4,
                sampleEvery = ifelse(model \%in\% c("Bozic", "Exp"), 1,
                             0.025),
                initSize = 500, s = 0.1, sh = -1,
                K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                endTimeEvery = ifelse(model \%in\% c("Bozic", "Exp"), -9,
                                               5 * sampleEvery),
                finalTime = 0.25 * 25 * 365, onlyCancer = TRUE,
                max.memory = 2000, max.wall.time = 200,
                verbosity = 0)

 oncoSimulPop(Nindiv, poset, model = "Bozic", numPassengers = 30, mu = 1e-6,
                detectionSize = 1e8, detectionDrivers = 4,
                sampleEvery = ifelse(model \%in\% c("Bozic", "Exp"), 1,
                             0.025),
                initSize = 500, s = 0.1, sh = -1,
                K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                endTimeEvery = ifelse(model \%in\% c("Bozic", "Exp"), -9,
                                               5 * sampleEvery),
                finalTime = 0.25 * 25 * 365, onlyCancer = TRUE,
                max.memory = 2000, max.wall.time = 200,
                verbosity = 0, mc.cores = detectCores())


}

\arguments{
  
  \item{Nindiv}{Number of individuals or number of different
    trajectories to simulate.}
  \item{poset}{
    
    The poset that specifies the order restrictions. See \code{\link{poset}}.

}
\item{model}{
  One of "Bozic", "Exp", "McFarlandLog" (the last one can be abbreviated
  to "McFL").

}
\item{numPassengers}{
  
  The number of passenger genes.The total number of genes (drivers plus
  passengers) must be smaller than 64.
  
  All driver genes should be included in the poset (even if they depend
  on no one and no one depends on them), and will be numbered from 1 to
  the total number of driver genes. Thus, passenger genes will be
  numbered from (number of driver genes + 1):(number of drivers + number
  of passengers).

}
\item{mu}{
  Mutation rate.

}
\item{detectionSize}{
  What is the minimal number of cells for cancer to be detected.

}
\item{detectionDrivers}{
  The minimal number of drivers present in any clone for cancer to be detected.

}
\item{sampleEvery}{
  
  How often the whole population is sampled. This is not the same as the
  interval between successive samples that keep stored (for that, see
  \code{keepEvery}).

  For very fast growing clones, you might need to have a small value
  here to minimize possible numerical problems (such as huge increase in
  population size between two successive samples that can then lead to
  problems for random number generators). Likewise, for models with
  density dependence (such as McF) this value should be very small.

}
\item{initSize}{
  Initial population size.

}
\item{s}{
  Selection coefficient for drivers.

}
\item{sh}{
  Selection coefficient for drivers with restrictions not satisfied. A
  value of 0 means there are no penalties for a driver appearing in a
  clone when its restrictions are not satisfied.

  To specify "sh=Inf" (in Diaz-Uriarte, 2014) use sh = -1.

}
\item{K}{
  Initial population equilibrium size in the McFarland models.

}
\item{keepEvery}{
  Time interval between successive whole population samples that are
  actually stored. This must be larger or equal to sampleEvery. If keepEvery is
  not a multiple integer of sampleEvery, the keepEvery in use will be the
  smallest multiple integer of keepEvery larger than the specified
  keepEvery.

  If you want nice plots, set \code{sampleEvery} and \code{keepEvery} to
  small values (say, 1 or 0.5). Otherwise, you can use a
  \code{sampleEvery} of 1 but a \code{keepEvery} of 15, so that the
  return objects are not huge.

}

\item{endTimeEvery}{
  If endTimeEvery is > 0, even if conditions for finishing a simulation
  (number of drivers or population size) are met at time \emph{t}, the
  simulation will run at least until \emph{t + endTimeEvery} and
  conditions will be checked again. Only if conditions for finishing a
  simulation are still met, will the simulation end.

  The reason for this parameter is to ensure that, say, a clone with a
  certain number of drivers that would cause the simulation to end has
  not just appeared but then gone extinct shortly after. Beware, though,
  that in simulations with very fast growth, setting large endTimeEvery
  can result in the simulations taking a long time to finish or hitting
  the wall time limit.}



\item{finalTime}{
  What is the maximum number of time units that the simulation can run.

}
\item{onlyCancer}{
  Return only simulations that reach cancer?

  If set to TRUE, only simulations that satisfy the detectionDrivers or
  the detectionSize will be returned (i.e., the simulation will be
  repeated until one which meets that cancer requirements is
  met). Otherwise, the simulation is returned regardless of final
  population size or number of drivers in any clone and this includes
  simulations where the population goes extinct.

}
\item{max.memory}{
  The largest size (in MB) of the matrix of Populations by Time. If it
  creating it would use more than this amount of memory, it is not
  created. This prevents you from accidentally passing parameters that
  will return an enormous object.

}
\item{max.wall.time}{
  Maximum wall time for each individual simulation run. If the
  simulation is not done in this time, it is aborted.

}
\item{verbosity}{
  If 0, run as silently as possible. Otherwise, increasing values of
  verbosity provide progressively more information about intermediate
  steps, possible numerical notes/warnings from the C++ code, etc.

}

\item{mc.cores}{Number of cores to use when simulating more than one
  individual (i.e., when calling oncoSimulPop).}

}
\details{

  The basic simulation algorithm implemented is the BNB one of Mather et
  al., 2012, where I have added modifications to fitness based on the
  restrictions in the order of mutations.

  Full details about the algorithm are provided in Mather et al.,
  2012. The evolutionary models, including references, and the rest of
  the parameters are explained in Diaz-Uriarte, 2014, especially in the
  Supplementary Material. The model called "Bozic" is based on Bozic et
  al., 2010, and the model called "McFarland" in McFarland et al., 2013.


  oncoSimulPop simply calls oncoSimulIndiv multiple times. When run on
  POSIX systems, it can use multiple cores (via mclapply).


  The \code{summary} methods for these classes return some of the return
  values (see next) as a one-row (for class oncosimul) or multiple row
  (for class oncosimulpop) data frame.  The \code{print} methods for
  these classes simply print the summary.

}
\value{
  
  For \code{oncoSimulIndiv} a list, of class "oncosimul", with the
  following components:
  
    \item{pops.by.time }{A matrix of the population sizes of the clones,
    with clones in columns and time in row. Not all clones are shown here,
    only those that were present in at least on of the keepEvery samples.}
  
  \item{NumClones   }{Total number of clones in the above matrix.  }

  \item{TotalPopSize}{Total population size at the end.}

  \item{Genotypes}{A matrix of genotypes. For each of the clones in the
    pops.by.time matrix, its genotype, with a 0 if the gene is not
    mutated and a 1 if it is mutated.}
  
  \item{MaxNumDrivers}{The largest number of mutated driver genes ever
    seen in the simulation in any clone.}

  \item{MaxDriversLast}{The largest number of mutated drivers in any
    clone     at the   end of the simulation.}

  \item{NumDriversLargestPop}{The number of mutated driver genes in the
    clone with largest population size. }

  \item{LargestClone}{Population size of the clone with largest number
    of population size.}

  \item{PropLargestPopLast}{Ratio of LargestClone/TotalPopSize}

  \item{FinalTime}{The time (in time units) at the end of the
    simulation.}

  \item{NumIter}{The number of iterations of the BNB algorithm.}

  \item{HittedWallTime}{TRUE if we reached the limit of max.wall.time. FALSE
    otherwise.}

  \item{TotalPresentDrivers}{The total number of mutated driver genes,
    whether or not in the same clone. The number of elements in
    \code{OccurringDrivers}, below.}

  \item{CountByDriver}{A vector of length number of drivers, with the
    count of the number of clones that have that driver mutated.}

  \item{OccurringDrivers}{The actual number of drivers mutated.}

  \item{PerSampleStats}{A 5 column matrix with a row for each sampling
    period. The columns are: total population size, population size of the
    largest clone, the ratio of the two, the largest number of drivers in
    any clone, and the number of drivers in the clone with the largest
    population size.}
  
  \item{other}{A list that contains statistics for an estimate of the
    simulation error when using the McFarland model. The relevant value
    is errorMF, which is -99 unless in the McFarland model. For the
    McFarland model it is the largest difference of successive death
    rates.}
  

  For \code{oncoSimulPop} a list of length \code{Nindiv}, and of class
  \code{"oncosimulpop"}, where each element of the list is itself a
  list, of class \code{oncosimul}, with components as described above.
  

}
\references{

  Bozic, I., et al., (2010). Accumulation of driver and passenger
  mutations during tumor progression. \emph{ Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{107},
  18545--18550.
  
  Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
             accumulation of mutations during tumor progression: effects
             of passengers, evolutionary models, and sampling, \emph{BMC
             Bioinformatics\/},
             \bold{16}(41). \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}.

  Inferring restrictions in the temporal order
  of mutations during tumor progression: effects of passenger mutations,
  evolutionary models, and
  sampling. \url{http://dx.doi.org/10.1101/005587}
  
  McFarland, C.~D. et al. (2013). Impact of deleterious passenger
  mutations on cancer progression.  \emph{Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{110}(8),
  2910--5.

  Mather, W.~H., Hasty, J., and Tsimring, L.~S. (2012). Fast stochastic
  algorithm for simulating evolutionary population dynamics.
  \emph{Bioinformatics (Oxford, England)\/}, \bold{28}(9), 1230--1238.

  

}

\author{Ramon Diaz-Uriarte}


\seealso{
  \code{\link{plot.oncosimul}}, \code{\link{examplePosets}},
  \code{\link{samplePop}}
  

}
\examples{
## use poset p701
data(examplePosets)
p701 <- examplePosets[["p701"]]

## Bozic Model

b1 <- oncoSimulIndiv(p701)
summary(b1)

plot(b1, addtot = TRUE)

## McFarland; use a small sampleEvery, but also a reasonable
##   keepEvery.
## We also modify mutation rate to values similar to those in the
##   original paper.
## Note that detectionSize will play no role
## finalTime is large, since this is a slower process
## initSize is set to 4000 so the default K is larger and we are likely
## to reach cancer. Alternatively, set K = 2000.

m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 5)
plot(m1, addtot = TRUE, log = "")




## Simulating 4 individual trajectories
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).


p1 <- oncoSimulPop(4, p701,
                   keepEvery = 10,
                   mc.cores = 2)


}
\keyword{misc}
\keyword{iteration}



