\name{oncoSimulIndiv}
\alias{oncoSimulIndiv}
\alias{oncoSimulPop}
\alias{print.oncosimul}
\alias{print.oncosimulpop}
\alias{summary.oncosimul}
\alias{summary.oncosimulpop}
\alias{oncoSimulSample}


\title{
  Simulate tumor progression for one or more individuals, optionally
  returning just a sample in time.
}

\description{
  Simulate tumor progression including possible restrictions in the
  order of driver mutations. Optionally add passenger
  mutations. Simulation is done using the BNB algorithm of Mather et
  al., 2012.
}
\usage{
 oncoSimulIndiv(fp, model = "Exp", 
                numPassengers = 0, mu = 1e-6, muEF = NULL,
                detectionSize = 1e8, detectionDrivers = 4,
                detectionProb = NA,                  
                sampleEvery = ifelse(model \%in\% c("Bozic", "Exp"), 1,
                             0.025),
                initSize = 500, s = 0.1, sh = -1,
                K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                minDetectDrvCloneSz = "auto",
                extraTime = 0,
                finalTime = 0.25 * 25 * 365, onlyCancer = TRUE,
                keepPhylog = FALSE,
                mutationPropGrowth = ifelse(model == "Bozic",
                                                       FALSE, TRUE),
                max.memory = 2000, max.wall.time = 200,
                max.num.tries = 500,
                errorHitWallTime = TRUE,
                errorHitMaxTries = TRUE,
                verbosity = 0,
                initMutant = NULL,
                seed = NULL)

oncoSimulPop(Nindiv, fp, model = "Exp", numPassengers = 0, mu = 1e-6,
                muEF = NULL,
                detectionSize = 1e8, detectionDrivers = 4,
                detectionProb = NA,                  
                sampleEvery = ifelse(model \%in\% c("Bozic", "Exp"), 1,
                             0.025),
                initSize = 500, s = 0.1, sh = -1,
                K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                minDetectDrvCloneSz = "auto",
                extraTime = 0,
                finalTime = 0.25 * 25 * 365, onlyCancer = TRUE,
                keepPhylog = FALSE,
                mutationPropGrowth = ifelse(model == "Bozic",
                                                       FALSE, TRUE),
                max.memory = 2000, max.wall.time = 200,
                max.num.tries = 500,
                errorHitWallTime = TRUE,
                errorHitMaxTries = TRUE,
                initMutant = NULL,
                verbosity = 0,
                mc.cores = detectCores(),
                seed = "auto")


oncoSimulSample(Nindiv,
                fp,
                model = "Exp",
                numPassengers = 0,
                mu = 1e-6,
                muEF = NULL,
                detectionSize = round(runif(Nindiv, 1e5, 1e8)),

                detectionDrivers = {
                                if(inherits(fp, "fitnessEffects")) {
                                    if(length(fp$drv)) {
                                        nd <- (2: round(0.75 * length(fp$drv)))
                                    } else {
                                        nd <- 9e6
                                    }
                                } else {
                                    nd <- (2 : round(0.75 * max(fp)))
                                }
                                if (length(nd) == 1) 
                                    nd <- c(nd, nd)
                                sample(nd, Nindiv,
                                       replace = TRUE)
                            },
                detectionProb = NA,
                sampleEvery = ifelse(model \%in\% c("Bozic", "Exp"), 1,
                    0.025),
                initSize = 500,
                s = 0.1,
                sh = -1,
                K = initSize/(exp(1) - 1),
                minDetectDrvCloneSz = "auto",
                extraTime = 0,
                finalTime = 0.25 * 25 * 365,
                onlyCancer = TRUE, keepPhylog = FALSE,
                mutationPropGrowth = ifelse(model == "Bozic",
                                                       FALSE, TRUE),
                max.memory = 2000,
                max.wall.time.total = 600,
                max.num.tries.total = 500 * Nindiv,
                typeSample = "whole",
                thresholdWhole = 0.5,
                initMutant = NULL,
                verbosity  = 1,
                showProgress = FALSE,
                seed = "auto")


}

\arguments{
  
  \item{Nindiv}{Number of individuals or number of different
    trajectories to simulate.}

  \item{fp}{
    
    Either a poset that specifies the order restrictions (see
    \code{\link{poset}} if you want to use the specification as in
    v.1. Otherwise, a fitnessEffects object (see
    \code{\link{allFitnessEffects}}).

    Other arguments below (s, sh, numPassengers) make sense only if you
    use a poset, as they are included in the fitnessEffects object.

}

\item{model}{
  One of "Bozic", "Exp", "McFarlandLog" (the last one can be abbreviated
  to "McFL"). The default is "Exp".

}

\item{numPassengers}{

  This has no effect if you use the \code{\link{allFitnessEffects}}
  specification.
  
  If you use the specification of v.1., the number of passenger
  genes. Note that using v.1 the total number of genes (drivers plus
  passengers) must be smaller than 64.
  
  All driver genes should be included in the poset (even if they depend
  on no one and no one depends on them), and will be numbered from 1 to
  the total number of driver genes. Thus, passenger genes will be
  numbered from (number of driver genes + 1):(number of drivers + number
  of passengers).

}

\item{mu}{
  Mutation rate. Can be a single value or a named vector. If a single
  value, all genes will have the same mutation rate. If a named vector,
  the entries in the vector specify the gene-specific mutation rate. If
  you pass a vector, it must be named, and it must have entries for all
  the genes in the fitness specification. Passing a vector is only
  available when using fitnessEffects objects for fitness specification.

  See also \code{mutationPropGrowth}.

}

\item{muEF}{

  Mutator effects. A mutatorEffects object as obtained from
  \code{\link{allMutatorEffects}}. This specifies how mutations in
  certain genes change the mutation rate over all the genome. Therefore,
  this allows you to specify mutator phenotypes: models where mutation
  of one (or more) gene(s) leads to an increase in the mutation
  rate. This is only available for version 2 (and above) specifications.
  
  
  All the genes specified in \code{muEF} MUST be included in
  \code{fp}. If you want to have genes that have no direct effect on
  fitness, but that affect mutation rate, you MUST specify them in
  \code{fp}, for instance as \code{noIntGenes} with an effect of 0.
  
  If you use mutator effects you must also use fitnessEffects in
  \code{fp}.

}

\item{detectionSize}{ What is the minimal number of cells for cancer
to be detected.  For \code{oncoSimulSample} this can be a vector.

If set to NA, \code{detectionSize} plays no role in stopping the
simulations.

} 

\item{detectionDrivers}{ The minimal number of drivers (not modules,
drivers, whether or not they are from the same module) present in any
clone for cancer to be detected. For \code{oncoSimulSample} this can be
a vector.

For \code{oncoSimulSample}, if there are drivers (either because you are
using a v.1 object or because you are using a fitnessEffects object with
a \code{drvNames} component ---see \code{\link{allFitnessEffects}}---)
the default is a vector of drivers from a uniform between 2 and 0.75 the
total number of drivers. If there are no drivers (because you are using
a fitnessEffects object without a \code{drvNames}, either because you
specified it explicitly or because all of the genes are in the
noIntGenes component) the simulations should not stop based on the
number of drivers (and, thus, the default is set to 9e6).

If set to NA, \code{detectionDrivers} plays no role in stopping the
simulations.

}

\item{detectionProb}{ Vector of arguments for the mechanism where
  probability of detection depends on size. If \code{NA}, this mechanism
  is not used. If \sQuote{default}, the vector will be populated with
  default values.  Otherwise, a named vector with some of the following
  named elements (see \sQuote{Details}):
  \itemize{

    \item{PDBaseline}{Baseline size subtracted to total population size
    to compute the probability of detection. If not given explicitly,
    the default is \code{1.2 * initSize}.
    }

    \item{p2}{The probability of detection at population size
    \code{n2}. If you specificy \code{p2} you must also specify
    \code{n2} and you must not specify \code{cPDetect}. The fault is 0.1.}

    \item{n2}{The population size at which probability of detection is
      \code{p2}. The default is \code{1.1*initSize}.}

    \item{cPDetect}{The change in probability of detection with size. If
      you specify it, you should not specify either of \code{p2} or
      \code{n2}. See \sQuote{Details}.  }

    \item{checkSizePEvery}{Time between successive checks for the
    probability of exiting as a function of population size. If not
    given explicitly, the default is 20. See \sQuote{Details}.}

}
If you only provide some of the elements (except for the pair \code{p2},
\code{n2}, where you must provide both if you provide any), the rest
will be filled with default values.

This option can not be used with v.1 objects.

}


%% \item{PDBaseline}{For the mechanism where probability of detection
%%   depends on size, the baseline. See \sQuote{Details}.
%% }

%% \item{p2}{For the mechanism where probability of detection
%%   depends on size, the probability of detection at size \code{n2}. See \sQuote{Details}.}

%% \item{n2}{For the mechanism where probability of detection depends on
%%   size, the size at which probability of detection is \code{p2}. See
%%   \sQuote{Details}. }

%% \item{cPDetect}{ For the mechanism where probability of detection
%%   depends on size, the change in probability of detection with size. If
%%   you specify it, you should not specify either of \code{p2} or
%%   \code{n2}. See \sQuote{Details}.  }

%% \item{checkSizePEvery}{ For the mechanism where probability of detection
%%   depends on size, time between successive checks. See \sQuote{Details}.  }

\item{sampleEvery}{
  
  How often the whole population is sampled. This is not the same as the
  interval between successive samples that are kept or stored (for that,
  see \code{keepEvery}).

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
    Only relevant if using a poset as this is included in the
  fitnessEffects object.
}
\item{sh}{
  Selection coefficient for drivers with restrictions not satisfied. A
  value of 0 means there are no penalties for a driver appearing in a
  clone when its restrictions are not satisfied.

  To specify "sh=Inf" (in Diaz-Uriarte, 2015) use sh = -1.

  Only relevant if using a poset as this is included in the
  fitnessEffects object.

}
\item{K}{
  Initial population equilibrium size in the McFarland models.

}

\item{keepEvery}{

  Time interval between successive whole population samples that are
actually stored. This must be larger or equal to \code{sampleEvery}. If
keepEvery is not a multiple integer of sampleEvery, the interval between
successive samples that are stored will be the smallest multiple integer
of \code{sampleEvery} that is larger than or equal to \code{keepEvery}.

  %% keepEvery in use will be the
  %% smallest multiple integer of sampleEvery larger than the specified
  %% keepEvery.

  If you want nice plots, set \code{sampleEvery} and \code{keepEvery} to
  small values (say, 5 or 2). Otherwise, you can use a
  \code{sampleEvery} of 1 but a \code{keepEvery} of 15, so that the
  return objects are not huge and the code runs a lot faster.


  Setting \code{keepEvery = NA} means we only keep the very last
  sample. This is useful if you only care about the final state of the
  simulation, not its complete history.
}

%% \item{endTimeEvery}{
%%   If endTimeEvery is > 0, even if conditions for finishing a simulation
%%   (number of drivers or population size) are met at time \emph{t}, the
%%   simulation will run at least until \emph{t + endTimeEvery} and
%%   conditions will be checked again. Only if conditions for finishing a
%%   simulation are still met, will the simulation end.

%%   The reason for this parameter is to ensure that, say, a clone with a
%%   certain number of drivers that would cause the simulation to end has
%%   not just appeared but then gone extinct shortly after. Beware, though,
%%   that in simulations with very fast growth, setting large endTimeEvery
%%   can result in the simulations taking a long time to finish or hitting
%%   the wall time limit.}

\item{minDetectDrvCloneSz}{

  A value of 0 or larger than 0 (by default equal to
  \code{initSize} in the McFarland model). If larger than 0, when
  checking if we are done with a simulation, we verify that the sum of
  the population sizes of all clones that have a number of mutated
  drivers larger or equal to \code{detectionDrivers} is larger or equal
  to this \code{minDetectDrvCloneSz}.

  The reason for this parameter is to ensure that, say, a clone with a
  certain number of drivers that would cause the simulation to end has
  not just appeared and is present in only one individual that might
  then immediately go extinct. This can be relevant in secenarios such
  as the McFarland model.

  See also \code{extraTime}.
}

\item{extraTime}{
  A value larger than zero waits those many additional time periods
  before exiting after having reached the exit condition (population
  size, number of drivers).

  The reason for this setting is to prevent the McFL models from always
  exiting at a time when one clone is increasing its size quickly (see
  \code{minDetectDrvCloneSz}). By setting an \code{extraTime} larger than 0,
  we can sample at points when we are at the plateau.
  

}

\item{finalTime}{
  What is the maximum number of time units that the simulation can
  run. Set to NA to disable this limit.

}
\item{onlyCancer}{
  Return only simulations that reach cancer?

  If set to TRUE, only simulations that satisfy the
  \code{detectionDrivers} or the \code{detectionSize} requirements or
  that are "detected" because of the \code{detectionProb} mechanism will
  be returned: the simulation will be repeated, within the limits set by
  \code{max.num.tries} and \code{max.wall.time} (and, for
  \code{oncoSimulSample} also \code{max.num.tries.total} and
  \code{max.wall.time.total}), until one which meets the
  \code{detectionDrivers} or \code{detectionSize} or one which is
  detected stochastically under \code{detectionProb} is
  obtained.

  If \code{onlyCancer = FALSE} the simulation is returned regardless of
  final population size or number of drivers in any clone and this
  includes simulations where the population goes extinct.

}

\item{keepPhylog}{
  If TRUE, keep track of when and from which clone each clone is
  created. See also \code{\link{plotClonePhylog}}. 
}



\item{mutationPropGrowth}{
  If TRUE, make mutation rate proportional to growth rate, so clones
  that grow faster also mutate faster. Thus, $mutation_rate = mu *
  birth_rate$. This is a simple way of approximating that mutation happens
  at cell division (it is not strictly making mutation happen at cell division,
  since mutation is not strictly coupled with division).
  Of course, this only makes sense in models where birth rate changes.
}

%% \item{seed}{The seed to use for the C++ random number generator. If not
%%   passed explicitly (the default) then chosen from R.}

\item{initMutant}{For v.2: a string with the mutations of the initial
  mutant, if any. This is the same format as for
  \code{\link{evalGenotype}}. The default (if you pass nothing) is to
  start the simulation from the wildtype genotype with nothing mutated.
  For v.1 we no longer accept initMutant: it will be ignored. }

%% \item{initMutant}{For v.2, a string with the mutations of the initial
%%   mutant, if any. This is the same format as for
%%   \code{\link{evalGenotype}}. For v.1, the single mutation of the
%%   initial clone for the simulations. The default (if you pass nothing)
%%   is to start the simulation from the wildtype genotype with nothing
%%   mutated.}


\item{max.num.tries}{Only applies when \code{onlyCancer = TRUE}. What is
  the maximum number of times, for an individual simulation, we can
  repeat the simulation for it to reach cancer? There are certain
  parameter settings where reaching cancer is extremely unlikely and you
  might not want to run forever in those cases.  }


\item{max.num.tries.total}{Only applies when \code{onlyCancer = TRUE}
  and for \code{oncoSimulSample}. What is the maximum number of times,
  over all simulations for all individuals in a population sample, that
  we can repeat the simulations so that cancer is reached for all
  individuals?  The idea is to set a limit on the average minimal
  probability of reaching cancer for a set of simulations to be
  accepted.}

\item{max.wall.time}{
  Maximum wall time for each individual simulation run. If the
  simulation is not done in this time, it is aborted.

}

\item{max.wall.time.total}{

  Maximum wall time for all the simulations (when using
  \code{oncoSimulSample}), in seconds. If the simulation is not
  completed in this time, it is aborted. To prevent problems from a
  single individual simulation going wild, this limit is also enforced
  per simulation (so the run can be aborted directly from C++).

}

\item{errorHitMaxTries}{ If TRUE (the default) a simulation that reaches
  the maximum number of repetitions allowed is considered not to have
  succesfully finished and, thus, an error, and no output from it will
  be reported. This is often what you want.%%  But if you want to see
  %% what is happening for some set of parameters that seem to never
  %% finish, set this to FALSE, and check what the output looks like.
  See  \code{Details}.
}

\item{errorHitWallTime}{
  If TRUE (the default) a simulation that reaches the maximum wall time
  is considered not to have succesfully finished and, thus, an error,
  and no output from it will be reported. This is often what you
  want. %% But if you want to see what is happening for some set of
  %% parameters that seem to never finish, set this to FALSE, and check
  %% what the output looks like.
  See \code{Details}.
}



\item{max.memory}{
  The largest size (in MB) of the matrix of Populations by Time. If it
  creating it would use more than this amount of memory, it is not
  created. This prevents you from accidentally passing parameters that
  will return an enormous object.

}






\item{verbosity}{ If 0, run silently. Iincreasing values of verbosity
  provide progressively more information about intermediate steps,
  possible numerical notes/warnings from the C++ code, etc. Values less
  than 0 supress some default notes: use with care.

}

  \item{typeSample}{
    "singleCell" (or "single") for single cell sampling, where the
    probability of sampling a cell (a clone) is directly proportional to
    its population size.  "wholeTumor" (or "whole") for whole tumor
    sampling (i.e., this is similar to a biopsy being the entire
    tumor). See \code{\link{samplePop}}.
  }

  \item{thresholdWhole}{
    In whole tumor sampling, whether a gene is detected as mutated depends
    on thresholdWhole: a gene is considered mutated if it is altered in at
    least thresholdWhole proportion of the cells in that individual. See \code{\link{samplePop}}.
  }
  

\item{mc.cores}{Number of cores to use when simulating more than one
  individual (i.e., when calling oncoSimulPop).}

\item{showProgress}{If TRUE, provide information, during exection,
of the individual done, and the number of attempts and time used.}

\item{seed}{The seed for the C++ PRNG. You can pass a value. If you set
  it to NULL, then a seed will be generated in R and passed to C++. If
  you set it to "auto", then if you are using v.1, the behavior is the
  same as if you set it to NULL (a seed will be generated in R and
  passed to C++) but if you are using v.2, a random seed will be
  produced in C++. %% using the randutils code from M.E.O'Neill.
  If you   need reproducibility, either pass a value or set it to NULL (setting
  it to NULL will make the C++ seed reproducible if you use the same
  seed in R via \code{set.seed}). However, even using the same value of
  \code{seed} is unlikely to give the exact same results between
  platforms and compilers.  Moreover, note that the defaults for
  \code{seed} are not the same in \code{oncoSimulIndiv},
  \code{oncoSimulPop} and \code{oncoSimulSample}.

  When using oncoSimulPop, if you want reproducibility, you might want
  to, in addition to setting \code{seed = NULL}, also do
  \code{RNGkind("L'Ecuyer-CMRG")} as we use
  \code{\link[parallel]{mclapply}}; look at the vignette of
  \pkg{parallel}.

}



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


  Changing options \code{errorHitMaxTries} and \code{errorHitWallTime}
  can be useful when conducting many simulations, as in the call to
  \code{oncoSimulPop}: setting them to TRUE means nothing is recorded
  for those simulations where ending conditions are not reached but
  setting them to FALSE would allow you to record the output; this would
  potentially result in a mixture where some simulations would not have
  reached the ending condition, but this might sometimes be what you
  want. Note, however, that \code{oncoSimulSample} always has both them
  to TRUE, as it could not be otherwise.



  \code{GenotypesWDistinctOrderEff} provides the information about
  order effects that is missing from \code{Genotypes}. When there are
  order effects, the \code{Genotypes} matrix can contain genotypes
  that are not distinguishable. Suppose there are two genes, the first
  and the second. In the \code{Genotype} output you can get two
  columns where there is a 1 in both genes: those two columns
  correspond to the two possible orders (first gene mutated first, or
  first gene mutated after the
  second). \code{GenotypesWDistinctOrderEff} disambiguates this. The
  same is done by \code{GenotypesLabels}; this is easier to decode for
  a human (a string of gene labels) but a little bit harder to parse
  automatically.  Note that when you use the default print method for
  this object, you get, among others, a two-column display with the
  \code{GenotypeLabels} information. When order matters, a genotype
  shown as ``x > y _ z'' means that a mutation in ``x'' happened
  before a mutation in ``y''; there is also a mutation in ``z'' (which
  could have happened before or after either of ``x'' or ``y''), but
  ``z'' is a gene for which order does not matter. When order does not
  matter, a comma "," separates the identifiers of mutated genes.


  Detection of cancer can be a deterministic process, where cancer is
  always detected (and, thus, simulation ended) when certain conditions
  are met (\code{detectionSize},
  \code{detectionDrivers}). Alternatively, it can be stochastic process
  where probability of detection depends on size. Every so often (see
  below) we assess population size, and detect cancer or not
  probabilistically (comparing the probability of detection for that
  size with a random uniform number). Probability of detection changes
  with population size according to the function
  
  \deqn{ 1 - e^{-cPDetect (populationsize - PDBaseline)}}.

  You can pass \code{cPDetect} manually (you will need to set \code{n2}
  and \code{p2} to NA). However, it might be more intuitive to specify
  the pair \code{n2}, \code{p2}, such that the probability of detection
  is \eqn{p2} for population size \eqn{n2} (and from that pair we solve
  for the value of \code{cPDetect}). How often do we check? That is
  controlled by \code{checkSizePEvery}, the (minimal) time between
  successive checks (from among the sampling times given by
  \code{sampleEvery}: the interval between successive assessments will
  be the smallest multiple integer of \code{sampleEvery} that is
  larger than \code{checkSizePEvery} ---see vignette for details).
  \code{checkSizePEvery} has, by default, a different (and much larger)
  value than \code{sampleEvery} both to allow to examine the effects of
  sampling, and to avoid many costly random number generations.

  Please note that \code{detectionProb} is NOT available with version 1
  objects.

  %% : we assess probability of exiting at every sampling time that (the
  %% closest possible to) a multiple of \code{checkSizePEvery}.

}
\value{
  
  For \code{oncoSimulIndiv} a list, of class "oncosimul", with the
  following components:
  
    \item{pops.by.time }{A matrix of the population sizes of the clones,
    with clones in columns and time in row. Not all clones are shown here,
    only those that were present in at least on of the keepEvery samples.}
  
  \item{NumClones   }{Total number of clones in the above matrix. This
  is not the total number of distinct clones that have appeared over all
  simulations (which is likely to be larger or much larger).}

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
    simulation error when using the McFarland model as well as other
    statistics. For the McFarland model, the relevant value is errorMF,
    which is -99 unless in the McFarland model. For the McFarland model
    it is the largest difference of successive death rates. The entries
    names minDMratio and minBMratio are the smallest ratio, over all
    simulations, of death rate to mutation rate or birth rate to
    mutation rate. The BNB algorithm thrives when those are large.}
  

  For \code{oncoSimulPop} a list of length \code{Nindiv}, and of class
  \code{"oncosimulpop"}, where each element of the list is itself a
  list, of class \code{oncosimul}, with components as described above.


  In v.2, the output is of both class "oncosimul" and "oncosimul2". The
  oncoSimulIndiv return object differs in

  \item{GenotypesWDistinctOrderEff}{
    A list
  of vectors, where each vector corresponds to a genotype in the
  \code{Genotypes}, showing (where it matters) the order of mutations.
    Each vector shows the genotypes, with the numeric
  codes, showing explicitly the order when it matters. So if you have
  genes 1, 2, 7 for which order relationships are given, and genes 3, 4,
  5, 6 for which other interactions exist, any mutations in 1, 2, 7 are
  shown first, and in the order they occurred, before showing the rest
  of the mutations.  See details.}
  
 


  \item{GenotypesLabels}{The genotypes, as character vectors with the
  original labels provided (i.e., not the integer codes). As before,
  mutated genes, for those where order matters, come first, and are
  separated by the rest by a "_". See details.}

  \item{OccurringDrivers}{This is the same as in v.1, but we use the
  labels, not the numeric id codes. Of course, if you entered integers
  as labels for the genes, you will see numbers (however, as a character
  string).}

%% Fixme: mention probCancer is not a binomial, but a negative binomial.


}
\references{

  Bozic, I., et al., (2010). Accumulation of driver and passenger
  mutations during tumor progression. \emph{ Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{107},
  18545--18550.
  
  Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
  accumulation of mutations during tumor progression: effects of
  passengers, evolutionary models, and sampling
  \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}

  Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
  Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
  
  McFarland, C.~D. et al. (2013). Impact of deleterious passenger
  mutations on cancer progression.  \emph{Proceedings of the National
  Academy of Sciences of the United States of America\/}, \bold{110}(8),
  2910--5.

  Mather, W.~H., Hasty, J., and Tsimring, L.~S. (2012). Fast stochastic
  algorithm for simulating evolutionary population dynamics.
  \emph{Bioinformatics (Oxford, England)\/}, \bold{28}(9), 1230--1238.

  

}


\note{Please, note that the meaning of the fitness effects in the
  McFarland model is not the same as in the original paper; the fitness
  coefficients are transformed to allow for a simpler fitness function
  as a product of terms. This differs with respect to v.1. See the
  vignette for details.}

\author{Ramon Diaz-Uriarte}


\seealso{
  \code{\link{plot.oncosimul}}, \code{\link{examplePosets}},
  \code{\link{samplePop}}, \code{\link{allFitnessEffects}}
  

}
\examples{

#################################
#####
##### Examples using v.1
#####
#################################


## use poset p701
data(examplePosets)
p701 <- examplePosets[["p701"]]

## Exp Model

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
                     keepEvery = 10,
                     onlyCancer = FALSE)
plot(m1, addtot = TRUE, log = "")




## Simulating 4 individual trajectories
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).


p1 <- oncoSimulPop(4, p701,
                   keepEvery = 10,
                   mc.cores = 2)
summary(p1)
samplePop(p1)

p2 <- oncoSimulSample(4, p701)


#########################################
######
###### Examples using v.2:
######
#########################################



#### A model similar to the one in McFarland. We use 2070 genes.

set.seed(456)
nd <- 70  
np <- 2000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drv = seq.int(nd))
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.02,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         onlyCancer = FALSE)
plot(mcf1s, addtot = TRUE, lwdClone = 0.6, log = "")
summary(mcf1s)
plot(mcf1s)


#### Order effects with modules, and 5 genes without interactions
#### with fitness effects from an exponential distribution

oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
oiI1 <- oncoSimulIndiv(oi, model = "Exp")
oiI1$GenotypesLabels
oiI1   ## note the order and separation by "_"


oiP1 <- oncoSimulPop(2, oi,
                     keepEvery = 10,
                     mc.cores = 2)
summary(oiP1) 


## Even if order exists, this cannot reflect it;
## G1 to G10 are d1, d2, f1..,f3, and the 5 genes without
## interaction
samplePop(oiP1)


oiS1 <- oncoSimulSample(2, oi)

## The output contains only the summary of the runs AND
## the sample:
oiS1 

## And their sizes do differ
object.size(oiS1)
object.size(oiP1)



######## Using a poset for pancreatic cancer from Gerstung et al.
###      (s and sh are made up for the example; only the structure
###       and names come from Gerstung et al.)


pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
plot(pancr)

### Use an exponential growth model

pancr1 <- oncoSimulIndiv(pancr, model = "Exp")
pancr1
summary(pancr1)
plot(pancr1)
pancr1$GenotypesLabels


## Pop and Sample
pancrPop <- oncoSimulPop(4, pancr,
                        keepEvery = 10,
                       mc.cores = 2)
summary(pancrPop)
pancrSPop <- samplePop(pancrPop)
pancrSPop

pancrSamp <- oncoSimulSample(2, pancr)
pancrSamp



## Using gene-specific mutation rates
muv <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
ni <- rep(0.01, 5)
names(ni) <- names(muv)
femuv <- allFitnessEffects(noIntGenes = ni)
oncoSimulIndiv(femuv, mu = muv)
                      

## Reinitialize the RNG
set.seed(NULL)                      

}
\keyword{misc}
\keyword{iteration}



