\name{createInterventions}
\alias{createInterventions}
\alias{adapt_interventions_to_cpp}
\alias{transform_intervention}
\alias{verify_interventions}
\alias{check_double_id}
\alias{check_what_happens}


\title{
  Function that checks and creates an specification for interventions. 
}

\description{
  This functions checks that the user has specified correctly the
  interventions and also makes some modifications in the specification,
  so the "core" of the code that runs the simulation can "understand" them 
  and execute them.
}

\usage{
    createInterventions(interventions, genotFitness, frequencyType = "auto")
}

\arguments{

    \item{interventions}{
        Interventions must be a list of lists, where each 
        "sub-list" must have the following fields:

            * ID: The identifier of the intervention, must be unique.
            * Trigger: The situation in the simulation that triggers/activates 
                        the intervention.
            * What Happens: "What happens" in the simulation. Basically, once the trigger 
                            is satisfied, this defines how the population is affected by 
                            the intervention.
            * Periodicity: Defines the periodicity of the intervention.
            * Repetitions: Defines the maximum repetitions of each intervention in the simulation.
    }

    \item{genotFitness}{
        Object that \code{allFitnessEffects} returns, it is necessary to call that 
        function before creating interventions. Also, when calling \code{allFitnessEffects}
        frequencyDependentFitness must be TRUE.
    }

    \item{frequencyType}{
        If you want to specify the frequency type of the simulation, by default is set to "auto"
    }
}

\details{
    N/A
}

\value{
    Returns the same list of list that the user specifies, 
    but with the following changes:

    First, it transforms the arguments that refer to the genotipes, for example:
    n_A is the actual population of A in the simulation for a T given. 
    But in the C++ part, "A" receives a Genotype ID, in this case 1, so n_A in the simulation
    is n_1. (For more info run \code{allFitnessEffects} with parameter \code{frequencyDependentFitness = TRUE}, 
    then, check the data that returns, specificly, the field \code{$full_FDF_spec}. 
    There you have more info about those transformations).

    Then, it checks that all fields of the sub-lists are correctly specified.

    Finally, it returns the list of interventions with the modifications needed for the code to 
    interpret it correctly. 
}

\examples{

    #first we create and the populations to simulate.
    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1.001 + (0*n_A)",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")
    
    # now we especify intervention to drastically reduce A population
    # depending on the T of the simulation
    list_of_interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "(T >= 5)",
            WhatHappens   = "n_A = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    # we transform the intervention to somthing the simulation can process
    final_interventions <- createInterventions(interventions = list_of_interventions, afd3)

    # we run the simulations passing the interventions as an argument
    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    sampleEvery = 0.01,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = final_interventions
                    )

    # you can also make the intervention depend on the total population
    list_of_interventions1 <- list(
        list(ID           = "intOverTotPop",
            Trigger       = "(N >= 5000)",
            WhatHappens   = "N = N * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    # or depend over the population of a genotype 
    list_of_interventions2 <- list(
        list(ID           = "intOverTotPop",
            Trigger       = "(n_A >= 5000)",
            WhatHappens   = "n_ = n_B * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    # or mix it all together using logic conectors
    list_of_interventions2 <- list(
        list(ID           = "intOverTotPop",
            Trigger       = "((n_A >= 200) and (N >= 2000)) or (T > 20)",
            WhatHappens   = "n_ = n_B * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )
}