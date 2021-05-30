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
        If you want to specify the frequency type of the simulation
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
}