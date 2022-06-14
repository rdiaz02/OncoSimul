\name{createUserVars}
\alias{createUserVars}
\alias{createRules}
\alias{adapt_rules_to_cpp}
\alias{transform_rule}
\alias{verify_user_vars}
\alias{verify_rules}
\alias{check_same_name}
\alias{check_double_rule_id}
\alias{check_acttion}


\title{
  Functions that check and create specifications for user variables and rules. 
}

\description{
  This functions check that the user has specified correctly the
  user variables and rules and also makes some modifications in the specification,
  so the "core" of the code that runs the simulation can "understand" them 
  and execute them.
}

\usage{
    createUserVars(userVars)

    createRules(rules, genotFitness, frequencyType = "auto")
}

\arguments{

    \item{userVars}{
        userVars must be a list of lists, where each 
        "sub-list" must have the following fields:

            * name: The name of the variable, must be unique.
            * Value: initial numeric value of the variable.
    }

    \item{rules}{
        rules must be a list of lists, where each 
        "sub-list" must have the following fields:

            * ID: The identifier of the rule, must be unique.
            * Condition: boolean expression that, if true, determines the execution of the rule.
            * Condition: expression thatdetermines the variables that will be modified when the condition is true, 
            it can be arbitrarily complex using other simulation parameters such as N, T and genotype populations and rates.
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
    For \code{createUserVars}, the same list that the user specifies, after checking that all the parameters are correctly specified.
    For \code{createRules} the same list of list that the user specifies, 
    but with the following changes:

    First, it transforms the arguments that refer to the genotipes, for example:
    n_A is the actual population of A in the simulation for a T given. 
    But in the C++ part, "A" receives a Genotype ID, in this case 1, so n_A in the simulation
    is n_1. (For more info run \code{allFitnessEffects} with parameter \code{frequencyDependentFitness = TRUE}, 
    then, check the data that returns, specificly, the field \code{$full_FDF_spec}. 
    There you have more info about those transformations).

    Then, it checks that all fields of the sub-lists are correctly specified.

    Finally, it returns the list of rules with the modifications needed for the code to 
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
    
    # now we specify some user variables 
    userVars <- list(
        list(Name           = "user_var1",
            Value       = 0
        ),
        list(Name           = "user_var2",
            Value       = 3
        ),
        list(Name           = "user_var3",
            Value       = 2.5
        )
    )

    # we call the function to check the specification of the variables
    userVars <- createUserVars(userVars = userVars)

    # we determine the rules that modify the variables
    rules <- list(
        list(ID = "rule_1",
            Condition = "T > 20",
            Action = "user_var_1 = 1"
        ),list(ID = "rule_2",
            Condition = "T > 30",
            Action = "user_var_2 = 2; user_var3 = 2*N"
        ),list(ID = "rule_3",
            Condition = "T > 40",
            Action = "user_var_3 = 3;user_var_2 = n_A*n_B"
        )
    )

    # we call the function to check the specification of the rules
    rules <- createRules(rules = rules, afd3)

    # we run the simulations passing theese lists as arguments
    ep3 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    sampleEvery = 0.01,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            userVars = userVars,
                    rules = rules
                    )

    # you can also make the rules depend on the total population
     rules <- list(
            list(ID = "rule_1",
                Condition = "N > 5000",
                Action = "user_var_1 = 1"
            ),list(ID = "rule_2",
                Condition = "N <= 5000",
                Action = "user_var_1 = 2"
            ),list(ID = "rule_3",
                Condition = "N > 4000",
                Action = "user_var_2 = 1;user_var_3 = 1"
            ),list(ID = "rule_4",
                Condition = "N <= 4000",
                Action = "user_var_2 = 2;user_var_3 = 3"
            )
        )
    

    # or depend on the population of each genotype
    rules <- list(
            list(ID = "rule_1",
                Condition = "n_B > 300",
                Action = "user_var_1 = 1"
            ),list(ID = "rule_2",
                Condition = "n_B > 400",
                Action = "user_var_1 = 2"
            ),list(ID = "rule_3",
                Condition = "n_B <= 300",
                Action = "user_var_1 = 3"
            ),list(ID = "rule_4",
                Condition = "n_B <= 200",
                Action = "user_var_1 = 4"
            )
        )

    # or depend on other previously defined user vars
    rules <- list(
        list(ID = "rule_3",
            Condition = "T > 10",
            Action = "user_var_1 = 1"
        ),list(ID = "rule_1",
            Condition = "user_var_1 = 0",
            Action = "user_var_2 = 1"
        ),list(ID = "rule_2",
            Condition = "user_var_1 = 1",
            Action = "user_var_2 = 2"
        )
    )
    

    # or mix it all together using logic conectors
    rules <- list(
        list(ID = "rule_3",
            Condition = "T > 10 and N < 5000",
            Action = "user_var_1 = 1"
        ),list(ID = "rule_1",
            Condition = "user_var_1 = 0 and n_B > 1000",
            Action = "user_var_2 = 1"
        )
    )
}