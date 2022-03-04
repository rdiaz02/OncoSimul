inittime <- Sys.time()
cat(paste("\n Starting user variable tests", date(), "\n"))
############################################################################################################################
############################################################################################################################
############################################################################################################################

test_that("1. A user varibale is created correctly",{
    userVars <- list(
        list(Name = "user_var_1",
            Value = 0
        )
    )

    testthat::expect_output(createUserVars(userVars), "Checking user variable: user_var_1" )
    userVars <- createUserVars(userVars)

    # we check that the Name and value atributesare correct
    testthat::expect_equal(userVars[[1]]$Name, "user_var_1")
    testthat::expect_equal(userVars[[1]]$Value, 0)
})

test_that("2. Two user variables cannot have the same Name (check_double_name)", {
    userVars <- list(
        list(Name = "user_var_1",
             Value = 0
        ),list(Name = "user_var_1",
               Value = 2
        )   
    )
    testthat::expect_error(createUserVars(userVars), "Check the user variables, there are 2 or more that have same Names")    
})

test_that("3. A rule is created correctly",{
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B > 0)",
                    ".9 + 0.4 * (n_A > 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")
    
    rules <- list(
        list(ID = "rule_1",
            Condition = "T = 50",
            Action = "user_var_1 = 2"
        )
    )

    testthat::expect_output(createRules(rules, afd3), "Checking rule: rule_1" )
    rules <- createRules(rules, afd3)

    # we check that the transformation of the Action and Condition atributes is correct
    testthat::expect_equal(rules[[1]]$Condition, "T = 50")
    testthat::expect_equal(rules[[1]]$Action, "user_var_1 = 2")
})

test_that("4. Two rules cannot have the same ID (check_double_id)", {
    rules <- list(
        list(ID = "rule_1",
            Condition = "n_B >= 5",
            Action = "user_var_1 = 1"
        ),list(ID = "rule_1",
            Condition = "n_A >= 5",
            Action = "user_var_1 = 2"
        )
    )
    testthat::expect_error(createRules(rules, afd3), "Check the rules, there are 2 or more that have same IDs")    
})

# test_that("5. Rules change user vars corectly (depending on T)", {
#     userVars <- list(
#         list(Name = "user_var_1",
#             Value = 0
#         ),
#         list(Name = "user_var_2",
#             Value = 0
#         ),
#         list(Name = "user_var_3",
#             Value = 0
#         )
#     )

#     userVars <- createUserVars(userVars)

#     rules <- list(
#         list(ID = "rule_1",
#             Condition = "T > 20",
#             Action = "user_var_1 = 1"
#         ),list(ID = "rule_2",
#             Condition = "T > 30",
#             Action = "user_var_2 = 2"
#         ),list(ID = "rule_3",
#             Condition = "T > 40",
#             Action = "user_var_2 = 1"
#         )
#     )

#     rules <- createRules(rules, afd3)

#     sfd3 <- oncoSimulIndiv( afd3,
#                             model = "McFLD",
#                             onlyCancer = FALSE,
#                             finalTime = 100,
#                             mu = 1e-4,
#                             initSize = 5000,
#                             sampleEvery = 0.001,
#                             userVars = userVars,
#                             rules = rules)

#     testthat::expect_equal(interventions[[1]]$WhatHappens, 1)
#     testthat::expect_equal(interventions[[1]]$WhatHappens, 2)
#     testthat::expect_equal(interventions[[1]]$WhatHappens, 3)
# })

cat(paste("\n Ending interventions tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)