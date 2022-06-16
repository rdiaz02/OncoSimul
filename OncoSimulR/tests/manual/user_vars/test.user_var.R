cat(paste("\n Starting user vars tests at", date()))
############################################################################################################################
############################################################################################################################
############################################################################################################################

## 5. Rules change user vars corectly (depending on T)
gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c("1",
                "1 + 0.2 * (n_B > 0)",
                ".9 + 0.4 * (n_A > 0)"
                ))
afd3 <- allFitnessEffects(genotFitness = gffd3,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "abs")
userVars <- list(
    list(Name = "user_var_1",
        Value = 0
    ),
    list(Name = "user_var_2",
        Value = 0
    ),
    list(Name = "user_var_3",
        Value = 0
    )
)

userVars <- createUserVars(userVars)

rules <- list(
    list(ID = "rule_1",
        Condition = "T > 20",
        Action = "user_var_1 = 1"
    ),list(ID = "rule_2",
        Condition = "T > 30",
        Action = "user_var_2 = 2"
    ),list(ID = "rule_3",
        Condition = "T > 40",
        Action = "user_var_3 = 3;user_var_2 = 1"
    )
)

rules <- createRules(rules, afd3)

sfd3 <- oncoSimulIndiv( afd3,
                        model = "McFLD",
                        onlyCancer = FALSE,
                        finalTime = 100,
                        mu = 1e-4,
                        initSize = 5000,
                        sampleEvery = 0.001,
                        userVars = userVars,
                        rules = rules)

for(line in sfd3$other$userVarValues){
    if(line[4] < 20){
        testthat::expect_equal(line[1], 0)
        testthat::expect_equal(line[2], 0)
        testthat::expect_equal(line[3], 0)
    }else if(line[4] < 30){
        testthat::expect_equal(line[1], 1)
        testthat::expect_equal(line[2], 0)
        testthat::expect_equal(line[3], 0)
    }else if(line[4] < 40){
        testthat::expect_equal(line[1], 1)
        testthat::expect_equal(line[2], 2)
        testthat::expect_equal(line[3], 0)
    }else{
        testthat::expect_equal(line[1], 1)
        testthat::expect_equal(line[2], 1)
        testthat::expect_equal(line[3], 3)
    }
}

# 6. Rules change user vars corectly (depending on N)
gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c("1",
                "1 + 0.2 * (n_B > 0)",
                ".9 + 0.4 * (n_A > 0)"
                ))
afd3 <- allFitnessEffects(genotFitness = gffd3,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "abs")

userVars <- list(
    list(Name = "user_var_1",
        Value = 0
    ),
    list(Name = "user_var_2",
        Value = 0
    ),
    list(Name = "user_var_3",
        Value = 0
    )
)

userVars <- createUserVars(userVars)

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

rules <- createRules(rules, afd3)

sfd3 <- oncoSimulIndiv( afd3,
                        model = "McFLD",
                        onlyCancer = FALSE,
                        finalTime = 100,
                        mu = 1e-4,
                        initSize = 5000,
                        sampleEvery = 0.001,
                        userVars = userVars,
                        rules = rules)

for(i in 2:nrow(sfd3$pops.by.time)){
    line <- sfd3$pops.by.time[i,]
    population <- line[2] + line[3] + line[4]
    time <- line[1]
    varIndex <- 0
    for(aux in sfd3$other$userVarValues){
        varIndex <- varIndex + 1
        if(aux[4] == time){
            break
        }
    }
    vars <- sfd3$other$userVarValues[[varIndex]]
    if(population <= 4000){
        testthat::expect_equal(vars[1], 2)
        testthat::expect_equal(vars[2], 2)
        testthat::expect_equal(vars[3], 3)
    }else if(population <= 5000){
        testthat::expect_equal(vars[1], 2)
        testthat::expect_equal(vars[2], 1)
        testthat::expect_equal(vars[3], 1)
    }else{
        testthat::expect_equal(vars[1], 1)
        testthat::expect_equal(vars[2], 1)
        testthat::expect_equal(vars[3], 1)
    }
}

# 7. Rules change user vars corectly (depending on n_x)
gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c("1",
                "1 + 0.2 * (n_B > 0)",
                ".9 + 0.4 * (n_A > 0)"
                ))
afd3 <- allFitnessEffects(genotFitness = gffd3,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "abs")
userVars <- list(
    list(Name = "user_var_1",
        Value = 0
    )
)

userVars <- createUserVars(userVars)

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

rules <- createRules(rules, afd3)

sfd3 <- oncoSimulIndiv( afd3,
                        model = "McFLD",
                        onlyCancer = FALSE,
                        finalTime = 100,
                        mu = 1e-4,
                        initSize = 5000,
                        sampleEvery = 0.001,
                        userVars = userVars,
                        rules = rules)

for(i in 2:nrow(sfd3$pops.by.time)){
    line <- sfd3$pops.by.time[i,]
    population <- line[4]
    time <- line[1]
    varIndex <- 0
    for(aux in sfd3$other$userVarValues){
        varIndex <- varIndex + 1
        if(aux[2] == time){
            break
        }
    }
    vars <- sfd3$other$userVarValues[[varIndex]]
    if(population <= 200){
        testthat::expect_equal(vars[1], 4)
    }else if(population <= 300){
        testthat::expect_equal(vars[1], 3)
    }else if(population <= 400){
        testthat::expect_equal(vars[1], 1)
    }else{
        testthat::expect_equal(vars[1], 2)
    }
}

# 8. Rules change user vars corectly (depending on other user vars)
gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                Fitness = c("1",
                "1 + 0.2 * (n_B > 0)",
                ".9 + 0.4 * (n_A > 0)"
                ))
afd3 <- allFitnessEffects(genotFitness = gffd3,
                        frequencyDependentFitness = TRUE,
                        frequencyType = "abs")
userVars <- list(
    list(Name = "user_var_1",
        Value = 0
    ),
    list(Name = "user_var_2",
        Value = 0
    )
)

userVars <- createUserVars(userVars)

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

rules <- createRules(rules, afd3)

sfd3 <- oncoSimulIndiv( afd3,
                        model = "McFLD",
                        onlyCancer = FALSE,
                        finalTime = 100,
                        mu = 1e-4,
                        initSize = 5000,
                        sampleEvery = 0.001,
                        userVars = userVars,
                        rules = rules)


for(line in sfd3$other$userVarValues){
    if(line[3] == 0){
        testthat::expect_equal(line[1], 0)
        testthat::expect_equal(line[2], 0)
    }else if(line[3] <= 10){
        testthat::expect_equal(line[1], 0)
        testthat::expect_equal(line[2], 1)
    }else{
        testthat::expect_equal(line[1], 1)
        testthat::expect_equal(line[2], 2)
    }
}

cat(paste("\n Ending user var tests at", date()))