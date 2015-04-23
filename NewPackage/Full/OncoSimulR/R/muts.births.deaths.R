## this can help checking if we are OK assuming that mutation is much
## smaller than birth and rate, as assumed by BNB

muts.births.deaths <- function(mutation, selection, drivers = 30) {

    

    maxMut <- drivers * max(mutation)
    mutatedDrivers <- 0:drivers
    leftToMutate <- drivers - mutatedDrivers

    ## Bozic: in a given clone, death and total mutation a ratio of number
    ## of drivers left to mutate or already mutated
    ## minDeathBozic <- 1 ## not needed
    ## maxBirthBozic <- minBirthBozic <- 1 ## not needed
    totalMutBozic <- max(mutation) * leftToMutate
    deathBozic <- (1 - max(selection))^mutatedDrivers
    ratioDeathMutationBozic <- deathBozic/totalMutBozic
    worstRatioDeathMutationBozic <- min(ratioDeathMutationBozic)
    worstRatioBirthMutationBozic <- 1/max(totalMutBozic)

    ## check:
    mdb <-
        (1/max(mutation)) *
            min( ((1-max(selection))^(0:drivers))/(drivers:0) ) 
    stopifnot(all.equal(mdb, worstRatioDeathMutationBozic))
    

    ## McFl
    minDeathMcFl <- 1
    totalMutMcFl <- totalMutBozic

       
    ## next is when all drivers hit but not satisfied and we set the max
    ## s in those cases to selection for drivers, at most.
    ## In fact is often much smaller.
    birthMcFl <- 1/((1 + max(selection))^mutatedDrivers)
    ratioBirthMutationMcFl <- birthMcFl/totalMutMcFl
    worstRatioBirthMutationMcFl <- min(ratioBirthMutationMcFl)
    worstRatioDeathMutationMcFl <- minDeathMcFl/max(totalMutMcFl)

    
    ## Exp
    ## For Exp never an issue, as long as total mutation rate is less than
    ## 1: Mutation rate is birth * mutation

    minDeathExp <- 1
    maxBirthExp <- (1 + max(selection))^mutatedDrivers
    minBirthExp <- (1 - max(selection))^mutatedDrivers
    mutExpMax <- max(mutation) * leftToMutate * maxBirthExp
    mutExpMin <- max(mutation) * leftToMutate * minBirthExp

    maxMutExp <- max(mutExpMax, mutExpMin)
    worstRatioDeathMutationExp <- 1/maxMutExp
    ##next are of course identical
    worstRatioBirthMutationExpMax <- min(maxBirthExp/mutExpMax)
    worstRatioBirthMutationExpMin <- min(minBirthExp/mutExpMin)



    ## check: both have 1 in the numerator, if McFl with defaults
    stopifnot(all.equal(worstRatioDeathMutationMcFl,
                        worstRatioBirthMutationBozic))


    ## check:
    stopifnot(all.equal(worstRatioBirthMutationMcFl,
                        worstRatioDeathMutationExp))

    ## check: both are just 1/(mutation * drivers)
    stopifnot(all.equal(worstRatioDeathMutationMcFl,
                        worstRatioBirthMutationExpMax))
    stopifnot(all.equal(worstRatioBirthMutationBozic,
                        worstRatioBirthMutationExpMax))
    stopifnot(all.equal(worstRatioBirthMutationBozic,
                        worstRatioDeathMutationMcFl))
    stopifnot(all.equal(worstRatioDeathMutationMcFl,
                        1/(drivers * max(mutation) ) ))

    ## message(paste("\n The value for ratio of death to mutation in ",
    ##               " McFl are approximate; they can be slightly smaller ",
    ##               " because of fluctuations in population size.",
    ##               " I have arbitrarily made them 10% smaller.",
    ##               " The above calculations can be badly off if you ",
    ##               " modify the default K."))

    ## adjust
    
    worstRatioDeathMutationMcFl <- (0.80 * minDeathMcFl)/max(totalMutMcFl)
    
    df <- data.frame(Model = c("Exp", "Bozic", "McFL"),
                     worstRatioBirthMutation = c(worstRatioBirthMutationExpMax,
                         worstRatioBirthMutationBozic,
                         worstRatioBirthMutationMcFl),
                     worstRatioDeathMutation = c(worstRatioDeathMutationExp,
                         worstRatioDeathMutationBozic,
                            worstRatioDeathMutationMcFl)
                     )
    

    ## Former, naive way

    ## minMut <- drivers * min(mutation)
    ## maxMutExp <- maxBirthExp * maxMut
    ## print(maxBirthExp)

    
    ## minBirthExp <- (1 - max(selection))^drivers
    ## minDeathExp <- 1
    ## maxBirthExp <- (1 + max(selection))^drivers

    ## minBirthMcFl <- 1/((1 + max(selection))^drivers)
    ## minDeathMcFl <- 1
    ## maxBirthMcFl <- (1 + max(selection))^drivers
    
    ## maxMut <- drivers * max(mutation)
    ## ## print(maxMut)

    
    ## df <- data.frame(Model = c("Exp", "Bozic", "McFL"),
    ##                  MinBirth = c(minBirthExp, minBirthBozic, minBirthMcFl),
    ##                  MinDeath = c(minDeathExp, minDeathBozic, minDeathMcFl),
    ##                  MaxBirth = c(maxBirthExp, 1, maxBirthMcFl),
    ##                  MaxMutation = c(maxMutExp, maxMut, maxMut)
    ##                  )

    ## df$RatioBirthMutation <- df$MinBirth/df$MaxMutation
    ## df$RatioBirthMutation[1] <- 1/maxMut ## always this
    ## df$RatioDeathMutation <- df$MinDeath/df$MaxMutation
    ## ## cat(paste("\n Exp: min birth ", minBirthExp,
    ## ##           "  min death ", minDeathExp, " max mut ", maxMut))

    ## ## cat(paste("\n Bozic: min birth ", minBirthBozic,
    ## ##           "  min death ", minDeathBozic, " max mut ", maxMut))

    ## ## cat(paste("\n McFL: min birth ", minBirthMcFl,
    ## ##           "  min death ", minDeathMcFl, " max mut ", maxMut))
    return(df)
}
