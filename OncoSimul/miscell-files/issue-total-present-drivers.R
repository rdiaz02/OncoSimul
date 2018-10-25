RNGkind("Mersenne-Twister")

set.seed(1)
lni <- 5 ## no fitness effects genes
fni <- 50 ## fitness effects genes
no <- 1e2
ft <- 2
s3 <- 2.0
mu <- 1e-5 ## easier to see
## noInt have no fitness effects, but can accumulate mutations
ni <- rep(0, lni)
## Those with fitness effects in one module, so
## neither fitness nor mut. rate blow up
gn <- paste(paste0("a", 1:fni), collapse = ", ")
f3 <- allFitnessEffects(epistasis = c("A" = s3),
                        geneToModule = c("A" = gn),
                        noIntGenes = ni
                        , drvNames = c("a2", "a5")
                        )
s3.ng <- oncoSimulIndiv(f3,
                        mu = mu,
                        mutationPropGrowth = FALSE,
                        finalTime =ft,
                        initSize = no,
                        onlyCancer = FALSE,
                        verbosity = 4,
                        seed = NULL)
s3.ng








