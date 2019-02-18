## Several checks of remaining things


## Binary objects that differ


library(OncoSimulR)


data(examplesFitnessEffects)
ex1 <- examplesFitnessEffects
rm(examplesFitnessEffects)

## now, load the object from the original version
load("~/caca/OncoSimul/OncoSimulR/data/examplesFitnessEffects.RData")

identical(names(ex1), names(examplesFitnessEffects))


for(i in 1:length(ex1)) {
    cat("\n doing i ", i, " \n")
    print(identical(ex1[[i]], examplesFitnessEffects[[i]]))
}

names(ex1[[1]])
names(examplesFitnessEffects[[1]])

## yes, with the new version there are more components in the objects.

for(i in 1:length(ex1)) {
    cat("\n doing i ", i, " \n")
    print(identical(ex1[[i]][c(1:7, 9:13)],
                    examplesFitnessEffects[[i]][c(1:7, 9:13)]))
}

## Comparing igraph objects is always a pain.
for(i in 1:length(ex1)) {
    print(identical(
        as_adj(ex1[[i]][[8]], edges = TRUE, sparse = FALSE),
        as_adj(examplesFitnessEffects[[i]][[8]], edges = TRUE, sparse = FALSE)
    ))
}

##################################

rm(list = ls())

load("~/Proyectos/OncoSimul/OncoSimulR/inst/testdata_fee.RData")
fee0 <- fee
feex0 <- feex
rm(fee)
rm(feex)

load("~/caca/OncoSimul/OncoSimulR/inst/testdata_fee.RData")

identical(fee[c(1:16)],
          fee0[c(1:16)])

identical(feex,
          feex0)





##################################

rm(list = ls())

load("~/Proyectos/OncoSimul/miscell-files/fee.RData")
fee0 <- fee
rm(fee)

load("~/caca/OncoSimul/miscell-files/fee.RData")

identical(fee[c(1:16)],
          fee0[c(1:16)])

## the "x" object from the original fee.RData is no longer present.
## Can we use the former object?

## But test-fix.R, under miscell-files, is an old file that does not even
## tests anything as such. All of that is in test.Z-fixation.R now, using
## testdata_fee.RData



