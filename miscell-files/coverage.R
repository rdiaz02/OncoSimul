## Run the coverage script

library(covr)
library(shiny)
library(DT)
## library(OncoSimulR)
## library(help = OncoSimulR)

setwd("../OncoSimulR")

## cov3 <- package_coverage(type = "all", combine_types = FALSE, quiet = FALSE)
## What matters most is tests after all. And
## running all and not combining (or one of them) screws up
## navigating to the code
cov4 <- package_coverage(type = "tests", combine_types = TRUE, quiet = FALSE)
## cov4 <- package_coverage()
save(file = "../miscell-files/coverage-results.RData", cov4)
## zero_coverage(cov4) 
shine(cov4)
