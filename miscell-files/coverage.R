## Run the coverage script


## library(shiny)
## library(DT)
## library(OncoSimulR)
## ## library(help = OncoSimulR)



library(covr)
setwd("../OncoSimulR")

## cov3 <- package_coverage(type = "all", combine_types = FALSE, quiet = FALSE)
## What matters most is tests after all. And
## running all and not combining (or one of them) screws up
## navigating to the code


## ## Make sure C++-11 code also covered.
## ## See https://github.com/jimhester/covr/issues/198#issuecomment-228568000
## options(covr.flags = c(CPPFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
##                        CXX1XFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
##                        CXXFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
##                        CFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
##                        LDFLAGS = "--coverage -fno-elide-constructors"))

## cov_O0 <- package_coverage(type = "tests", quiet = FALSE)




## But I run with -O3, functions are inlined, etc, so paths could
## differ. Do it also this way. And this deals correctly (almost, not
## fully) with the .h files and their structs.

options(covr.flags = c(CPPFLAGS = "-O3 -g --coverage",
                       CXX1XFLAGS = "-O3 -g --coverage",
                       CXXFLAGS = "-O3 -g --coverage",
                       CFLAGS = "-O3 -g --coverage",
                       LDFLAGS = "--coverage"))

system.time(cov_O3 <- package_coverage(type = "tests", quiet = FALSE))



save(file = "../miscell-files/coverage-results-O3.RData", cov_O3)


options(covr.flags = c(CPPFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
                       CXX1XFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
                       CXXFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
                       CFLAGS = "-O0 -g --coverage -fno-default-inline -fno-inline -fno-elide-constructors",
                       LDFLAGS = "--coverage -fno-elide-constructors"))

cov_O0 <- package_coverage(type = "tests", quiet = FALSE)

save(file = "../miscell-files/coverage-results-O0.RData", cov_O0)


## save(file = "../miscell-files/coverage-results.RData", cov_O0, cov_O3)


## cov3 <- package_coverage(type = "all", combine_types = FALSE, quiet = FALSE)
## cov2 <- package_coverage(type = "all", combine_types = TRUE, quiet = FALSE)

## cov4 <- package_coverage()
## save(file = "../miscell-files/coverage-results.RData", cov4, cov3, cov2)
## zero_coverage(cov4) 
## shine(cov_O0)

report(cov_O3)

