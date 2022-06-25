
## Copyright 2013-2021 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

# we load some functions that are needed to compile
# source("./R/new-restrict.R")

# this function create the interventions, verifies its correct  specification and returns those interventions
# so they can be processed correctly by C++
createInterventions <- function(interventions, genotFitness, frequencyType = "auto"){
    if((frequencyType == "abs") || (frequencyType == "auto")){
        return (adapt_interventions_to_cpp(verify_interventions(interventions), frequencyType, genotFitness))
    } else if(frequencyType == "rel"){
        stop("Frequency type as relative is not implemented... yet.")
    } else {
        stop("You have specified the freqType wrong. Review it please.")
    }
}

# this intervention transforms the genotype specification from the user-specified to the C++ one.
# transformIntervention will give more detail about how this works. This "transformation" is done for 
# trigger and what_happens attributes from the intervention.
adapt_interventions_to_cpp <- function(interventions, frequencyType, genotFitness) {

    for(i in 1:length(interventions)){
        interventions[[i]]$Trigger <- transform_intervention(as.character(interventions[[i]]$Trigger), genotFitness)
        interventions[[i]]$WhatHappens <- transform_intervention(as.character(interventions[[i]]$WhatHappens), genotFitness)
    }
    
    return(interventions)
}

# this function is the one in charge to transform the genotypes. The user will specify genotypes as A, B, "A,B,C"... etc.
# but in the C++ side, that processes operations, there is no "A", because when fitness is specified, there is an that changes n_A
# to n_1 or n_A_B_C to n_1_2_3. The function that does those "transformations", is all_orders_fv, we just reverse engineered this function
# and borrowed some funcionality so when the user specifies interventions that involve genotypes fitness's, the transformation can be made
transform_intervention <- function(sentence, genotInfo){
    
    prefix <- "n_"
    prefixre <- "^n_"

    # use functions that handle the change of identifiers, from A -> 1, B->2, A,B -> 1_2
    # functionality "borrowed" by all_orders_fv
    x <- genotInfo$fitnessLandscape[, -ncol(genotInfo$fitnessLandscape), drop = FALSE]
    pasted <- apply(x, 1, function(z) paste(sort(which(z == 1)), collapse = "_"))
    npasted <- apply(x, 1, function(z) paste(sort(colnames(x)[which(z == 1)]), collapse = "_"))
    flvars <- paste0(prefix, pasted)
    names(flvars) <- npasted
    flvars2 <- flvars
    
    # now we have mapped all the original genotypes names with the id'd ones 
    names(flvars2) <- paste0(prefix, names(flvars))
    # we obtain all the posible combinations of the genotypes
    rflvars2 <- rev(flvars2)
    full_rflvars <- all_orders_fv(rflvars2, prefix, prefixre)
    
    # to finish, we replace the previous names with the new ones
    return(stringr::str_replace_all(sentence,
                                    stringr::fixed(full_rflvars)))
}

## this function checks for inconsistencies in the specification of interventions
verify_interventions <- function(interventions){

    ## we check if there are interventions with the same ID.
    check_double_id(interventions)
    
    for(i in 1:length(interventions)){
        ## check that the interventions are lists
        print(paste0("Checking intervention: ", interventions[[i]]["ID"]))
        if(is.list(interventions[[i]]) == FALSE){
            stop("Type should be a list of lists")
        }
        
        ## check that exists ID, Trigger, WhatHappens and TimeSensitive in the list 
        if(!exists("ID", interventions[[i]])){
            stop("Attribute ID was not specified.")
        } else if(!exists("Trigger", interventions[[i]])){
            stop("Attribute Trigger was not specified.")
        } else if(!exists("WhatHappens", interventions[[i]])){
            stop("Attribute WhatHappens was not specified.")
        } else if(!exists("Repetitions", interventions[[i]]) || !exists("Periodicity", interventions[[i]]) ){
            stop("Either repetitions or Periodicity must be specified for the intervention")
        }
        
        ## check the data-type of the interventions[[i]]s
        if(!is.character(interventions[[i]]$Trigger) || !is.character(interventions[[i]]$WhatHappens) || !is.character(interventions[[i]]$ID)){
            stop("Interventions are wrongfully specified. Exit...")
        }

         ## case where user specifies negative periodicity of repetitions
        if (!exists("Periodicity", interventions[[i]]) && interventions[[i]]$Periodicity < 0){
            stop("The periodicity is negative in the intervention")
        }
        if(!exists("Repetitions", interventions[[i]]) && interventions[[i]]$Repetitions < 0){
            stop("The repetitions are negative in the intervention")
        }
        
        ## contemplate the case where only repetitions or only periodicity is specified.
        if(exists("Repetitions", interventions[[i]]) && !(exists("Periodicity", interventions[[i]])) ){
            interventions[[i]]$Periodicity = -1
        } else if(exists("Periodicity", interventions[[i]]) && !(exists("Repetitions", interventions[[i]])) ){ # If the user specifies the periodicity instead of the reps, then TimeSensitive should be setted to "Yes"
            interventions[[i]]$Repetitions = -0.5
        }

        ## in case user specifies periodicity or reps to infinity, we "cast" it to INT_MAX so C++ understands
        if (exists("Periodicity", interventions[[i]]) && interventions[[i]]$Periodicity == Inf){
            interventions[[i]]$Periodicity <- -1
        }
        if(exists("Repetitions", interventions[[i]]) && interventions[[i]]$Repetitions == Inf){
            interventions[[i]]$Repetitions <- .Machine$integer.max
        }

        check_what_happens(interventions[[i]]$WhatHappens)
    }

    return(interventions)
}

# this function checks that there are no interventions with the same ID.
check_double_id <- function(interventions){

    if(!is.list(interventions)){
        stop("This should be a list")
    }

    i <- 1
    buffer <- list()
    while(i <= length(interventions)){
        index_l <- i
        buffer[[i]] <- interventions[[i]]
        j <- 1
        while(j <= length(buffer)){
            if(j!= index_l){
                if(buffer[[j]]$ID == interventions[[i]]$ID){
                    stop("Check the interventions, there are 2 or more that have same IDs")
                }
            }
            j <- j + 1
        } 
        i <- i + 1
    }
}

# check that the what_happens is correctly specified.
check_what_happens <- function(what_happens){
    # what happens has this form:
    # <genotype_to_apply_some_operation> = <some_operation>
    # we need to assure that the left part is just before the "="

    string1 <- what_happens
    str_split = strsplit(string1, "\\s+")[[1]]
    # now we check that str_split[[2]] is "="

    if(str_split[[2]] != "="){
        stop("The specification of WhatHappens is wrong.\n It should be: 
        <genotype_to_apply_some_operation or total_population> = <some_operation>\n Exiting.")
    }

    flag = FALSE

    for(s in str_split){
        if(s == "="){
            if(flag == TRUE){
                stop("The specification of WhatHappens is wrong.\n It should be: 
        <genotype_to_apply_some_operation or total_population> = <some_operation>\n Exiting.")
            } else if(flag == FALSE){
                flag = TRUE
            }
        }
    }
}

