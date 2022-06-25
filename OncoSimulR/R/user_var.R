
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

# this function create the user_vars, varifies they are correct and returns those user_vars
# so they can be processed correctly by C++
createUserVars <- function(userVars){
    return (verify_user_vars(userVars))
}

# this function create the rules, verifies its correct  specification and returns those rules
# so they can be processed correctly by C++
createRules <- function(rules, genotFitness, frequencyType = "auto"){
    if((frequencyType == "abs") || (frequencyType == "auto")){
        return (adapt_rules_to_cpp(verify_rules(rules), frequencyType, genotFitness))
    } else if(frequencyType == "rel"){
        stop("Frequency type as relative is not implemented... yet.")
    } else {
        stop("You have specified the freqType wrong. Review it please.")
    }
}

# this transforms the genotype specification from the user-specified to the C++ one.
# transformRule will give more detail about how this works. This "transformation" is done for 
# condition and action attributes from the rule.
adapt_rules_to_cpp <- function(rules, frequencyType, genotFitness) {

    for(i in 1:length(rules)){
        rules[[i]]$Condition <- transform_rule(as.character(rules[[i]]$Condition), genotFitness)
        rules[[i]]$Action <- transform_rule(as.character(rules[[i]]$Action), genotFitness)
    }
    
    return(rules)
    
}

# this function is the one in charge to transform the genotypes. The user will specify genotypes as A, B, "A,B,C"... etc.
# but in the C++ side, that processes operations, there is no "A", because when fitness is specified, there is an that changes n_A
# to n_1 or n_A_B_C to n_1_2_3. The function that does those "transformations", is all_orders_fv, we just reverse engineered this function
# and borrowed some funcionality so when the user specifies rules that involve genotypes fitness's, the transformation can be made
transform_rule <- function(sentence, genotInfo){
    
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

## this function checks for inconsistencies in the specification of user variables
verify_user_vars <- function(userVars){

    ## we check if there are user variables with the same name.
    check_same_name(userVars)
    
    for(i in 1:length(userVars)){
        ## check that the user variables are lists
        print(paste0("Checking user variable: ", userVars[[i]]["Name"]))
        if(is.list(userVars[[i]]) == FALSE){
            stop("Type should be a list of lists")
        }
        
        ## check that exists Name, Value in the list 
        if(!exists("Name", userVars[[i]])){
            stop("Attribute Name was not specified.")
        } else if(!exists("Value", userVars[[i]])){
            stop("Attribute Value was not specified.")
        }
        
        ## check the data-type of the userVars[[i]]s
        if(!is.character(userVars[[i]]$Name) || !is.double(userVars[[i]]$Value)){
            stop("UserVars are wrongfully specified. Exit...")
        }   
    }

    return(userVars)
}


## this function checks for inconsistencies in the specification of rules
verify_rules <- function(rules){

    ## we check if there are rules with the same ID.
    check_double_rule_id(rules)
    
    for(i in 1:length(rules)){
        ## check that the rules are lists
        print(paste0("Checking rule: ", rules[[i]]["ID"]))
        if(is.list(rules[[i]]) == FALSE){
            stop("Type should be a list of lists")
        }
        
        ## check that exists ID, Condition, Action in the list 
        if(!exists("ID", rules[[i]])){
            stop("Attribute ID was not specified.")
        } else if(!exists("Condition", rules[[i]])){
            stop("Attribute Condition was not specified.")
        } else if(!exists("Action", rules[[i]])){
            stop("Attribute Action was not specified.")
        } 
        
        ## check the data-type of the rules[[i]]s
        if(!is.character(rules[[i]]$Condition) || !is.character(rules[[i]]$Action) || !is.character(rules[[i]]$ID)){
            stop("Rules are wrongfully specified. Exit...")
        }

        check_action(rules[[i]]$Action)
    }

    return(rules)
}


# this function checks that there are no user variables with the same name.
check_same_name <- function(userVars){

    if(!is.list(userVars)){
        stop("This should be a list")
    }

    i <- 1
    buffer <- list()
    while(i <= length(userVars)){
        index_l <- i
        buffer[[i]] <- userVars[[i]]
        j <- 1
        while(j <= length(buffer)){
            if(j!= index_l){
                if(buffer[[j]]$Name == userVars[[i]]$Name){
                    stop("Check the user variables, there are 2 or more that have same Names")
                }
            }
            j <- j + 1
        } 
        i <- i + 1
    }
}


# this function checks that there are no rules with the same ID.
check_double_rule_id <- function(rules){

    if(!is.list(rules)){
        stop("This should be a list")
    }

    i <- 1
    buffer <- list()
    while(i <= length(rules)){
        index_l <- i
        buffer[[i]] <- rules[[i]]
        j <- 1
        while(j <= length(buffer)){
            if(j!= index_l){
                if(buffer[[j]]$ID == rules[[i]]$ID){
                    stop("Check the rules, there are 2 or more that have same IDs")
                }
            }
            j <- j + 1
        } 
        i <- i + 1
    }
}

# check that the action is correctly specified.
check_action <- function(action){
    # TODO: finish
    # Actions have this form:
    # <variable to modify> = <some_operation or value>;<variable to modify> = <some_operation or value>;...
    # we need to assure that the left part is just before the "="

    string1 <- action
    str_split = strsplit(string1, ";")[[1]]
    # now we check that str_split[[2]] is "="

    for(act in str_split){
        string2 <- act
        str_split2 = strsplit(string1, "\\s+")[[1]]

        if(str_split2[[2]] != "="){
            stop("The specification of Action is wrong.\n It should be: 
            <variable to modify> = <some_operation or value>;<variable to modify> = <some_operation or value>;...\n Exiting.")
        }

        flag = FALSE

        for(s in str_split){
            if(s == "="){
                if(flag == TRUE){
                    stop("The specification of Action is wrong.\n It should be: 
                    <variable to modify> = <some_operation or value>;<variable to modify> = <some_operation or value>;...\n Exiting.")
                } else if(flag == FALSE){
                    flag = TRUE
                }
            }
        }
    }
    
}

