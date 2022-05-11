//     Copyright 2013-2021 Ramon Diaz-Uriarte

//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _USER_VAR_H__
#define _USER_VAR_H__


#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"
#include "new_restrict.h"
#include <cfloat>
#include <limits>
#include <Rcpp.h>
#include <iostream>
#include <random>
#include <set>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <exception>
#include "exprtk.h"
#include <cmath>

// we declare the needed symbols for the table
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

// We define what a rule is
typedef struct {
    std::string id; // identifier of the rule
    std::string condition; // condition that must be true for the rulle to execute
    std::string action; // action to be performed when the condition is satisfied
} Rule;

//Define a structure with all info asociated with user_vars
typedef struct{
    std::map<std::string, double> userVars; // maps a user variable name to its current value
    std::vector<Rule> rules; // rules that apply to the group of user variables
    std::map<std::string, double> mapGenoToPop; // variable that maps a genotype to its population
}UserVarsInfo;

// function that creates the UserVarsInfo structure
UserVarsInfo createUserVarsInfo(Rcpp::List Rules,
                                Rcpp::List userVars,
                                const fitnessEffectsAll& fitnessEffects, 
                                const std::vector<spParamsP>& popParams, 
                                std::vector<Genotype> Genotypes);

// function that creates a rule in memory
Rule createRule(std::string id, 
                std::string condition, 
                std::string action);

// Function that add an user_var to the map of user_vars
UserVarsInfo addUserVar(UserVarsInfo uvif, std::string name, double value);

// Function that add a rule to the array of rules
UserVarsInfo addRule(UserVarsInfo uvif, Rule r);

// function that destroys the rule in memory
UserVarsInfo destroyRule(UserVarsInfo uvif, Rule r);

// function that executes the whole list of rules for the user variables
void executeRules (UserVarsInfo& uvif,
                    double &currentTime,
                    std::map<std::string, double> birthMap,
                    std::map<std::string, double> deathMap,
                    std::map<std::string, double> mutationMap);

// function that compares two rules
int compareRules(Rule r1, Rule r2);

void printUserVarsInfo(UserVarsInfo uvif);
void printRule(Rule r);

#endif