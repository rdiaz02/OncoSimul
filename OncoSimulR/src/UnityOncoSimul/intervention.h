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

#ifndef _INTERVENTIONS_H__
#define _INTERVENTIONS_H__


#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"
#include "new_restrict.h"
#include "user_var.h"
#include "multivariate_hypergeometric.h"
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

#define NOT_PERIODICITY -1.0
#define NOT_REPS -0.5

// we declare the needed symbols for the table
typedef exprtk::symbol_table<double> symbol_table_t;
typedef exprtk::expression<double> expression_t;
typedef exprtk::parser<double> parser_t;

// We define what an intervention is
typedef struct {
    std::string id; // identifier of the intervention
    std::string trigger; // condition for the intervention to be executed
    std::string what_happens; // what the "action" will be in case the trigger "triggers"
    int repetitions; // how many repetitions to apply on the set
    float periodicity; // periodicity of the intervention (each 10 u.t.)
    float lastTimeExecuted; // last time from the current time it was executed
} Intervention;

//Define a structure with all info asociated with interventions
typedef struct{
    std::vector<Intervention> interventions; // array or list of interventions
    std::map<std::string, double> mapGenoToPop; // variable that maps a genotype to its population
}InterventionsInfo;

// function that creates the InterventionsInfo structure
InterventionsInfo createInterventionsInfo(Rcpp::List interventions, 
                                          const fitnessEffectsAll& fitnessEffects, 
                                          const std::vector<spParamsP>& popParams, 
                                          std::vector<Genotype> Genotypes);

// function that creates an intervention in memory
Intervention createIntervention(std::string id, 
                                std::string trigger, 
                                std::string what_happens, 
                                float periodicity, 
                                int repetitions);

// Function that add an intervention to the array of interventions
InterventionsInfo addIntervention(InterventionsInfo iif, Intervention i);

// function that destroys the intervention in memory
InterventionsInfo destroyIntervention(InterventionsInfo iif, Intervention i);

// function that executes the whole list of interventions
bool executeInterventions(InterventionsInfo& iif, 
                         UserVarsInfo& uvif,
                         double &totPopSize, 
                         double &currentTime, 
                         const fitnessEffectsAll& fitnessEffects, 
                         std::vector<Genotype> Genotypes, 
                         std::vector<spParamsP>& popParams);

// function that compares two interventions
int compareInterventions(Intervention i1, Intervention i2);

void printInterventionsInfo(InterventionsInfo iif);
void printIntervention(Intervention i);



//  These were the private functions
bool isValidEquation(std::string equation);
void parseWhatHappens(InterventionsInfo& iif, 
                     UserVarsInfo& uvif,
                     Intervention intervention,
                     double &totPopSize, double currentTime);
bool updatePopulations(InterventionsInfo& iif, 
                       const fitnessEffectsAll& fitnessEffects, 
                       const std::vector<Genotype>& Genotypes, 
                       std::vector<spParamsP>& popParams);
// functions for debugging
//void printIntervention(Intervention i);
//void printInterventionsInfo(InterventionsInfo iif);
// function that applies hypergeometric progressions to the reduction of the population 
void reduceTotalPopulation(InterventionsInfo& iif, double target);


#endif
