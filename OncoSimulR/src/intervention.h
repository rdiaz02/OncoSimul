#ifndef _INTERVENTIONS_H__
#define _INTERVENTIONS_H__


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
#include "exprtk.h"

// We define what an intervention might be like
typedef struct {
    std::string id; // identifier of the intervention
    std::string trigger; // condition for the intervention to be executed
    std::string what_happens; // what the "action" will be in case the trigger "triggers"
    int repetitions; // how many repetitions to apply on the set
    float periodicity; // periodicity of the intervention (each 10 u.t.)
    float lastTimeExecuted; // last time from the current time it was executed
    std::string flagTimeSensitiveIntervention; // flag that determines if the intervention dependes on some periodicity "Yes(Y)" o "No(N)"
} Intervention;

//Define a structure with all info asociated with interventions
typedef struct{
    std::vector<Intervention> interventions; // array or list of interventions
    std::map<std::string, double> mapGenoToPop; // variable that maps a genotype to its population
}InterventionsInfo;

// function that creates the InterventionsInfo structure
InterventionsInfo createInterventionsInfo(Rcpp::List interventions, fitnessEffectsAll& fitnessEffects, const std::vector<spParamsP>& popParams, const Genotype ge);

// function that creates an intervention in memory
Intervention createIntervention(std::string id, std::string trigger, std::string what_happens, float duration, int repetitions);

// Function that add an intervention to the array of interventions
InterventionsInfo addIntervention(InterventionsInfo iif, Intervention i);

// function that destroys the intervention in memory
InterventionsInfo destroyIntervention(InterventionsInfo iif, Intervention i);

// function that executes the whole list of interventions
bool executeInterventions(Rcpp::List interventions, int &totPopSize, double &currentTime, fitnessEffectsAll& fitnessEffects, const Genotype& ge, const std::vector<Genotype>& Genotypes);

// function that applies hypergeometric progressions to the reduction of the population 
void reducePopulation(InterventionsInfo * iif, double target);

// function that compares two interventions
int compareInterventions(Intervention i1, Intervention i2);

// function that tracks and updates the simulation time depending on the intervention (and its trigger)
// for example if we have a time sensitive trigger like (each 20 seconds, reduce the poplation by 20%) this function
// will track the current time and the previous time the intervention was executed to then, update the time.
// void updateTimeOfSimulation(double &currentTime, Intervention * currentIntervention);

// function that maps every genotype with their own population
//std::map<std::string, int> mapGenesToPop(const std::vector<spParamsP>& popParams);

//bool parseWhatHappens(InterventionsInfo * iif, Intervention intervention, double totPopSize, double currentTime);

bool updatePopulations(InterventionsInfo * iif /*, fitnessEffectsAll& fitnessEffects, const std::vector<spParamsP>& popParams*/);

#endif