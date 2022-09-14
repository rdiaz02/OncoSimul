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

#include "intervention.h"



InterventionsInfo addIntervention(InterventionsInfo iif, Intervention i){
    //TODO: controlar que no exista ya una intervención con las mismas caracteristicas

    for(size_t k=0; k<iif.interventions.size(); k++){
        if(compareInterventions(iif.interventions[k], i) == 0){
            Rcpp::Rcout << "There are two interventions that are the same:";
            printIntervention(iif.interventions[k]);
            printIntervention(i);
            return iif;
        }
    }

    iif.interventions.push_back(i);
    
    return iif;
}

// function that creates an intervention in memory
Intervention createIntervention(std::string id, 
                                std::string trigger, 
                                std::string what_happens, 
                                float periodicity, 
                                int repetitions){
    Intervention i;
    i.id = id;
    i.trigger = trigger;
    i.what_happens = what_happens;
    i.periodicity = periodicity;
    i.repetitions = repetitions;
    i.lastTimeExecuted = 0.0;
    return i;
}

InterventionsInfo createInterventionsInfo(Rcpp::List interventions, 
                                          const fitnessEffectsAll& fitnessEffects, 
                                          const std::vector<spParamsP>& popParams, 
                                          std::vector<Genotype> Genotypes){
   
    // we declare the variables needed
    InterventionsInfo iif;
    Intervention iv;
    int totalEntries;
    //we use auxiliar variables to store the values from R
    std::vector<Rcpp::List> vectorOfList;
    std::vector<std::string> auxIDIntervention;
    std::vector<std::string> auxTriggerIntervention;
    std::vector<std::string> auxWhatHappensIntervention;
    std::vector<int> auxRepetitionsIntervention;
    std::vector<float> auxPeriodicityIntervention;

    for(int i=0; i<interventions.length(); i++){
        vectorOfList.push_back(Rcpp::as<Rcpp::List>(interventions[i]));
    }

    for(int i=0; i<interventions.length(); i++){
        Rcpp::List listItem = vectorOfList[i];

        auxIDIntervention.push_back(Rcpp::as<std::string>(listItem["ID"]));
        auxTriggerIntervention.push_back(Rcpp::as<std::string>(listItem["Trigger"]));
        auxWhatHappensIntervention.push_back(Rcpp::as<std::string>(listItem["WhatHappens"]));
        auxRepetitionsIntervention.push_back(Rcpp::as<int>(listItem["Repetitions"]));
        auxPeriodicityIntervention.push_back(Rcpp::as<float>(listItem["Periodicity"]));
    }

    totalEntries = interventions.length();

    //now we dump the info in the structs created
    for(int i=0; i<totalEntries; i++){
        iv = createIntervention(auxIDIntervention.at(i), auxTriggerIntervention.at(i), auxWhatHappensIntervention.at(i), auxPeriodicityIntervention.at(i), auxRepetitionsIntervention.at(i));
        iif = addIntervention(iif, iv);
    }

    // mapping for the genes and their population
    iif.mapGenoToPop = evalFVars(fitnessEffects, Genotypes, popParams, true);

    return iif;
}

int compareInterventions(Intervention i1, Intervention i2){

    // it is not allowed to have 2 interventions with the same ID.
    if(i1.id == i2.id){
        return 1;
    } 

    if (i1.trigger != i2.trigger) {
        return -1;
    }
    
    if (i1.what_happens != i2.what_happens){
        return -1;
    }
    
    if(i1.repetitions != i2.repetitions){
        return -1;
    } 
    
    if(i1.periodicity != i2.periodicity){
        return -1;
    }
    
    if(i1.lastTimeExecuted != i2.lastTimeExecuted){
        return -1;
    }
    // if they are equal in all aspects, then returns 0
    return 0;
}


InterventionsInfo destroyIntervention(InterventionsInfo iif, Intervention i){
    
    for (size_t z = 0; z < iif.interventions.size(); z++){
        if(compareInterventions(iif.interventions.at(z), i) == 0){
            iif.interventions.erase(iif.interventions.begin() + z);
            return iif;
        }
    }

    return iif;
}
