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
void reduceTotalPopulation(InterventionsInfo& iif, 
                           double target, 
                           double totPopSize);

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

bool executeInterventions(InterventionsInfo& iif,
                         UserVarsInfo& uvif, 
                         double &totPopSize, 
                         double &currentTime, 
                         const fitnessEffectsAll& fitnessEffects, 
                         std::vector<Genotype> Genotypes, 
                         std::vector<spParamsP>& popParams){

    // Now we add all the info needed for the symbol table so exprtk can operate 
    symbol_table_t symbol_table;

    for(auto& iterator : uvif.userVars) {
        symbol_table.add_variable(iterator.first, iterator.second);
    }

    for(auto& iterator : iif.mapGenoToPop) {
        symbol_table.add_variable(iterator.first, iterator.second);
    } 

    double N = totPopSize;

    if(N == 0){
        throw std::invalid_argument("Total Population = 0. There is nothing to intervene.\n");
    }

    double T = currentTime;

    symbol_table.add_constant("N", N);//We reserve N to total population size
    symbol_table.add_constant("T", T); //Pass current time to exprtk
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    bool interventionDone = false;

    //we iterate over the user-specified interventions 
    for(auto& intervention: iif.interventions){
        parser_t parser_trigger;
        //if parser fails to compile, throws exception
        if (!parser_trigger.compile(intervention.trigger, expression)){
            // error control, just in case the parsing it's not correct
            Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

            for (std::size_t i = 0; i < parser_trigger.error_count(); ++i){
                typedef exprtk::parser_error::type error_t;
                error_t error = parser_trigger.get_error(i);
                // FIXMEmaybe: Use warning or error to capture it easily in tests?
                REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                    i,
                    error.token.position,
                    exprtk::parser_error::to_str(error.mode).c_str(),
                    error.diagnostic.c_str(),
                    intervention.trigger.c_str());
            }
            std::string errorMessage = "The expression was imposible to parse.";
            throw std::invalid_argument(errorMessage);
        } else {
            //a trigger is just a TRUE/FALSE condition
            if(expression.value() == 1){
                parser_t parser_wh;
                if(intervention.repetitions >= 0 && intervention.periodicity == NOT_PERIODICITY){ // case where interventions are based only in repetitions
                    //if parser fails to compile, throws exception
                    if (!parser_wh.compile(intervention.what_happens, expression)){
                        // error control, just in case the parsing it's not correct
                        Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

                        for (std::size_t i = 0; i < parser_wh.error_count(); ++i){
                            typedef exprtk::parser_error::type error_t;
                            error_t error = parser_wh.get_error(i);
                            // FIXMEmaybe: Use warning or error to capture it easily in tests?
                            REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                                i,
                                error.token.position,
                                exprtk::parser_error::to_str(error.mode).c_str(),
                                error.diagnostic.c_str(),
                                intervention.what_happens.c_str());
                        }
                        std::string errorMessage = "The expression was imposible to parse.";
                        throw std::invalid_argument(errorMessage);
                    } else 
                        parseWhatHappens(iif, uvif, intervention, N, T);
                    // we reduce by one the number of interventions
                    intervention.repetitions--;
                    // we update the last time it was executed (debugging purposes)
                    intervention.lastTimeExecuted = T;
                    // we update interventionDone flag
                    interventionDone = true;

                } else if(intervention.repetitions >= 0 && intervention.periodicity > 0) { // case there is periodicity but also repetitions
                    if((T - intervention.lastTimeExecuted) >= intervention.periodicity){ // with condition satisfied we execute the intervention

                        if (!parser_wh.compile(intervention.what_happens, expression)){
                            Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

                            for (std::size_t i = 0; i < parser_wh.error_count(); ++i){
                                typedef exprtk::parser_error::type error_t;
                                error_t error = parser_wh.get_error(i);
                                // FIXMEmaybe: Use warning or error to capture it easily in tests?
                                REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                                    i,
                                    error.token.position,
                                    exprtk::parser_error::to_str(error.mode).c_str(),
                                    error.diagnostic.c_str(),
                                    intervention.what_happens.c_str());
                            }
                            std::string errorMessage = "The expression was imposible to parse.";
                            throw std::invalid_argument(errorMessage);
                        } else 
                            parseWhatHappens(iif, uvif, intervention, N, T);
                        // update new lastTimeExecuted
                        intervention.lastTimeExecuted = T;
                        //  update amount of repetitions
                        intervention.repetitions--;
                        // we update interventionDone flag
                        interventionDone = true;
                    } 
                } else if (intervention.periodicity > 0 && intervention.repetitions == NOT_REPS) { // case where only periodicty is specified
                    if((T - intervention.lastTimeExecuted) >= intervention.periodicity){
                        if (!parser_wh.compile(intervention.what_happens, expression)){
                            Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

                            for (std::size_t i = 0; i < parser_wh.error_count(); ++i){
                                typedef exprtk::parser_error::type error_t;
                                error_t error = parser_wh.get_error(i);
                                // FIXMEmaybe: Use warning or error to capture it easily in tests?
                                REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                                    i,
                                    error.token.position,
                                    exprtk::parser_error::to_str(error.mode).c_str(),
                                    error.diagnostic.c_str(),
                                    intervention.what_happens.c_str());
                            }
                            std::string errorMessage = "The expression was imposible to parse.";
                            throw std::invalid_argument(errorMessage);
                        } else 
                            parseWhatHappens(iif, uvif, intervention, N, T);
                        // update new lastTimeExecuted
                        intervention.lastTimeExecuted = T;
                        // we update interventionDone flag
                        interventionDone = true;
                    }
                }
            }
        }
    }

    // Now with the structure storing all the changes, we need to store the data in their own
    // original structures where the data was sourced from 
    // once the structure is updated, we update the structures that store the info while the simulation is running
    if(interventionDone){
        if(!updatePopulations(iif, fitnessEffects, Genotypes, popParams)){
                throw std::runtime_error("There was an issue updating the populations while intervening.\n");
        }
        return true;
    }
    return false;   
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////// PRIVATE FUNCTIONS ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void parseWhatHappens(InterventionsInfo& iif,
                     UserVarsInfo& uvif,
                     Intervention intervention, 
                     double &totPopSize, 
                     double currentTime){
    
    // now we need to parse the "what_happens" intervention
    // TODO: raise exception, malformed what_happens specification
    if(!isValidEquation(intervention.what_happens)){
        throw std::runtime_error("The intervention was wrongfully specified.\n");
    }

    bool totalPopFlag = false;

    symbol_table_t symbol_table;

    for(auto& iterator : uvif.userVars) {
        symbol_table.add_variable(iterator.first, iterator.second);
    }
    
    for(auto& iterator : iif.mapGenoToPop) {
        symbol_table.add_variable(iterator.first, iterator.second);
    } 

    double N = totPopSize;
    double T = currentTime;

    symbol_table.add_constant("N", N);//We reserve N to total population size
    symbol_table.add_constant("T", T); //Pass current time to exprtk
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);
    
    //Variables needed for tokeninzing
    std::vector <std::string> tokens;
    std::stringstream check1(intervention.what_happens); 
    std::string intermediate; 

    // Tokenizing whathappens '=' 
    while(getline(check1, intermediate, '=')){ 
        tokens.push_back(intermediate); 
    } 
    // now we see if the operation affects total population or genotype popolation
    std::string leftMostWhatHappens = tokens[0];
    leftMostWhatHappens.erase(std::remove(leftMostWhatHappens.begin(), leftMostWhatHappens.end(), ' '), leftMostWhatHappens.end());
    if(leftMostWhatHappens.compare("N") == 0)
        totalPopFlag = true;

    // we use the right-most side of the equation
    std::string rightMostWhatHappens = tokens[1];
    parser_t parser_rwh;
    if (!parser_rwh.compile(rightMostWhatHappens, expression)){
        // error control, just in case the parsing it's not correct
        Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

        for (std::size_t i = 0; i < parser_rwh.error_count(); ++i){
            typedef exprtk::parser_error::type error_t;
            error_t error = parser_rwh.get_error(i);
            // FIXMEmaybe: Use warning or error to capture it easily in tests?
            REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                i,
                error.token.position,
                exprtk::parser_error::to_str(error.mode).c_str(),
                error.diagnostic.c_str(),
                rightMostWhatHappens.c_str());
        }
        std::string errorMessage = "The expression was imposible to parse.";
        throw std::invalid_argument(errorMessage);
    } else{
        // value cant have decimals (can't exist a half-cell)
        double res = floor(expression.value());

        // once the value is calculated, we must assure if the operation is for the total population
        // or for some specific-genotype
        if (totalPopFlag && (res > N)) {
            // TODO: Throw exception of some kind, this CANNOT happen by any means
            throw std::runtime_error("You have specified an intervention that is not allowed.");
        } 

        // check that user does not create a WhatHappens that creates population: n_A = n_A * 2
        if(!totalPopFlag){
            try { //check with the left part of the equation finding the value for the n_*
                const double& value = iif.mapGenoToPop.at(leftMostWhatHappens);
                if(res > value){
                    Rcpp::Rcerr << "In intervention:" << intervention.id << " with WhatHappens: " << intervention.what_happens << ". You cannot intervene to generate more population.";
                    throw std::runtime_error("You have specified an intervention that is not allowed.");
                }
            }
            catch (const std::out_of_range&) {
                Rcpp::Rcout << "Key \"" << leftMostWhatHappens.c_str() << "\" not found" << std::endl;
            }
        }

        if(totalPopFlag && res == N){ // this case is absurd, but it might happen, we just return.
            return;
        } else if(totalPopFlag && (res < N)){// reduce total amount of population using hipergeometric distribution
            reduceTotalPopulation(iif, res, totPopSize);
        } else { // update new value for genotype
            std::map<std::string, double>::iterator it = iif.mapGenoToPop.find(leftMostWhatHappens); 
            if(it != iif.mapGenoToPop.end()){
                it->second = res; 
            }
        }
    }
}

// This function is needed if what we are trying to descrease is the whole population, and not just the population of 1 genotype
// nn by default is equal 1
// n array con las poblaciones de cada genotipo y su ncols(diferentes tipos de genotipos en la población)
// target is the target size, to which number the population would get reduced to.
void reduceTotalPopulation(InterventionsInfo& iif, double target, double totPopSize){
    
    // first we take all the population from the structure, and we create a vector with its populations
    std::vector<double> populations;
    Rcpp::NumericMatrix rcpp_mhgeo_distribution; 
    Rcpp::NumericMatrix rcpp_populations_matrix;
    Rcpp::NumericVector rcpp_populations;
    Rcpp::NumericVector rcpp_target;
    //auxiliar variable to check the structure (and the atribute mapGenoToPop) is well populated
    double totalPop = 0.0;
    for(auto map : iif.mapGenoToPop){
        totalPop += map.second;
        populations.push_back(map.second);
    }
    // we convert the vector to something R can understand
    rcpp_populations = Rcpp::wrap(populations);
    rcpp_populations.attr("dim") = Dimension(1, populations.size());

    rcpp_populations_matrix = Rcpp::wrap(rcpp_populations);
    rcpp_target = Rcpp::wrap(target);
    
    // then, we specify the total genotypes of the populations and we obtain a distribution
    rcpp_mhgeo_distribution = my_rmvhyper(1, rcpp_populations_matrix, rcpp_target);

    populations = Rcpp::as<std::vector<double>>(rcpp_mhgeo_distribution);

    //quick check before creating the matrix CANT BE DONE HERE (totPopSize is not updated until crude is filled)
    //if(totPopSize != totalPop){
    //    throw std::runtime_error("TotalPop != totPopSize, exiting...");
    //}

    int i=0;
    for(auto &map : iif.mapGenoToPop){
        //check if it goes out of bounds
        map.second = populations[i];
        i++;
    }
    
}

// this function updates the population by genotypes
bool updatePopulations(InterventionsInfo& iif, 
                       const fitnessEffectsAll& fitnessEffects, 
                       const std::vector<Genotype>& Genotypes, 
                       std::vector<spParamsP>& popParams){

    
    std::map<std::string, std::string> fvarsmap = fitnessEffects.fitnessLandscape.flfVarsBmap;
    if(fvarsmap.empty()){
        fvarsmap = fitnessEffects.fitnessLandscape.flfVarsDmap;
    }
    if(fvarsmap.empty()){
        return false;
    }
    std::string freqType = fitnessEffects.frequencyType;

    for(auto map : iif.mapGenoToPop){
        // find the population associated to a genotype
        for(const auto& iterator : fvarsmap){
            if(map.first == iterator.second){
                std::vector<int> genotype = stringVectorToIntVector(iterator.first);//genotype (as int vector)
                int position = findPositionInGenotypes(Genotypes, genotype);
                //just to make sure...
                if(position != 0){
                    int realPos = position - 1;
                    if(freqType == "abs"){
                        popParams[realPos].popSize = map.second;
                        // maybe add "rel" in the future?
                    } else {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

void printIntervention(Intervention i){

    Rcpp::Rcout << i.id << " info:\n";
    Rcpp::Rcout << "\t Trigger: " << i.trigger << "\n";
    Rcpp::Rcout << "\t What Happens: " << i.what_happens << "\n";
    Rcpp::Rcout << "\t Repetitions: " << i.repetitions << "\n";
    Rcpp::Rcout << "\t Periodicity: " << i.periodicity << "\n";
    Rcpp::Rcout << "\t Last Time Executed: " << i.lastTimeExecuted << "\n";
}

void printInterventionsInfo(InterventionsInfo iif){

    for(auto intervention: iif.interventions){
        printIntervention(intervention);
        // print the info associated with genotypes and their population 
    }

    for(auto map : iif.mapGenoToPop) {
            Rcpp::Rcout << "Genotype: " << map.first << " Population: " << map.second << "\n";
    }

}

// private function that checks that the equations specified by the user are correctly specified
// if no "=" or two "0" are found, returns false, if just one "=" is found, the returns true.
bool isValidEquation(std::string equation){
    bool flag = false;
    for(std::string::size_type i=0; i<equation.size(); i++){
        if (flag == true){
            if(equation[i] == '='){
                return false;
            }
        }
        if(equation[i] == '='){
            flag = true;
        }
    }

    return flag;
}



