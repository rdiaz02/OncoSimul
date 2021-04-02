#include "intervention.h"

bool isValidEquation(std::string equation);
bool parseWhatHappens(InterventionsInfo *iif, Intervention intervention, double totPopSize, double currentTime);

InterventionsInfo addIntervention(InterventionsInfo iif, Intervention i){
    //TODO: controlar que no exista ya una intervención con las mismas caracteristicas
    iif.interventions.push_back(i);
    
    return iif;
}

Intervention createIntervention(std::string id, std::string trigger, std::string what_happens, float periodicity, int repetitions, std::string flagTimeSensitiveIntervention){
    Intervention i;
    i.id = id;
    i.trigger = trigger;
    i.what_happens = what_happens;
    i.periodicity = periodicity;
    i.repetitions = repetitions;
    i.flagTimeSensitiveIntervention = flagTimeSensitiveIntervention;
    return i;
}

InterventionsInfo createInterventionsInfo(Rcpp::List interventions, fitnessEffectsAll& fitnessEffects, const std::vector<spParamsP>& popParams, const std::vector<Genotype>& Genotypes){
   
    // we declare the variables needed
    InterventionsInfo iif;
    Intervention iv;
    int totalEntries;
    //we use auxiliar variables to store the values from R
    std::vector<std::string> auxIDIntervention = Rcpp::as<std::vector<std::string> >(interventions["ID"]);
    std::vector<std::string> auxTriggerIntervention = Rcpp::as<std::vector<std::string> >(interventions["Trigger"]);
    std::vector<std::string> auxWhatHappensIntervention = Rcpp::as<std::vector<std::string> >(interventions["WhatHappens"]);
    std::vector<int> auxRepetitionsIntervention = Rcpp::as<std::vector<int> >(interventions["Repetitions"]);
    std::vector<float> auxPeriodicty = Rcpp::as<std::vector<float> >(interventions["Periodicity"]);
    std::vector<std::string> auxFlagTimeSensitiveIntervention = Rcpp::as<std::vector<std::string> >(interventions["TimeSensitive"]);
    totalEntries = auxIDIntervention.size();

    //TODO: necesitamos chequear si la intervención tiene algún tipo de periocidad o no

    //now we dump the info in the structs created
    for(int i=0; i<totalEntries; i++){
        iv = createIntervention(auxIDIntervention.at(i), auxTriggerIntervention.at(i), auxWhatHappensIntervention.at(i), auxPeriodicty.at(i), auxRepetitionsIntervention.at(i), auxFlagTimeSensitiveIntervention.at(i));
        iif = addIntervention(iif, iv);
    }

    // mapping for the genes and their population
    iif.mapGenoToPop = evalFVars(fitnessEffects, Genotypes, popParams);

    return iif;
}

int compareInterventions(Intervention i1, Intervention i2){

    // comparamos todos los casos, si falla en alguno retorna negativo
    if(i1.id != i2.id){
        return -1;
    } else if (i1.trigger != i2.trigger) {
        return -1;
    } else if (i1.what_happens != i2.what_happens){
        return -1;
    }  else if(i1.repetitions != i2.repetitions){
        return -1;
    } else if(i1.periodicity != i2.periodicity){
        return -1;
    } else if(i1.lastTimeExecuted != i2.lastTimeExecuted){
        return -1;
    } else if (i1.flagTimeSensitiveIntervention == i2.flagTimeSensitiveIntervention){
        return -1;
    }
    // if they are equal in all aspects, then returns 0
    return 0;
}

InterventionsInfo destroyIntervention(InterventionsInfo iif, Intervention i){
    
    for (int z=0; z<iif.interventions.size(); z++){
        if(compareInterventions(iif.interventions.at(z), i) == 0){
            iif.interventions.erase(iif.interventions.begin() + z);
            return iif;
        }
    }

    printf("The intervention you are trying to destroy hasn't been found.\n");
    return iif;
}

//TODO: aclarar los argumentos que se van a pasar a la función, creo que faltan...
bool execute_interventions(Rcpp::List interventions, int &totPopSize, double &currentTime, fitnessEffectsAll& fitnessEffects, const Genotype& ge, const std::vector<Genotype>& Genotypes, const std::vector<spParamsP>& popParams){
    //create the structure with all the information of the interventions
    InterventionsInfo iif = createInterventionsInfo(interventions, fitnessEffects, popParams, Genotypes);

    // TODO: revisar si se refiera a los genotipos por letras o por números
    // Now we add all the info needed for the symbol table so exprtk can operate 
     // se declaran los objetos necesarios para hacer la tabla de símbolos
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;

    symbol_table_t symbol_table;
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

    //iteramos sobre las intervenciones
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
            parser_t parser_wh;
            if(expression.value()){
                if(intervention.repetitions > 0){ // caso de que sean intervenciones basadas en repeticiones
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
                    } else if(!parseWhatHappens(&iif, intervention, N, T)){
                        printf("Something went wrong.\n");
                        return false;
                    }
                    // we reduce by one the number of interventions
                    intervention.lastTimeExecuted = T;
                    intervention.repetitions--;

                } else if(intervention.flagTimeSensitiveIntervention == "Yes" || intervention.flagTimeSensitiveIntervention == "Y") { // case there are time-based interventions (each 5 seconds, do this)
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
                        } else if(!parseWhatHappens(&iif, intervention, N, T)){
                            printf("Something went wrong.\n");
                            return false;
                        }
                        // update new lastTimeExecuted
                        intervention.lastTimeExecuted = T;
                    } 
                } else { // case where just an intervention needs to be executed just one time
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
                    } else if(!parseWhatHappens(&iif, intervention, N, T)){
                        printf("Something went wrong.\n");
                        return false;
                    }
                }
            }
        }
    }

    // TODO: Now with the structure storing all the changes, we need to store the data in their own
    // original structures where the data was sourced from (comes from Fitness or spParams?).
    

    return true;
}

bool parseWhatHappens(InterventionsInfo * iif, Intervention intervention, double totPopSize, double currentTime){
    
    // now we need to parse the "what_happens" intervention
    // TODO: raise exception, malformed what_happens specification
    if(!isValidEquation(intervention.what_happens)){
        // por el momento printf y return
        printf("The intervention was wrongfully specified.\n");
        return false;
    }

     // se declaran los objetos necesarios para hacer la tabla de símbolos
    typedef exprtk::symbol_table<double> symbol_table_t;
    typedef exprtk::expression<double> expression_t;
    typedef exprtk::parser<double> parser_t;
    bool totalPopFlag = false;

    symbol_table_t symbol_table;
    for(auto& iterator : iif->mapGenoToPop) {
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
        // once the value is calculated, we must assure if the operation is for the total population
        // or for some specific-genotype
        double res = expression.value();
        if(totalPopFlag){// reduce total amount of population using hipergeometric distribution
            reducePopulation(iif, res);
        } else { // update new value for genotype
            std::map<std::string, double>::iterator it = iif->mapGenoToPop.find(leftMostWhatHappens); 
            if(it != iif->mapGenoToPop.end()){
                it->second = res; 
            }
        }
    }

    return true;
}


// This function is needed if what we are trying to descrease is the whole population, and not just the population of 1 genotype
// Necesitamos N como numero total de la población
// k pertenecen a la categoría A
// N-k pertenecen a la B
// medir la probabilidad de obtener x (0<=x<=k) elementos de una categoría A
// target es el tamaño objetivo, el n, la muestra que se saca de la población total
void reducePopulation(InterventionsInfo * iif, double target){
    //ERROR CONTROL
    if(iif == NULL){
        printf("La estructura InterventionsInfo está vacía en reducePopulation.\n");
        return;
    }

    
}

// funcion que actualiza las poblaciones una vez que las intervenciones han sido ejecutadas
bool updatePopulations(InterventionsInfo * iif /*, fitnessEffectsAll& fitnessEffects, const std::vector<spParamsP>& popParams*/){
    return true;
}

// functions for debugging
void printIntervention(Intervention i){

    std::cout << "Intervention " << i.id << " info:\n";
    std::cout << "\t Trigger: " << i.trigger << "\n";
    std::cout << "\t What Happens: " << i.what_happens << "\n";
    std::cout << "\t Repetitions: " << i.repetitions << "\n";
    std::cout << "\t Periodicity: " << i.periodicity << "\n";
    std::cout << "\t Last Time Executed: " << i.lastTimeExecuted << "\n";
    std::cout << "\t Is time sensitive? " << i.flagTimeSensitiveIntervention << "\n";
}

void printInterventionsInfo(InterventionsInfo iif){

    for(auto intervention: iif.interventions){
        printIntervention(intervention);
    }

    // print the info associated with genotypes and their population
    for(auto map : iif.mapGenoToPop) {
        std::cout << "Genotype: " << map.first << " Population: " << map.second;
    } 
}


//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// PRIVATE FUNCTIONS ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

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



