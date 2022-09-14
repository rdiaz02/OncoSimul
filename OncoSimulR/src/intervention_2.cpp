#include "intervention.h"




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

