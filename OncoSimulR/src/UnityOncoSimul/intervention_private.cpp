#include "intervention.h"


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
