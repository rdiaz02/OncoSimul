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

#include "user_var.h"

bool isValid(std::string equation);
void parseAction(UserVarsInfo& uvif, 
                    Rule rule,  
                    double currentTime,
                    std::map<std::string, double> birthMap,
                    std::map<std::string, double> deathMap,
                    std::map<std::string, double> mutationMap);

UserVarsInfo addRule(UserVarsInfo uvif, Rule r){

    for(size_t i = 0; i<uvif.rules.size(); i++){
        if(compareRules(uvif.rules[i], r) == 0){
            Rcpp::Rcout << "There are two rules that are the same:";
            printRule(uvif.rules[i]);
            printRule(r);
            return uvif;
        }
    }

    uvif.rules.push_back(r);
    
    return uvif;
}

// function that creates a rule in memory
Rule createRule(std::string id, 
                        std::string condition, 
                        std::string action){
    Rule r;
    r.id = id;
    r.condition = condition;
    r.action = action;
    return r;
}

UserVarsInfo createUserVarsInfo(Rcpp::List rules,
                                Rcpp::List userVars,
                                const fitnessEffectsAll& fitnessEffects, 
                                const std::vector<spParamsP>& popParams, 
                                std::vector<Genotype> Genotypes){
   
    // we declare the variables needed
    UserVarsInfo uvif;
    Rule r;
    int totalEntries;
    //we use auxiliar variables to store the values from R
    std::vector<Rcpp::List> vectorOfList;
    std::vector<std::string> auxIDRule;
    std::vector<std::string> auxConditionRule;
    std::vector<std::string> auxActionRule;
    std::string auxName;
    double auxValue;


    for(int i=0; i<rules.length(); i++){
        vectorOfList.push_back(Rcpp::as<Rcpp::List>(rules[i]));
    }

    for(int i=0; i<rules.length(); i++){
        Rcpp::List listItem = vectorOfList[i];

        auxIDRule.push_back(Rcpp::as<std::string>(listItem["ID"]));
        auxConditionRule.push_back(Rcpp::as<std::string>(listItem["Condition"]));
        auxActionRule.push_back(Rcpp::as<std::string>(listItem["Action"]));
    }

    totalEntries = rules.length();

    // Now we dump the info in the structs created
    for(int i=0; i<totalEntries; i++){
        r = createRule(auxIDRule.at(i), auxConditionRule.at(i), auxActionRule.at(i));
        uvif = addRule(uvif, r);
    }

    // We fill up de user variables map
    vectorOfList = {};
    for(int i=0; i<userVars.length(); i++){
        vectorOfList.push_back(Rcpp::as<Rcpp::List>(userVars[i]));
    }
    for(int i=0; i<userVars.length(); i++){
        Rcpp::List listItem = vectorOfList[i];
        auxName = Rcpp::as<std::string>(listItem["Name"]);
        auxValue = Rcpp::as<double>(listItem["Value"]);
        uvif.userVars.insert({auxName, auxValue});
    }

    // mapping for the genes and their population
    uvif.mapGenoToPop = evalFVars(fitnessEffects, Genotypes, popParams, true);
    
    return uvif;
}

int compareRules(Rule r1, Rule r2){

    // it is not allowed to have 2 rules with the same ID.
    if(r1.id == r2.id){
        return 1;
    } 

    //TODO: check if conditions and actions are equivalent
    if (r1.condition != r2.condition) {
        return -1;
    }
    
    if (r1.action != r2.action){
        return -1;
    }
    // if they are equal in all aspects, then returns 0
    return 0;
}

UserVarsInfo destroyRule(UserVarsInfo uvif, Rule r){
    
    for (size_t i = 0; i<uvif.rules.size(); i++){
        if(compareRules(uvif.rules.at(i), r) == 0){
            uvif.rules.erase(uvif.rules.begin() + i);
            return uvif;
        }
    }

    return uvif;
}

void executeRules(UserVarsInfo& uvif,
                         double &currentTime,
                         std::map<std::string, double> birthMap,
                         std::map<std::string, double> deathMap,
                         std::map<std::string, double> mutationMap){

    // Now we add all the info needed for the symbol table so exprtk can operate 
    symbol_table_t symbol_table;
    double N = 0;
    for(auto& iterator : uvif.userVars) {
        symbol_table.add_variable(iterator.first, iterator.second);
    } 

    for(auto& iterator : uvif.mapGenoToPop) {
        symbol_table.add_variable(iterator.first, iterator.second);
        N += iterator.second;
    } 

    if(N == 0){
        throw std::invalid_argument("Total Population = 0.\n");
    }

    std::string name = "";
    std::string prefix = "b_";
    for(auto& iterator : birthMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    } 

    prefix = "d_";
    for(auto& iterator : deathMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    } 

    prefix = "m_";
    for(auto& iterator : mutationMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    
    } 

    double T = currentTime;

    symbol_table.add_constant("N", N);//We reserve N to total population size
    symbol_table.add_constant("T", T); //Pass current time to exprtk
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);

    //we iterate over the user-specified rules 
    for(auto& rule: uvif.rules){
        parser_t parser_condition;
        //if parser fails to compile, throws exception
        if (!parser_condition.compile(rule.condition, expression)){
            // error control, just in case the parsing it's not correct
            Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

            for (std::size_t i = 0; i < parser_condition.error_count(); ++i){
                typedef exprtk::parser_error::type error_t;
                error_t error = parser_condition.get_error(i);
                // FIXMEmaybe: Use warning or error to capture it easily in tests?
                REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                    i,
                    error.token.position,
                    exprtk::parser_error::to_str(error.mode).c_str(),
                    error.diagnostic.c_str(),
                    rule.condition.c_str());
            }
            std::string errorMessage = "The expression was imposible to parse.";
            throw std::invalid_argument(errorMessage);
        } else {
            //a condition is just a TRUE/FALSE condition
            if(expression.value() == 1){
                parser_t parser_action;
                //if parser fails to compile, throws exception
                if (!parser_action.compile(rule.action, expression)){
                    // error control, just in case the parsing it's not correct
                    Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

                    for (std::size_t i = 0; i < parser_action.error_count(); ++i){
                        typedef exprtk::parser_error::type error_t;
                        error_t error = parser_action.get_error(i);
                        // FIXMEmaybe: Use warning or error to capture it easily in tests?
                        REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                            i,
                            error.token.position,
                            exprtk::parser_error::to_str(error.mode).c_str(),
                            error.diagnostic.c_str(),
                            rule.action.c_str());
                    }
                    std::string errorMessage = "The expression was imposible to parse.";
                    throw std::invalid_argument(errorMessage);
                } else 
                    parseAction(uvif, rule, T, birthMap, deathMap, mutationMap);
            }
        }
    } 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////// PRIVATE FUNCTIONS ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void parseAction(UserVarsInfo& uvif, 
                    Rule rule, 
                    double currentTime,
                    std::map<std::string, double> birthMap,
                    std::map<std::string, double> deathMap,
                    std::map<std::string, double> mutationMap){
    
    // now we need to parse the "action" rule

    symbol_table_t symbol_table;
    double N = 0;
    for(auto& iterator : uvif.userVars) {
        symbol_table.add_variable(iterator.first, iterator.second);
    } 

    for(auto& iterator : uvif.mapGenoToPop) {
        symbol_table.add_variable(iterator.first, iterator.second);
        N += iterator.second;
    }

    std::string name = "";
    std::string prefix = "b_";
    for(auto& iterator : birthMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    } 

    prefix = "d_";
    for(auto& iterator : deathMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    } 

    prefix = "m_";
    for(auto& iterator : mutationMap) {
        name = "";
        name.insert(0, iterator.first);
        name.insert(0, prefix);
        symbol_table.add_variable(name, iterator.second);
    
    } 

    double T = currentTime;

    symbol_table.add_constant("N", N);//We reserve N to total population size
    symbol_table.add_constant("T", T); //Pass current time to exprtk
    // TODO: add other general variables such as f_x and n_x 
    symbol_table.add_constants();

    expression_t expression;
    expression.register_symbol_table(symbol_table);
    
    //Variables needed for tokeninzing
    std::vector <std::string> tokens;
    std::vector <std::string> assingmentTokens;
    std::stringstream check1(rule.action); 
    std::stringstream check2; 
    std::string intermediate; 

    // Tokenizing action ';' 
    while(getline(check1, intermediate, ';')){ 
        tokens.push_back(intermediate); 
    }

    for(std::string token: tokens){
        assingmentTokens.clear();
        check2.clear();
        check2 << token; 
        // We check that each of the assgnments of the rule's action is correctly specified
        if(!isValid(token)){
            throw std::runtime_error("The rule was wrongfully specified.\n");
        }

        // Tokenizing action '=' 
        while(getline(check2, intermediate, '=')){ 
            assingmentTokens.push_back(intermediate); 
        }
        std::string leftMostAction = assingmentTokens[0];
        leftMostAction.erase(std::remove(leftMostAction.begin(), leftMostAction.end(), ' '), leftMostAction.end());

        // we use the right-most side of the equation
        std::string rightMostAction = assingmentTokens[1];
        parser_t parser_ra;
        if (!parser_ra.compile(rightMostAction, expression)){
            // error control, just in case the parsing it's not correct
            Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;

            for (std::size_t i = 0; i < parser_ra.error_count(); ++i){
                typedef exprtk::parser_error::type error_t;
                error_t error = parser_ra.get_error(i);
                // FIXMEmaybe: Use warning or error to capture it easily in tests?
                REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
                    i,
                    error.token.position,
                    exprtk::parser_error::to_str(error.mode).c_str(),
                    error.diagnostic.c_str(),
                    rightMostAction.c_str());
            }
            std::string errorMessage = "The expression was imposible to parse.";
            throw std::invalid_argument(errorMessage);
        } else{
            double res = expression.value();

            std::map<std::string, double>::iterator it = uvif.userVars.find(leftMostAction); 
            if (it != uvif.userVars.end())
                it->second = res;
            else{
                std::string errorMessage = "The rule varies a non existing user variable.";
                throw std::invalid_argument(errorMessage);
            }
        }
    }
}

void printRule(Rule r){

    Rcpp::Rcout << r.id << " info:\n";
    Rcpp::Rcout << "\t Condition: " << r.condition << "\n";
    Rcpp::Rcout << "\t Action: " << r.action << "\n";
}

void printUserVarsInfo(UserVarsInfo uvif){

    // print user variables current values 
    for(auto map : uvif.userVars) {
            Rcpp::Rcout << "Name: " << map.first << " Value: " << map.second << "\n";
    }

    // print user variable rules
    for(auto rule: uvif.rules){
        printRule(rule);
    }    
}

// private function that checks that the equations specified by the user are correctly specified
// if no "=" or two "=" are found, returns false, if just one "=" is found, the returns true.
bool isValid(std::string equation){
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



