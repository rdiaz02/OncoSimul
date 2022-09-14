#include "intervention.h"


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



