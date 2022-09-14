#include "intervention.h"


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
