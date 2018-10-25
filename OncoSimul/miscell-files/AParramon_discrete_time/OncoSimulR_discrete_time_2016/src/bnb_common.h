//     Copyright 2013, 2014, 2015 Ramon Diaz-Uriarte

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



#ifndef _BNB_COMMON_H_
#define _BNB_COMMON_H_

#include<Rcpp.h>
#include "common_classes.h"
#include "debug_common.h"
// #include "new_restrict.h" // for the TypeModel enum

// // Simple custom exception for exceptions that lead to re-runs.
// class rerunExcept: public std::runtime_error {
// public:
//   rerunExcept(const std::string &s) :
//     std::runtime_error(s) {}
// };


// struct spParamsP {
//   double popSize;
//   double birth;
//   double death;
//   double W;
//   double R;
//   double mutation; 
//   double timeLastUpdate;
//   std::multimap<double, int>::iterator pv;
//   double absfitness; //convenient for Beerenwinkel
//   int numMutablePos; //for mutator if need update of mutation
// };


inline void W_f_st(spParamsP& spP){
  spP.W = spP.death + spP.birth + spP.mutation;
}

inline void R_f_st(spParamsP& spP) {
  spP.R = sqrt( pow( spP.birth - spP.death, 2) + 
		( 2.0 * (spP.birth + spP.death) + 
		  spP.mutation) * spP.mutation );
}


inline double pE_f_st(double& pM, const spParamsP& spP){
  double pE = (spP.death * (1.0 - pM ) )/(spP.W - spP.death - spP.birth * pM );
  if( !std::isfinite(pE) ) {
    DP2(spP.death);  DP2(spP.birth); DP2(pM); DP2(spP.W);
    DP2(spP.mutation);
    std::string error_message = R"(pE.f: pE not finite.
      This is expected to happen when mutationPropGrowth = TRUE 
      and you have have an initMutant with death >> birth,
      as that inevitably leads to net birth rate of 0
      and mutation rate of 0)";
    throw std::range_error(error_message);
  }
  return pE;
}

inline double pB_f_st(const double& pE,
			     const spParamsP& spP) {
  return (spP.birth * pE)/spP.death; 
}

void mapTimes_updateP(std::multimap<double, int>& mapTimes,
			     std::vector<spParamsP>& popParams,
			     const int index,
			     const double time);


void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes);


void print_spP(const spParamsP& spP);

double pM_f_st(const double& t, const spParamsP& spP);

double ti_nextTime_tmax_2_st(const spParamsP& spP,
			     const double& currentTime,
			     const double& tSample,
			     int& ti_dbl_min,
			     int& ti_e3);

double Algo2_st(const spParamsP& spP,
		const double& ti,
		const int& mutationPropGrowth);

double Algo3_st(const spParamsP& spP, const double& t);

void precissionLoss();

void init_tmpP(spParamsP& tmpParam);

double returnMFE(double& e1,
		 // const double& K,
		 const std::string& typeFitness);

double returnMFE(double& e1,
		 // const double& K,
		 const TypeModel typeModel);

void computeMcFarlandError(double& e1,
			   double& n_0,
			   double& n_1,
			   double& tps_0,
			   double& tps_1,
			   const std::string& typeFitness,
			   const double& totPopSize,
			   const double& K);

void computeMcFarlandError(double& e1,
			   double& n_0,
			   double& n_1,
			   double& tps_0,
			   double& tps_1,
			   const TypeModel typeModel,
			   const double& totPopSize,
			   const double& K);



void updateRatesMcFarland(std::vector<spParamsP>& popParams,
				 double& adjust_fitness_MF,
				 const double& K,
			  const double& totPopSize);

void updateRatesMcFarlandLog(std::vector<spParamsP>& popParams,
				    double& adjust_fitness_MF,
				    const double& K,
			     const double& totPopSize);


void updateRatesMcFarland0(std::vector<spParamsP>& popParams,
				  double& adjust_fitness_MF,
				  const double& K,
				  const double& totPopSize,
				  const int& mutationPropGrowth,
			   const double& mu);

void updateRatesBeeren(std::vector<spParamsP>& popParams,
			      double& adjust_fitness_B,
			      const double& initSize,
			      const double& currentTime,
			      const double& alpha,
			      const double& totPopSize,
			      const int& mutationPropGrowth,
		       const double& mu);


void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
             const std::vector<double>& sampleTotPopSize,
             const std::vector<double>& sampleLargestPopSize,
             const std::vector<double>& sampleLargestPopProp,
             const std::vector<int>& sampleMaxNDr,
     const std::vector<int>& sampleNDrLargestPop);

#endif

