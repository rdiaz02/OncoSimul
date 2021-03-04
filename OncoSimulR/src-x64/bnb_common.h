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



#ifndef _BNB_COMMON_H_
#define _BNB_COMMON_H_

#include<Rcpp.h>
#include "common_classes.h"
#include "debug_common.h"
#include "new_restrict.h" // for the TypeModel enum


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


void print_mapTimes(std::multimap<double, int>& mapTimes);
  
void print_initMutant(const std::vector < std::vector<int> >& initMutant);

void print_Genotype(const Genotype& ge);

void mapTimes_updateP(std::multimap<double, int>& mapTimes,
		      std::vector<spParamsP>& popParams,
		      const int index,
		      const double time);


void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes);


void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
		 const std::vector<double>& sampleTotPopSize,
		 const std::vector<double>& sampleLargestPopSize,
		 const std::vector<double>& sampleLargestPopProp,
		 const std::vector<int>& sampleMaxNDr,
		 const std::vector<int>& sampleNDrLargestPop);



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

double returnMFE_new(double& en1,
		     const TypeModel typeModel);


void computeMcFarlandError_new(double& en1,
			       double& en1sc,
			       double& totPopSize_previous,
			       double& DA_previous,
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

void updateRatesFDFMcFarlandLog(std::vector<spParamsP>& popParams,
				const std::vector<Genotype>& Genotypes,
				const fitnessEffectsAll& fitnessEffects,
				double& adjust_fitness_MF,
				const double& K,
				const double& totPopSize,
				const double& currentTime);

void updateRatesMcFarlandLog_D(std::vector<spParamsP>& popParams,
			       double& adjust_fitness_MF,
			       const double& K,
			       const double& totPopSize);

void updateRatesFDFMcFarlandLog_D(std::vector<spParamsP>& popParams,
				  const std::vector<Genotype>& Genotypes,
				  const fitnessEffectsAll& fitnessEffects,
				  double& adjust_fitness_MF,
				  const double& K,
				  const double& totPopSize,
				  const double& currentTime);


void updateRatesFDFExp(std::vector<spParamsP>& popParams,
		       const std::vector<Genotype>& Genotypes,
		       const fitnessEffectsAll& fitnessEffects,
		       const double& currentTime);

void updateRatesFDFBozic(std::vector<spParamsP>& popParams,
			 const std::vector<Genotype>& Genotypes,
			 const fitnessEffectsAll& fitnessEffects,
			 const double& currentTime);


void updateBirthDeathRates(std::vector<spParamsP>& popParams,
			   const std::vector<Genotype>& Genotypes,
			   const fitnessEffectsAll& fitnessEffects,
			   double& adjust_fitness_MF,
			   const double& K,
			   const double& totPopSize,
			   const double& currentTime,
			   const TypeModel typeModel);

void detect_ti_duplicates(const std::multimap<double, int>& m,
			  const double ti,
			  const int spcies);



void message1(const int verbosity, const std::string message,
	      const int iteration, const double currentTime,
	      const unsigned int numSpecies,
	      const double totalPopulationSize,
	      const double timeNextPopSample,
	      const double minNextMutationTime);

void messageNewSpecies(const int verbosity,
		       const int iteration, 
		       const unsigned int numSpecies,
		       const int nextMutant);

void vvmessageNewSpecies(const int verbosity,
			 const unsigned int sp,
			 const Genotype& newGenotype,
			 const Genotype& parentGenotype,
			 const spParamsP& tmpParam,
			 const spParamsP& parentParam);

void messageSampling(const int verbosity,
		     const double tSample,
		     const double finalTime,
		     std::vector<spParamsP>& popParams);

void messagePostSampling(const int verbosity,
			 std::vector<spParamsP>& popParams,
			 const double totPopSize);

double setDummyMutationRate(std::vector<double> mu, const int verbosity);


void initPops(
	     unsigned int& numSpecies,
	     double& totPopSize,
	     int& outNS_i,
	     double& lastStoredSample,
	     std::vector<Genotype>& Genotypes,
	     std::vector<spParamsP>& popParams,
	     std::vector<Genotype>& genot_out,
	     std::vector<double>& popSizes_out,
	     std::vector<int>& index_out,
	     std::vector<double>& time_out,
	     std::vector<double>& sampleTotPopSize,
	     std::vector<double>& sampleLargestPopSize,
	     std::vector<int>& sampleMaxNDr,
	     std::vector<int>& sampleNDrLargestPop,
	     POM& pom,
	     std::mt19937& ran_gen,
	     const std::vector<std::vector<int> >& initMutant,
	     const std::vector<double>& initSize,
	     const fitnessEffectsAll& fitnessEffects,
	     const std::vector<double>& mu,
	     const fitnessEffectsAll& muEF,
	     const std::vector<int>& full2mutator,
	     const std::map<int, std::string>& intName,
	     const fitness_as_genes& genesInFitness,	     
	     const double& dummyMutationRate,
	     const double& K,
	     const double& death,
	     const double& currentTime,	     
	     const double& keepEvery,
	     const int& mutationPropGrowth,	     
	     const TypeModel typeModel,
	     const int& verbosity
	      );

#endif
