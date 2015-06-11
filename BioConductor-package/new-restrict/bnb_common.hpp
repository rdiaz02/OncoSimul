#ifndef _BNB_COMMON_H_
#define _BNB_COMMON_H_

#include<Rcpp.h>
#include"debug_common.hpp"


// Simple custom exception for exceptions that lead to re-runs.
class rerunExcept: public std::runtime_error {
public:
  rerunExcept(const std::string &s) :
    std::runtime_error(s) {}
};


struct spParamsP {
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  double timeLastUpdate;
  std::multimap<double, int>::iterator pv;
  double absfitness; //convenient for Beerenwinkel
  int numMutablePos; //for mutator if need update of mutation
};


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
    //  DP2(spP.death); DP2(pM); 
    // print_spP(spP); 
    // for the running out of mutations
    // set a very tiny mutation prob, and if by chance it happens,
    // then mutate back to one self, so do not mutate
    // o.w., work out the math carefully
    throw std::range_error("pE.f: pE not finite");
  }
  return pE;
}

inline double pB_f_st(const double& pE,
			     const spParamsP& spP) {
  return (spP.birth * pE)/spP.death; 
}



inline void mapTimes_updateP(std::multimap<double, int>& mapTimes,
			     std::vector<spParamsP>& popParams,
			     const int index,
			     const double time) {
  // Update the map times <-> indices
  // First, remove previous entry, then insert.
  // But if we just created the species, nothing to remove from the map.
  if(popParams[index].timeLastUpdate > -1)
    mapTimes.erase(popParams[index].pv);
  popParams[index].pv = mapTimes.insert(std::make_pair(time, index));
}


inline void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes) {
  // we want minNextMutationTime and nextMutant
  nextMutant = mapTimes.begin()->second;
  minNextMutationTime = mapTimes.begin()->first;
}


inline void whichDrivers(int& totalPresentDrivers,
				std::string& strDrivers,
				const std::vector<int>& countByDriver){
  std::string comma = "";
  for(size_t i = 0; i < countByDriver.size(); ++i) {
    if(countByDriver[i] > 0) {
      strDrivers += (comma + std::to_string(i + 1)); 
      comma = ", ";
      ++totalPresentDrivers;
    }
  }
  if(totalPresentDrivers == 0) strDrivers = "NA";
}

// inline void whichDrivers(int& totalPresentDrivers,
// 				std::string& strDrivers,
// 				const std::vector<int>& countByDriver){
//   std::string comma = "";
//   DP2(countByDriver.size());
  
//   for(size_t i = 0; i < countByDriver.size(); ++i) {
//     DP1(" v 1");
//     DP2(i);
//     DP2( i + 1 );
//     DP1(" v 2");
//     int ii = i + 1;
//     DP2(ii);
//     DP1(" v 3");
//     std::string uu =  std::to_string(ii);
//     DP1(" v 4");
//     DP1(" v 5");
//     DP1(" v 6");
//     if(countByDriver[i] > 0) {
//       DP1("v 7");
//       DP2(countByDriver[i]);
//       DP1("v8");
//       DP1("v9");
//       int iii = i + 1;
//       DP1("v10");
//       std::string uuu = std::to_string(iii);
//       DP1("v11");
//       DP1("v12");
//       strDrivers += (comma + std::to_string(i + 1)); 
//       comma = ", ";
//       ++totalPresentDrivers;
//     }
//   }
//   if(totalPresentDrivers == 0) strDrivers = "NA";
// }

// inline void whichDrivers(int& totalPresentDrivers,
// 			 std::string& strDrivers,
// 			 const std::vector<int>& countByDriver){
//   DP1("w 1");
//   std::string comma = "";
//   std::string uu;
//   int ui;
//   for(size_t i = 0; i < countByDriver.size(); ++i) {
//     DP2(i);
//     DP2( i + 1);
//     ui = static_cast<int>(i) + 1;
//     DP1("after ui");
//     DP2(ui);
    
//     uu = std::to_string(ui);
//     DP1("w 2");
//     DP2(i);
//     if(countByDriver[i] > 0) {
//       DP1("w 3");
//       DP2(countByDriver[i]);
//        DP1("w 3 a");
//        std::string uu = std::to_string(i + 1);
//       strDrivers += (comma + std::to_string(i + 1));
//        DP1("w 3 aa");
//       comma = ", ";
//        DP1("w 3 ab");
//       ++totalPresentDrivers;
//       DP1("w3b");
//     }
//   }
//   DP1("w 4");
//   if(totalPresentDrivers == 0) strDrivers = "NA";
// }



inline void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
			       const std::vector<double>& sampleTotPopSize,
			       const std::vector<double>& sampleLargestPopSize,
			       const std::vector<double>& sampleLargestPopProp,
			       const std::vector<int>& sampleMaxNDr,
			       const std::vector<int>& sampleNDrLargestPop){

  for(size_t i = 0; i < sampleTotPopSize.size(); ++i) {
    perSampleStats(i, 0) = sampleTotPopSize[i];
    perSampleStats(i, 1) = sampleLargestPopSize[i]; // Never used in R FIXME: remove!!
    perSampleStats(i, 2) = sampleLargestPopProp[i]; // Never used in R
    perSampleStats(i, 3) = static_cast<double>(sampleMaxNDr[i]);
    perSampleStats(i, 4) = static_cast<double>(sampleNDrLargestPop[i]);
  }
}


void print_spP(const spParamsP& spP);

double pM_f_st(const double& t, const spParamsP& spP);

double ti_nextTime_tmax_2_st(const spParamsP& spP,
			     const double& currentTime,
			     const double& tSample,
			     int& ti_dbl_min,
			     int& ti_e3);

double Algo2_st(const spParamsP& spP,
		const double& ti);

double Algo3_st(const spParamsP& spP, const double& t);

void precissionLoss();

void init_tmpP(spParamsP& tmpParam);

double returnMFE(double& e1,
			const double& K,
		 const std::string& typeFitness);

void computeMcFarlandError(double& e1,
				  double& n_0,
				  double& n_1,
				  double& tps_0,
				  double& tps_1,
				  const std::string& typeFitness,
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
				  const int& mutatorGenotype,
			   const double& mu);

void updateRatesBeeren(std::vector<spParamsP>& popParams,
			      double& adjust_fitness_B,
			      const double& initSize,
			      const double& currentTime,
			      const double& alpha,
			      const double& totPopSize,
			      const int& mutatorGenotype,
		       const double& mu);


#endif

