//     Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

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


#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"

#include <limits>
#include <iostream>
#include <random>
#include <bitset>
#include <set>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <ctime>
#include <sys/time.h> 

#include <stdexcept>

using namespace Rcpp;
using std::vector;


// To track if mutation is really much smaller than birth/death
#define MIN_RATIO_MUTS
#ifdef MIN_RATIO_MUTS
// There is really no need for these to be globals?
// Unless I wanted to use them inside some function. So leave as globals.
double g_min_birth_mut_ratio = DBL_MAX;
double g_min_death_mut_ratio = DBL_MAX;
double g_tmp1 = DBL_MAX;
#endif


typedef std::bitset<64> Genotype64;


// Format of restrictTable
// - mutations in columns
// - first row, the number
// - second row, the number of dependencies
//   - rest of rows, the id of the dependency
//   - past number of dependencies: a -9
// In fact, the first row is redundant. Leave it, just in case.


// Genotypes: the first (or column 0) genotype is the all ceros.
// would not be needed, but makes Algo5 a lot simpler.

// But mutatedPos start at 0. 
// Will need to add 1 when plotting and analyzing with R.

// Ojo: typeCBN is going to be an int.
// But from R we pass a string, and that determined the integer.


static void fitness(spParamsP& tmpP,
		    const spParamsP& parentP,
		    const int& mutatedPos, 
		    Rcpp::IntegerMatrix restrictTable,
		    const std::string& typeCBN,
		    const Genotype64& newGenotype,
		    // const double& birthRate, 
		    const double& s,
		    // const double& death,
		    const int& numDrivers,
		    const std::string& typeFitness,
		    // const double& genTime,
		    // const double& adjust_fitness_B,
		    const double& sh){
  //const double& adjust_fitness_MF) {

  using namespace Rcpp;
  // Two pieces: split into two functions??
  //    - checking restrictions
  //    - returning actual fitness according


  int numDependencies;
  int sumDriversMet = 0;
  int sumDriversNoMet = 0;
  int sumDependenciesMet = 0;

  // set appropriate defaults. Change only needed stuff.
  tmpP.birth = parentP.birth;
  tmpP.death = parentP.death;
  tmpP.absfitness = parentP.absfitness;





  //      **** Are driver constraints met? ***


  // Two cases: same s, sh, sp or different ones. If same, return three
  // integers: sumDriversMet, sumDriversNoMet, sumPassengers.  If
  // different, return three vectors, filled with the non-zero
  // entries. These vectors then are combined as dictated by the fintness
  // functions.

  // If same single s, sh, sp: function takes three integers. O.w. it
  // takes three integer vectors.

  

  if(mutatedPos >= numDrivers) { //the new mutation is a passenger
    return;
  } else {
    for(int m = 0; m < numDrivers; ++m) {
      if( newGenotype[m] ) { // this m is mutated
	const Rcpp::IntegerMatrix::Column thisRestrict = 
	  restrictTable(_, m);
	numDependencies = thisRestrict[1];
	if(!numDependencies) { // this driver has no dependencies
	  sumDriversMet++;
#ifdef DEBUGZ
	  Rcpp::Rcout << "\n No dependencies:  ";
	  DP2(sumDriversMet);
#endif

	}
	else {
	  sumDependenciesMet = 0;
	  for(int i = 2; i < (2 + numDependencies); i++) {
	    sumDependenciesMet += newGenotype[ thisRestrict[i] ];
	  }
	  if( ( (typeCBN == "Multiple") && (sumDependenciesMet) ) ||
	      ( (typeCBN == "CBN") && (sumDependenciesMet == numDependencies) )) {
	    sumDriversMet++;   
	  } else {
	    sumDriversNoMet++;
	  }
	}
      }
    }
  }

#ifdef DEBUGZ
  DP2(sumDriversMet);
  DP2(sumDriversNoMet);
  DP2(sh);
  DP2(typeFitness);
#endif

  // if sh < 0 : we do not allow any unment dependencies.  
  // if sh = 0: no penalty for unmet dependencies


  // FIXME: why not just pass the birth and death rates, and combine them
  // in arbitrary ways? Might even allow to pass on death and birth rates
  // from R. Only need care when any are density dependent.


  // Beware: doing it this way with Bozic1 is kind of questionable because
  // if birth = 1, there is no immediate extinction. In fact, there never
  // is.
  if((sh < 0) && sumDriversNoMet) {
    tmpP.absfitness = 0.0;
    tmpP.death = 1.0;
    tmpP.birth = 0.0; // this is what really matters so that
    // the pop does not get added.
    // Line with comment "fitness is 0"
  } else {
    if(typeFitness == "bozic1") {
      tmpP.death = pow( 1.0 - s, sumDriversMet) * 
	pow( 1.0 + sh, sumDriversNoMet);
      tmpP.birth = 1.0;
    // } else if (typeFitness == "bozic2") {
    //   double pp = pow( 1.0 - s, sumDriversMet) * 
    // 	pow( 1.0 + sh, sumDriversNoMet);
    //   tmpP.birth = (1.0/genTime) * (1.0 - 0.5 * pp );
    //   tmpP.death = (0.5/genTime) * pp;
    // } else if(typeFitness == "beerenwinkel") {
    //   // like Datta et al., 2013
    //   tmpP.absfitness = pow(1.0 + s, sumDriversMet) * 
    // 	pow( 1.0 - sh, sumDriversNoMet);
    //   tmpP.birth = adjust_fitness_B * tmpP.absfitness;
    // } else if(typeFitness == "mcfarland0") {
    //   tmpP.absfitness = pow(1.0 + s, sumDriversMet) / 
    // 	pow( 1.0 + sh, sumDriversNoMet);
    //   tmpP.birth = adjust_fitness_MF * tmpP.absfitness;
    // } else if(typeFitness == "mcfarland") {
    //   tmpP.birth = pow(1.0 + s, sumDriversMet) / 
    // 	pow( 1.0 + sh, sumDriversNoMet);
    } else if(typeFitness == "mcfarlandlog") {
      tmpP.birth = pow(1.0 + s, sumDriversMet) / 
	pow( 1.0 + sh, sumDriversNoMet);
    } else if (typeFitness == "exp") { 
      // Also like Datta et al., 2013 An additional driver gene mutation
      // increases a cell’s fitness by a factor of (1+sd), whereas an
      // additional housekeeper gene mutation decreases fitness by a
      // factor of (1-sh) and the effect of multiple mutations is
      // multiplicative
      tmpP.birth = pow(1.0 + s, sumDriversMet) * 
	pow( 1.0 - sh, sumDriversNoMet);

#ifdef DEBUGZ
      double posi = pow(1.0 + s, sumDriversMet);
      double negi = pow( 1.0 - sh, sumDriversNoMet);
      DP2(posi);
      DP2(negi);
#endif

    } // else if (typeFitness == "log") {
    //   tmpP.birth = birthRate+ s * log1p(sumDriversMet) - 
    // 	sh * log(1 + sumDriversNoMet);
    // } else { // linear
    //   tmpP.birth = birthRate + s * static_cast<double>(sumDriversMet) - 
    // 	sh * static_cast<double>(sumDriversNoMet);
    // } 
  }
}
// Notice: if restriction is 3 -> 4 -> 5
// and one has 5 and 4, only 4 is unmet. Beware of that.
// So we talk about the immediate dependency or restriction.
// Not the whole transitive closure.

// When birth == 0, popSize should become 0 immediately. 
// No evaluation through random numbers, etc.
// This is how we do it.

// How small can they get?
// d1 <- function(s, mut) { (1 - s)^mut}
// d2 <- function(s, mut) {  (0.5/4) * ((1 - s)^mut) }


// limited benchmarks suggest the following is slower
// static inline void new_sp_bitset2(unsigned int& sp, const Genotype64& newGenotype,
// 			   const std::vector<Genotype64>& Genotypes) {
//   sp = std::distance(Genotypes.begin(),
// 		     std::find(Genotypes.begin(), 
// 			       Genotypes.end(), newGenotype));
// }


static inline void new_sp_bitset(unsigned int& sp, const Genotype64& newGenotype,
			  const std::vector<Genotype64>& Genotypes) {
  sp = 0;

  for(sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] )
      break;
  }
}



static void getMutatedPos_bitset(int& mutatedPos, int& numMutablePosParent,
				 //gsl_rng *r,
				 std::mt19937& ran_generator, 
				 std::vector<int>& mutablePos,
				 const Genotype64& nextMutantGenotype,
				 // const int& nextMutant,
				 // const std::vector<Genotype64>& Genotypes,
				 const int& numGenes) {
  // We want mutatedPos and numMutablePosParent

  // Note: impossible to have a second recorded mutation in
  // the same gene.  

  // Remember numMutablePosParent is the number of mutable positions in
  // the parent!  so after mutation is one less, but we do not decrease it
  // here.
  
  numMutablePosParent = 0;
  for(int i = 0; i < numGenes; ++i) {
    if( !nextMutantGenotype.test(i) ) { 
      mutablePos[numMutablePosParent] = i;
      ++numMutablePosParent;
    }
  }
  
  if(numMutablePosParent > 1) {
    std::uniform_int_distribution<int> unif(0, numMutablePosParent - 1);
    mutatedPos = mutablePos[unif(ran_generator)];
  } else {
    mutatedPos = mutablePos[0];
  } 

  // if(numMutablePosParent > 1) {
  //   mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePosParent)];
  // } else {
  //   mutatedPos = mutablePos[0];
  // } 


#ifdef DEBUGV
      Rcpp::Rcout << "\n numMutablePosParent = " << numMutablePosParent;
      Rcpp::Rcout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

  // if(numMutablePos > 1) {
  //   mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
  // } else if (numMutablePos == 1) {
  //   mutatedPos = mutablePos[0];
  // } else {
  //   // Should never happen, as mutation = 0 if no mutable positions.
  //   throw std::out_of_range("Algo5: run out of mutable places!!??");
  // }

}


static void remove_zero_sp_v7(std::vector<int>& sp_to_remove,
			      std::vector<Genotype64>& Genotypes,
			      std::vector<spParamsP>& popParams,
			      std::multimap<double, int>& mapTimes) {
  //  here("entering remove_zero_sp_v7");
  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();
  std::vector<int>::reverse_iterator r = sp_to_remove.rbegin();
  int remove_this;
  // for(r = sp_to_remove.rbegin(); r != sp_to_remove.rend(); ++r) {
  while(r != sp_to_remove.rend() ) {
    remove_this = *r;
    mapTimes.erase(popParams[remove_this].pv);
    popParams.erase(popParams_begin + remove_this);
    Genotypes.erase(Genotypes_begin + remove_this);
    ++r;
  }
  // here("exiting remove_zero_sp_v7");

}

static inline int count_NDrivers(const Genotype64& Genotype,
				 const int& NumDrivers) {
  int totalDr = 0;
  for(int i = 0; i < NumDrivers; ++i)
    totalDr += Genotype[i];
  return totalDr;
}

static void totPopSize_and_fill_out_crude_P(int& outNS_i,
					    double& totPopSize, 
					    double& lastStoredSample,
					    std::vector<Genotype64>& genot_out,
					    //std::vector<unsigned long long>& sp_id_out,
					    std::vector<double>& popSizes_out,
					    std::vector<int>& index_out,
					    std::vector<double>& time_out,
					    std::vector<double>& sampleTotPopSize,
					    std::vector<double>& sampleLargestPopSize,
					    std::vector<int>& sampleMaxNDr,
					    std::vector<int>& sampleNDrLargestPop,
					    bool& simulsDone,
					    bool& reachDetection,
					    int& lastMaxDr,
					    double& done_at,
					    const std::vector<Genotype64>& Genotypes,
					    const std::vector<spParamsP>& popParams, 
					    const double& currentTime,
					    const int& NumDrivers,
					    const double& keepEvery,
					    const double& detectionSize,
					    const double& finalTime,
					    // const double& endTimeEvery,
					    const int& detectionDrivers,
					    const int& verbosity,
					    const double& minDetectDrvCloneSz,
					    const double& extraTime,
					    const double& fatalPopSize = 1e15) {
  // Fill out, but also compute totPopSize
  // and return sample summaries for popsize, drivers.
  
  // This determines if we are done or not by checking popSize, number of
  // drivers, etc
  
  // static int lastMaxDr = 0; // preserves value across calls to Algo5 from R.
  // so can not use it.
  bool storeThis = false;
  totPopSize = 0.0;
  
   // DP2(lastMaxDr);
  // DP2(detectionDrivers);
  // DP2(currentTime);
  // DP2((lastStoredSample + endTimeEvery));
  // DP2(detectionSize);

  // this could all be part of popSize_over_m_dr, with a better name
  int tmp_ndr = 0;
  int max_ndr = 0;
  double popSizeOverDDr = 0.0;

  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
    tmp_ndr = count_NDrivers(Genotypes[i], NumDrivers);
    if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
    if(tmp_ndr >= detectionDrivers) popSizeOverDDr += popParams[i].popSize;
  }
  lastMaxDr = max_ndr;

  
  if (keepEvery < 0) {
    storeThis = false;
  } else if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }


  // if( (totPopSize >= detectionSize) ||
  //     ( (lastMaxDr >= detectionDrivers) &&
  //       (popSizeOverDDr >= minDetectDrvCloneSz) ) ) {
  //   simulsDone = true;
  //   reachDetection = true;
  // }

  if(extraTime > 0) {
    if(done_at <  0) {
      if( (totPopSize >= detectionSize) ||
	  ( (lastMaxDr >= detectionDrivers) &&
	    (popSizeOverDDr >= minDetectDrvCloneSz) ) ) {
	done_at = currentTime + extraTime;
      }
    } else if (currentTime >= done_at) {
	simulsDone = true;
	reachDetection = true; 
      }
  } else if( (totPopSize >= detectionSize) ||
	     ( (lastMaxDr >= detectionDrivers) &&
	       (popSizeOverDDr >= minDetectDrvCloneSz) ) ) {
    simulsDone = true;
    reachDetection = true; 
  }
  
 

    
  // This is no longer used.
  // // Beware: this can lead to never stopping if
  // // decreases in popSize or drivers

  // // Logic: if a period k you meet any condition, recheck again at k +
  // // endTimeEvery, and if conditions met exit. Prevents exiting if you
  // // reach the cancer state almost by chance. But this is way too
  // // paranoid. The idea is of application mainly for McF and Beeren
  // // models, so we do not bail out as soon as just a single cell with one
  // // new driver. But this makes things very slow.

  // // Thus, never pass an endTimeEvery > 0, but use detectionDrivers = 1 +
  // // intended final Drivers.

  // // FIXME
  // // Ideally, we would check, for McFL, that popsize of the pop with
  // // required number of drivers is at least, say, > initSize.
  // // But that is not trivial, as requires more accounting. Do later.

  
  // if(endTimeEvery > 0) {
  //   if(done_at <= 0 ) {
  //     if( (totPopSize >= detectionSize) ||
  // 	   (lastMaxDr >= detectionDrivers)  )
  // 	done_at = currentTime + endTimeEvery;
  //   } else if (currentTime >= done_at) {
  //     if( (totPopSize >= detectionSize) ||
  // 	  (lastMaxDr >= detectionDrivers)  ) {
  // 	simulsDone = true;
  // 	reachDetection = true;
  //     }
  //     else
  // 	done_at = -9;
  //   }
  // } else if( (totPopSize >= detectionSize) ||
  // 	     (lastMaxDr >= detectionDrivers) )  {
    
  //     simulsDone = true;
  //     reachDetection = true;
  // }

  
  if(totPopSize >= fatalPopSize) {
    Rcpp::Rcout << "\n\totPopSize > " << fatalPopSize
		<<". You are likely to loose precision and run into numerical issues\n";
       }
  
  if(simulsDone)
    storeThis = true;


  if( storeThis ) {
    lastStoredSample = currentTime;
    outNS_i++;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    
    time_out.push_back(currentTime);
    
    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);
      
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = count_NDrivers(Genotypes[i], NumDrivers);
      }
    }
    sampleTotPopSize.push_back(totPopSize);
    sampleLargestPopSize.push_back(l_pop_s);
    sampleMaxNDr.push_back(max_ndr);
    sampleNDrLargestPop.push_back(ndr_lp);
  } 


  
  // if( storeThis ) {
  //   lastStoredSample = currentTime;
  //   outNS_i++;
  //   int tmp_ndr = 0;
  //   int max_ndr = 0;
  //   int ndr_lp = 0;
  //   double l_pop_s = 0.0;
    
  //   time_out.push_back(currentTime);
    
  //   for(size_t i = 0; i < popParams.size(); ++i) {
  //     genot_out.push_back(Genotypes[i]);
  //     popSizes_out.push_back(popParams[i].popSize);
  //     index_out.push_back(outNS_i);
  //     // I have to repeat the counting of drivers here.
  //     tmp_ndr = count_NDrivers(Genotypes[i], NumDrivers); 
  //     if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
  //     if(popParams[i].popSize > l_pop_s) {
  // 	l_pop_s = popParams[i].popSize;
  // 	ndr_lp = tmp_ndr;
  // 	// ndr_lp = count_NDrivers(Genotypes[i], NumDrivers); 
  //     }
  //     // lastMaxDr = max_ndr; // and this should have been out of the
  //     // popParams.size() loop
  //   }
  //   // lastMaxDr = max_ndr;
  //   sampleTotPopSize.push_back(totPopSize);
  //   sampleLargestPopSize.push_back(l_pop_s);
  //   sampleMaxNDr.push_back(max_ndr);
  //   sampleNDrLargestPop.push_back(ndr_lp);
  // }//  else if (keepEvery < 0) {
  //   // FIXME keepEvery
  //   // must keep track of results to bail out

  //   // FIXME counting max drivers should be done always, like counting
  //   // totPopSize.
    
  //   int tmp_ndr = 0;
  //   int max_ndr = 0;
   
  //   for(size_t i = 0; i < popParams.size(); ++i) {
  //     tmp_ndr = count_NDrivers(Genotypes[i], NumDrivers);
  //     if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
  //     // lastMaxDr = max_ndr;
  //   }
  //   lastMaxDr = max_ndr;
  // }

  
    
  
  if( !std::isfinite(totPopSize) ) {
    throw std::range_error("totPopSize not finite");
  }
  if( std::isnan(totPopSize) ) {
    throw std::range_error("totPopSize is NaN");
  }
  
  if(totPopSize > (4.0 * 1e15)) {
    if(verbosity > 0)
      Rcpp::Rcout << "\nWARNING: popSize > 4e15. Likely loss of precission\n";
  }
}

// FIXME: I might want to return the actual drivers in each period
// and the actual drivers in the population with largest popsize
// Something like what we do now with whichDrivers
// and count_NumDrivers


// static inline void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
// 			       const std::vector<double>& sampleTotPopSize,
// 			       const std::vector<double>& sampleLargestPopSize,
// 			       const std::vector<double>& sampleLargestPopProp,
// 			       const std::vector<int>& sampleMaxNDr,
// 			       const std::vector<int>& sampleNDrLargestPop){

//   for(size_t i = 0; i < sampleTotPopSize.size(); ++i) {
//     perSampleStats(i, 0) = sampleTotPopSize[i];
//     perSampleStats(i, 1) = sampleLargestPopSize[i]; // Never used in R FIXME: remove!!
//     perSampleStats(i, 2) = sampleLargestPopProp[i]; // Never used in R
//     perSampleStats(i, 3) = static_cast<double>(sampleMaxNDr[i]);
//     perSampleStats(i, 4) = static_cast<double>(sampleNDrLargestPop[i]);
//   }
// }

inline void reshape_to_outNS(Rcpp::NumericMatrix& outNS,
			     const std::vector<unsigned long long>& uniqueGenotV,
			     const std::vector<unsigned long long>& genot_out_ul,
			     const std::vector<double>& popSizes_out,
			     const std::vector<int>& index_out,
			     const std::vector<double>& time_out){
  
  std::vector<unsigned long long>::const_iterator fbeg = uniqueGenotV.begin();
  std::vector<unsigned long long>::const_iterator fend = uniqueGenotV.end();

  int column;

  for(size_t i = 0; i < genot_out_ul.size(); ++i) {
    column = std::distance(fbeg, lower_bound(fbeg, fend, genot_out_ul[i]) );
    // here("   looping over i ");
    outNS(index_out[i], column + 1) =  popSizes_out[i];
  }

  for(size_t j = 0; j < time_out.size(); ++j)
    outNS(j, 0) = time_out[j];
}

static inline void find_unique_genotypes(std::set<unsigned long long>& uniqueGenotypes,
				  const std::vector<unsigned long long>& genot_out_l) {
  for(size_t i = 0; i < genot_out_l.size(); ++i) 
    uniqueGenotypes.insert( genot_out_l[i] );
}

static inline void genot_out_to_ullong(std::vector<unsigned long long>& go_l,
			       const std::vector<Genotype64>& go) {
  for(size_t i = 0; i < go.size(); ++i)
    go_l[i] = go[i].to_ullong();
}


static inline void uniqueGenotypes_to_vector(std::vector<unsigned long long>& ugV,
				      const std::set<unsigned long long>& uniqueGenotypes) {
  ugV.assign(uniqueGenotypes.begin(), uniqueGenotypes.end() );
}


static inline void create_returnGenotypes(Rcpp::IntegerMatrix& returnGenotypes,
					  const int& numGenes,
					  const std::vector<unsigned long long>& uniqueGenotypesV){
  // In C++, as the original were bitsets, pos 0 is at the right
  // In R, pos 0 is at the left

  for(size_t i = 0; i < uniqueGenotypesV.size(); ++i) {
    Genotype64 tmpbs(uniqueGenotypesV[i]);
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(j, i) = tmpbs[j];
    }
  }
}

// FIXME: change this, now that we keep a count of drivers?
static inline void count_NumDrivers(int& maxNumDrivers, 
				    std::vector<int>& countByDriver,
				    Rcpp::IntegerMatrix& returnGenotypes,
				    const int& numDrivers){
  //				    Rcpp::IntegerVector& totDrivers){

  // At the end, over all genotypes
  // Using a returnGenotypes object as input
  // Redundant? 
  maxNumDrivers = 0;
  int tmpdr = 0;
  
  for(int j = 0; j < returnGenotypes.ncol(); ++j) {
    tmpdr = 0;
    for(int i = 0; i < numDrivers; ++i) {
      tmpdr += returnGenotypes(i, j);
      countByDriver[i] += returnGenotypes(i, j);
    }
    // totDrivers(j) = tmpdr;
    if(tmpdr > maxNumDrivers) maxNumDrivers = tmpdr;
  }
}

static inline void whichDrivers(int& totalPresentDrivers,
				std::string& strDrivers,
				const std::vector<int>& countByDriver){
  std::string comma = "";
  for(size_t i = 0; i < countByDriver.size(); ++i) {
    if(countByDriver[i] > 0) {
#ifdef _WIN32  
      strDrivers += (comma + SSTR(i + 1));
#endif
#ifndef _WIN32
      strDrivers += (comma + std::to_string(i + 1)); //SSTR(i + 1));
#endif      
      comma = ", ";
      ++totalPresentDrivers;
    }
  }
  if(totalPresentDrivers == 0) strDrivers = "NA";
}


static void sample_all_pop_P(std::vector<int>& sp_to_remove,
			     std::vector<spParamsP>& popParams,
			     // next only used with DEBUGV
			     const std::vector<Genotype64>& Genotypes,
			     const double& tSample,
			     const int& mutationPropGrowth){

  // here("entering sample_all_pop_P");
  // currentTime = tSample;
  sp_to_remove.clear();
  // sp_to_remove.push_back(0);
  // sp_to_remove[0] = 0;

  for(size_t i = 0; i < popParams.size(); i++) {
    //STOPASSERT(popParams[i].Flag == false);
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
	      << "\n      Genotype = " << Genotypes[i]
	      << "\n      sp_id = " << Genotypes[i].to_ullong() // sp_id[i]  
	      << "\n      pre-update popSize = " 
	      << popParams[i].popSize 
	      << "\n      time of sample = " << tSample 
	      << "\n      popParams[i].timeLastUpdate = " 
	      << popParams[i].timeLastUpdate 
	      << ";\n     t for Algo2 = " 
	      << tSample - popParams[i].timeLastUpdate 
	      << " \n     species R " << popParams[i].R
	      << " \n     species W " << popParams[i].W
	      << " \n     species death " << popParams[i].death
	      << " \n     species birth " << popParams[i].birth;
    // << " \n     species nextMutationTime " 
    // << popParams[i].nextMutationTime;
#endif

    // Account for forceSampling. When 
    // forceSampling, popSize for at least one species
    // was updated in previous loop, so we skip that one
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize = 
	Algo2_st(popParams[i], tSample, mutationPropGrowth);
    }
    if( popParams[i].popSize <=  0.0 ) {
      // this i has never been non-zero in any sampling time
      // eh??

      // If it is 0 here, remove from _current_ population. Anything that
      // has had a non-zero size at sampling time is preserved (if it
      // needs to be preserved, because it is keepEvery time).

      // sp_to_remove[0]++;
      sp_to_remove.push_back(i);
      // sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
      Rcpp::Rcout << "\n\n     Removing species i = " << i 
		<< " with sp_id = " << Genotypes[i].to_ullong(); //sp_id[i];
#endif
    } // else {
    //   popParams[i].Flag = true;
    // }
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
  }
  // here("exiting sample_all_pop_P");
}


static void innerBNB(const int& numGenes,
		     const double& initSize,
		     const double& K,
		     // const double& alpha,
		     const std::string& typeCBN,
		     //		     const double& genTime,
		     const std::string& typeFitness,
		     const int& mutationPropGrowth,
		     const double& mu,
		     const double& sh,
		     const double& s,
		     const double& death,
		     // const double& birthRate,
		     const double& keepEvery,
		     const double& sampleEvery,		     
		     const int& numDrivers,
		     const int& initMutant,
		     const time_t& start_time,
		     const double& maxWallTime,
		     const double& finalTime,
		     const double& detectionSize,
		     // const double& endTimeEvery,
		     const int& detectionDrivers,
		     const double& minDetectDrvCloneSz,
		     const double& extraTime,
		     const int& verbosity,
		     double& totPopSize,
		     double& e1,
		     double& n_0,
		     double& n_1,
		     double& ratioForce,
		     double& currentTime,
		     int& speciesFS,
		     int& outNS_i,
		     int& iter,
		     std::vector<Genotype64>& genot_out,
		     std::vector<double>& popSizes_out,
		     std::vector<int>& index_out,
		     std::vector<double>& time_out,
		     std::vector<double>& sampleTotPopSize,
		     std::vector<double>& sampleLargestPopSize,
		     std::vector<int>& sampleMaxNDr,
		     std::vector<int>& sampleNDrLargestPop,
		     bool& reachDetection,
		     std::mt19937& ran_generator,
		     double& runningWallTime,
		     bool& hittedWallTime,
		     Rcpp::IntegerMatrix restrictTable) {
		     //bool& anyForceRerunIssues
  //  if(numRuns > 0) {

  double dummyMutationRate = std::max(mu/1000, 1e-13);
  // double dummyMutationRate = 1e-10;
  // ALWAYS initialize this here, or reinit or rezero
  genot_out.clear();
  popSizes_out.clear();
  index_out.clear();
  time_out.clear();
  totPopSize = 0.0;
  sampleTotPopSize.clear();
  currentTime = 0.0;
  iter = 0;

  outNS_i = -1;

  sampleTotPopSize.clear();
  sampleLargestPopSize.clear();
  sampleMaxNDr.clear();
  sampleNDrLargestPop.clear();
  // end of rezeroing.

  
  // }
  // anyForceRerunIssues = false;
  
  bool forceSample = false;
  bool simulsDone = false;
  double lastStoredSample = 0.0;


  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double timeNextPopSample;
  double tSample;

  int nextMutant;
  unsigned int numSpecies = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  //int indexMutatedPos = 0;

  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  //int timeL = 1000;
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);

  // those to update
  int to_update = 1; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  Genotypes[0].reset();
  
  std::vector<spParamsP> popParams(1);

      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 1 ";
      // print_spP(popParams[0]);
      // // end debug
  
  
  const int sp_per_period = 5000;
  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 01 ";
      // print_spP(popParams[0]);
      // // end debug

  
  spParamsP tmpParam; 
  init_tmpP(tmpParam);
  init_tmpP(popParams[0]);
  popParams[0].popSize = initSize;
  totPopSize = initSize;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10000 ";
      // print_spP(popParams[0]);
      // // end debug


  
  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

    // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;

  int ti_dbl_min = 0;
  int ti_e3 = 0;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10002 ";
      // print_spP(popParams[0]);
      // // end debug


  
  // // Beerenwinkel
  // double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  //McFarland
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  // for McFarland error
  e1 = 0.0;
  n_0 = 0.0;
  n_1 = 0.0;
  double tps_0, tps_1; 
  tps_0 = totPopSize;
  tps_1 = totPopSize;


      // // FIXME debug
      // Rcpp::Rcout << "\n popSize[0]  at 10004 ";
      // print_spP(popParams[0]);
      // // end debug

  
  
  int lastMaxDr = 0;
  double done_at = -9;

#ifdef MIN_RATIO_MUTS
  g_min_birth_mut_ratio = DBL_MAX;
  g_min_death_mut_ratio = DBL_MAX;
  g_tmp1 = DBL_MAX;
#endif

      // // FIXME debug
      // Rcpp::Rcout << " popSize[0]  at 1b ";
      // print_spP(popParams[0]);
      // // end debug

    // This long block, from here to X1, is ugly and a mess!
  // This is what takes longer to figure out whenever I change
  // anything. FIXME!!

  if(initMutant >= 0)
    throw std::invalid_argument("initMutant no longer allowed. But in R code.");
  // if(initMutant >= 0) {
  //   popParams[0].numMutablePos = numGenes - 1;
  //   Genotypes[0].set(initMutant);
  //   // if(typeFitness == "beerenwinkel") {
  //   //   popParams[0].death = 1.0; //note same is in McFarland.
  //   //   // But makes sense here; adjustment in beerenwinkel is via fitness
      
  //   //   // initialize to prevent birth/mutation warning with Beerenwinkel
  //   //   // when no mutator. O.w., the defaults
  //   //   if(!mutationPropGrowth)
  //   // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
  //   //   popParams[0].absfitness = 1.0 + s;
  //   //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
  //   // 			currentTime, alpha, initSize, 
  //   // 			mutationPropGrowth, mu);
  //   // } else if(typeFitness == "mcfarland0") {
  //   //   // death equal to birth of a non-mutant.
  //   //   popParams[0].death = log1p(totPopSize/K); // log(2.0), except rare cases
  //   //   if(!mutationPropGrowth)
  //   // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
  //   //   popParams[0].absfitness = 1.0 + s;
  //   //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
  //   // 			    totPopSize,
  //   // 			    mutationPropGrowth, mu);
  //   // } else if(typeFitness == "mcfarland") {
  //   //   popParams[0].death = totPopSize/K;
  //   //   popParams[0].birth = 1.0 + s;
  //   // } else if(typeFitness == "mcfarlandlog") {
  //   if(typeFitness == "mcfarlandlog") {      
  //     popParams[0].death = log1p(totPopSize/K);
  //     popParams[0].birth = 1.0 + s;
  //   } else if(typeFitness == "bozic1") {
  //     tmpParam.birth =  1.0;
  //     tmpParam.death = -99.9;
  //   // } else if (typeFitness == "bozic2") {
  //   //   tmpParam.birth =  -99;
  //   //   tmpParam.death = -99;
  //   } else if (typeFitness == "exp") {
  //     tmpParam.birth =  -99;
  //     tmpParam.death = death;
  //   } // else { // linear or log
  //   //   tmpParam.birth =  -99;
  //   //   tmpParam.death = death;
  //   // } 
  //   // if( (typeFitness != "beerenwinkel") && (typeFitness != "mcfarland0") 
  //   // 	&& (typeFitness != "mcfarland") && (typeFitness != "mcfarlandlog")) // wouldn't matter
  //   //   fitness(popParams[0], tmpParam, initMutant, restrictTable,
  //   // 	      typeCBN, Genotypes[0], birthRate, s, numDrivers, 
  //   // 	      typeFitness, genTime, adjust_fitness_B, sh,
  //   // 	      adjust_fitness_MF);

  //   if( typeFitness != "mcfarlandlog") // wouldn't matter
  //     fitness(popParams[0], tmpParam, initMutant, restrictTable,
  // 	      typeCBN, Genotypes[0], // birthRate,
  // 	      s, numDrivers, 
  // 	      typeFitness, // genTime, adjust_fitness_B,
  // 	      sh);
  //   // adjust_fitness_MF);
  //   // we pass as the parent the tmpParam; it better initialize
  //   // everything right, or that will blow. Reset to init
  //   init_tmpP(tmpParam);
  // } else {
    popParams[0].numMutablePos = numGenes;
    // if(typeFitness == "beerenwinkel") {
    //   popParams[0].death = 1.0;
    //   // initialize to prevent birth/mutation warning with Beerenwinkel
    //   // when no mutator. O.w., the defaults
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
    // 			currentTime, alpha, initSize, 
    // 			mutationPropGrowth, mu);
    // } else if(typeFitness == "mcfarland0") {
    //   popParams[0].death = log1p(totPopSize/K);
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
    // 			    totPopSize,
    // 			    mutationPropGrowth, mu);
    // } else if(typeFitness == "mcfarland") {
    //   popParams[0].birth = 1.0;
    //   popParams[0].death = totPopSize/K;
    //   // no need to call updateRates
    // } else if(typeFitness == "mcfarlandlog") {
    if(typeFitness == "mcfarlandlog") {      
      popParams[0].birth = 1.0;
      popParams[0].death = log1p(totPopSize/K);
      // no need to call updateRates
    } else if(typeFitness == "bozic1") {
      popParams[0].birth = 1.0;
      popParams[0].death = 1.0;
    // } else if (typeFitness == "bozic2") {
    //   popParams[0].birth = 0.5/genTime;
    //   popParams[0].death = 0.5/genTime;
    } else if (typeFitness == "exp") {
      popParams[0].birth = 1.0;
      popParams[0].death = death;
    } // else { // linear or log
    //   popParams[0].birth = birthRate;
    //   popParams[0].death = death;
    // }
    //  }


  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 2 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  
  // these lines (up to, and including, R_F_st)
  // not needed with mcfarland0 or beerenwinkel
  if(mutationPropGrowth)
    popParams[0].mutation = mu * popParams[0].birth * popParams[0].numMutablePos;
  else
    popParams[0].mutation = mu * popParams[0].numMutablePos;

  W_f_st(popParams[0]);
  R_f_st(popParams[0]);

  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 3 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  // X1: end of mess of initialization block

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  if( keepEvery > 0 ) {
    // We keep the first ONLY if we are storing more than one.
    outNS_i++;
    time_out.push_back(currentTime);
    
    genot_out.push_back(Genotypes[0]);
    popSizes_out.push_back(popParams[0].popSize);
    index_out.push_back(outNS_i);

    sampleTotPopSize.push_back(popParams[0].popSize);
    sampleLargestPopSize.push_back(popParams[0].popSize);
    sampleMaxNDr.push_back(count_NDrivers(Genotypes[0], numDrivers));
    sampleNDrLargestPop.push_back(sampleMaxNDr[0]);
  }
  // FIXME: why next line and not just genot_out.push_back(Genotypes[i]);
  // if keepEvery > 0? We do that already.
  // It is just ugly to get a 0 in that first genotype when keepEvery < 0
  // uniqueGenotypes.insert(Genotypes[0].to_ullong());
  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;

  
#ifdef DEBUGV
  Rcpp::Rcout << "\n the initial species\n";
  print_spP(popParams[0]);
#endif


  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 4 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  
  while(!simulsDone) {
    // Check how we are doing with time as first thing.
    runningWallTime = difftime(time(NULL), start_time);
    if( runningWallTime > maxWallTime ) {
      hittedWallTime = true;
      forceSample = true;
      simulsDone = true;
    }
    
    iter++;
    if(verbosity > 1) {
      if(! (iter % iterL) ) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	Rcpp::Rcout << "\n\n    ... iteration " << iter;
	Rcpp::Rcout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      Rcpp::Rcout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    Rcpp::Rcout << " DEBUGV\n";
    Rcpp::Rcout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    if(iter == 1) { // handle special case of first iter
      // // FIXME debug
      // Rcpp::Rcout << " popSize[0] ";
      // print_spP(popParams[0]);
      // // end debug
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample, 
					 ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      //popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(to_update == 1) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, 
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n For to_update = \n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		  << "\n       genotype =  " << Genotypes[u_1] 
		  << "\n       popSize = " << popParams[u_1].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth;
#endif

      } else if(to_update == 2) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	tmpdouble2 = ti_nextTime_tmax_2_st(popParams[u_2],
					   currentTime,
					   tSample, ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_2, tmpdouble2);
	popParams[u_1].timeLastUpdate = currentTime;
	popParams[u_2].timeLastUpdate = currentTime;

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		  << "\n       genotype =  " << Genotypes[u_1] 
		  << "\n       popSize = " << popParams[u_1].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth


		  << "\n\n**     Species  = " << u_2 
		  << "\n       genotype =  " << Genotypes[u_2] 
		  << "\n       popSize = " << popParams[u_2].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble2
		  << " \n     species R " << popParams[u_2].R
		  << " \n     species W " << popParams[u_2].W
		  << " \n     species death " << popParams[u_2].death
		  << " \n     species birth " << popParams[u_2].birth;

#endif

      } else { // we sampled, so update all
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;
	  
#ifdef DEBUGV
	  Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
	    // << "\n       popParams[i].nextMutationTime = " 
	    // << popParams[i].nextMutationTime
		    << " \n     species R " << popParams[i].R
		    << " \n     species W " << popParams[i].W
		    << " \n     species death " << popParams[i].death
		    << " \n     species birth " << popParams[i].birth;
#endif
	}
      }
    }
    if(forceSample) {
      // A VERY ugly hack. Resetting tSample to jump to sampling.
      tSample = currentTime;
      // Need this, o.w. would skip a sampling.
      timeNextPopSample = currentTime;
    } 

    
    // ******************** 5.3 and do we sample? *********** 
    // Find minimum to know if we need to sample the whole pop
    getMinNextMutationTime4(nextMutant, minNextMutationTime, 
			    mapTimes);
    
    if(verbosity >= 2) {
      Rcpp::Rcout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample 
		<< "; popParams.size() = " << popParams.size() << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);



      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
#ifdef DEBUGV
	//if(verbosity > -2) {
	// We always warn about this, since interaction with ti==0
	Rcpp::Rcout << "\n Forced sampling triggered for next loop: \n    " << 
	  " popParams[nextMutant].popSize = " << 
	  popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	Rcpp::Rcout << " when nextMutant = " << nextMutant <<
	  " at iteration " << iter << "\n";
	//}
#endif
      }      
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet 
      // know total number of different species.
      // This is a protection against things going wild. Should
      // not happen in regular usage.
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2;
#ifdef DEBUGV 
      	//if(verbosity > -2) // we always warn about this
 
	Rcpp::Rcout << "\n Forced sampling triggered for next loop "
		  << " when numSpecies = " << 
	  numSpecies << " at iteration " << iter << "\n";
#endif
      }
      // Why are these lines here instead of somewhere else?
      // Right before the if for sampling or not?
      // FIXME
      // runningWallTime = difftime(time(NULL), start_time);
      // if( runningWallTime > maxWallTime ) {
      // 	hittedWalllTime = true;
      // 	forceSample = true;
      // 	simulsDone = true;
      // }


      if(popParams[nextMutant].numMutablePos != 0) {
	// this is the usual case. The alternative is the dummy or null mutation

      
	// ************   5.5   ***************
	getMutatedPos_bitset(mutatedPos, numMutablePosParent, // r,
			     ran_generator,
			     mutablePos,
			     Genotypes[nextMutant], 
			     numGenes);
      
	// ************   5.6   ***************

	newGenotype = Genotypes[nextMutant];
	newGenotype.set(mutatedPos);
	// newGenotype[mutatedPos] = 1;
      
	new_sp_bitset(sp, newGenotype, Genotypes);

	if(sp == numSpecies) {// New species
	  ++numSpecies;
	  init_tmpP(tmpParam);

	  if(verbosity >= 2) {
	    Rcpp::Rcout <<"\n     Creating new species   " << (numSpecies - 1)
			<< "         from species "  <<   nextMutant;
	  }
	
	  tmpParam.popSize = 1;

	  fitness(tmpParam, popParams[nextMutant], mutatedPos, 
		  restrictTable,
		  typeCBN, newGenotype, // birthRate,
		  s,
		  numDrivers, typeFitness, // genTime, adjust_fitness_B,
		  sh); //, adjust_fitness_MF);
	

	  if(tmpParam.birth > 0.0) {
	    tmpParam.numMutablePos = numMutablePosParent - 1;
	    if(mutationPropGrowth)
	      tmpParam.mutation = mu * tmpParam.birth * tmpParam.numMutablePos;
	    //	    tmpParam.mutation = mu * tmpParam.birth * (numMutablePosParent - 1);
	    else
	      tmpParam.mutation = mu * tmpParam.numMutablePos;
	    //tmpParam.mutation = mu * (numMutablePosParent - 1);
	    if (tmpParam.mutation > 1 )
	      Rcpp::Rcout << "WARNING: mutation > 1\n";
	    if (numMutablePosParent == 1) {
	      Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	      tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
	    }
	    W_f_st(tmpParam);
	    R_f_st(tmpParam);
	    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.
	    popParams.push_back(tmpParam);
	    Genotypes.push_back(newGenotype);
	    to_update = 2;
#ifdef MIN_RATIO_MUTS
	    g_tmp1 = tmpParam.birth/tmpParam.mutation;
	    if(g_tmp1 < g_min_birth_mut_ratio) g_min_birth_mut_ratio = g_tmp1;
	  
	    g_tmp1 = tmpParam.death/tmpParam.mutation;
	    if(g_tmp1 < g_min_death_mut_ratio) g_min_death_mut_ratio = g_tmp1;	
#endif	  
	  } else {// fitness is 0, so we do not add it
	    --sp;
	    --numSpecies;
	    to_update = 1;
	  }
	  //#ifdef DEBUGV	
	  if(verbosity >= 3) {
	    Rcpp::Rcout << " \n\n\n Looking at NEW species " << sp << " at creation";
	    Rcpp::Rcout << "\n Genotype = " << newGenotype; //Genotypes[sp];
	    Rcpp::Rcout << "\n sp_id = " << newGenotype.to_ullong() ;
	    //Genotypes[sp].to_ullong();
	    Rcpp::Rcout << "\n birth of sp = " << tmpParam.birth;
	    Rcpp::Rcout << "\n death of sp = " << tmpParam.death;
	    Rcpp::Rcout << "\n s = " << s;
	    Rcpp::Rcout << "\n parent birth = " << popParams[nextMutant].birth;
	    Rcpp::Rcout << "\n parent death = " << popParams[nextMutant].death;
	    Rcpp::Rcout << "\n parent Genotype = " << Genotypes[nextMutant];
	    print_spP(tmpParam);
	  }
	  //#endif
	} else {	// A mutation to pre-existing species
#ifdef DEBUGW
	  if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0) {
	    DP2(currentTime);
	    DP2(sp);
	    DP2(popParams[sp].timeLastUpdate);
	    print_spP(popParams[sp]);
	    throw std::out_of_range("currentTime - timeLastUpdate out of range. Serious bug");
	  }
#endif
	
	  if(verbosity >= 2) {
	    Rcpp::Rcout <<"\n     Mutated to existing species " << sp 
			<< " (Genotype = " << Genotypes[sp] 
			<< "; sp_id = " << Genotypes[sp].to_ullong() << ")"
			<< "\n from species "  <<   nextMutant
			<< " (Genotypes = " << Genotypes[nextMutant] 
			<< "; sp_id = " << Genotypes[sp].to_ullong() << ")";
	  }

	  // FIXME00: the if can be removed??
	  if(popParams[sp].popSize > 0.0) {
	    popParams[sp].popSize = 1.0 + 
	      Algo2_st(popParams[sp], currentTime, mutationPropGrowth);
	    if(verbosity >= 2) {
	      Rcpp::Rcout << "\n New popSize = " << popParams[sp].popSize << "\n";
	    }
	  } else {
	    throw std::range_error("\n popSize == 0 but existing? \n");
	  }
	
#ifdef DEBUGW
	   // This is wrong!!! if we set it to -999999, then the time to
  	  // next mutation will not be properly updated.  In fact, the
  	  // mapTimes map becomes a mess because the former pv in the
  	  // popParams is not removed so we end up inserting another pair
  	  // for the same species.

	  // popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	  //popParams[sp].Flag = true;
	}
	//   ***************  5.7 ***************
	// u_2 irrelevant if to_update = 1;
	u_1 = nextMutant;
	u_2 = static_cast<int>(sp);
      } else { // the null or dummy mutation case
	// Rcpp::Rcout << "\n null mutation; before popSize" << std::endl;
	// DP2(popParams[nextMutant].popSize);
	++popParams[nextMutant].popSize;
	to_update = 1;
	u_1 = nextMutant;
	u_2 = -99;
	// FIXME: do this conditionally on flag
	Rcpp::Rcout << "Note: updating in null mutation\n";
	// Rcpp::Rcout << "\n null mutation; after popSize" << std::endl;
	// DP2(popParams[nextMutant].popSize);
	// Rcpp::Rcout << "\n done null mutation; after popSize ********" << std::endl;
      }
    }
      else { //       *********** We are sampling **********
      to_update = 3; //short_update = false;
      if(verbosity >= 2) {
	Rcpp::Rcout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  Rcpp::Rcout << " at time " << tSample << "\n";
	} else
	  Rcpp::Rcout <<". We reached finalTime " << finalTime << "\n";
      }

      currentTime = tSample;
      if(verbosity >= 3)
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() << "\n";

      sample_all_pop_P(sp_to_remove, 
		       popParams, Genotypes, tSample, mutationPropGrowth);
      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove.size())
	remove_zero_sp_v7(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();
      
      totPopSize_and_fill_out_crude_P(outNS_i, totPopSize, 
				      lastStoredSample,
				      genot_out, 
				      //sp_id_out,
				      popSizes_out, index_out,
				      time_out, 
				      sampleTotPopSize,sampleLargestPopSize,
				      sampleMaxNDr, sampleNDrLargestPop,
				      simulsDone,
				      reachDetection,
				      lastMaxDr,
				      done_at,
				      Genotypes, popParams, 
				      currentTime,
				      numDrivers,
				      keepEvery,
				      detectionSize,
				      finalTime,
				      //endTimeEvery,
				      detectionDrivers,
				      verbosity,
				      minDetectDrvCloneSz,
				      extraTime); //keepEvery is for thinning
      if(verbosity >= 3) {
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() 
		  << "\n totPopSize after sampling " << totPopSize << "\n";
      }
      
      computeMcFarlandError(e1, n_0, n_1, tps_0, tps_1, 
			    typeFitness, totPopSize, K); //, initSize);

      // Largest error in McFarlands' method
      // if( (typeFitness == "mcfarland0") ||
      // 	  (typeFitness == "mcfarland") || 
      // 	  (typeFitness == "mcfarlandlog") ) {
      // 	tps_1 = totPopSize;
      // 	if(typeFitness == "mcfarland")
      // 	  etmp = abs( tps_1 - (tps_0 + 1) );
      // 	else {
      // 	  if( (tps_0 + 1.0) > tps_1 ) 
      // 	    etmp = (K + tps_0 + 1.0)/(K + tps_1);
      // 	  else
      // 	    etmp = (K + tps_1)/(K + tps_0 + 1);
      // 	}
      // 	if(etmp > e1) {
      // 	  e1 = etmp;
      // 	  n_0 = tps_0;
      // 	  n_1 = tps_1;
      // 	}
      // 	tps_0 = tps_1;
      // }

      // It goes here: zz: not detectionSize,
      // but the keepEvery? or sampleUntilKeep. Yes, use that.
      // endingSampleEvery.

      // Use driver criterion here!!! 
      // if endingSampleEvery
      // if totPopSize >= detectionSize: 
      //        do not break unless 
      //        tSample %% endingSampleEvery 

      // All of this has to be in totPopSize_and_fill


      // this if not sampleUntilKeep
      // if( (totPopSize >= detectionSize) ||
      // 	  (totPopSize <= 0.0) || (tSample >= finalTime)) {	
      // 	simulsDone = true;
      // 	break; // skip last update if beerenwinkel
      // }       
      
      if(simulsDone)
	break; //skip last updateRates

      // if( (typeFitness == "beerenwinkel") ) {
      // 	updateRatesBeeren(popParams, adjust_fitness_B,
      // 			  initSize, currentTime, alpha, totPopSize,
      // 			  mutationPropGrowth, mu);
      // } else if( (typeFitness == "mcfarland0") ) {
      // 	updateRatesMcFarland0(popParams, adjust_fitness_MF,
      // 			     K, totPopSize,
      // 			     mutationPropGrowth, mu);
      // } else if( (typeFitness == "mcfarland") ) {
      // 	updateRatesMcFarland(popParams, adjust_fitness_MF,
      // 			     K, totPopSize);
      // } else if( (typeFitness == "mcfarlandlog") ) {
      if( (typeFitness == "mcfarlandlog") ) {	
	updateRatesMcFarlandLog(popParams, adjust_fitness_MF,
			     K, totPopSize);
      }
      
#ifdef MIN_RATIO_MUTS
      // could go inside sample_all_pop but here we are sure death, etc, current
      // But I catch them when they are created. Is this really needed?
      for(size_t i = 0; i < popParams.size(); i++) {
	g_tmp1 = popParams[i].birth/popParams[i].mutation;
	if(g_tmp1 < g_min_birth_mut_ratio) g_min_birth_mut_ratio = g_tmp1;
	
	g_tmp1 = popParams[i].death/popParams[i].mutation;
	if(g_tmp1 < g_min_death_mut_ratio) g_min_death_mut_ratio = g_tmp1;
      }
#endif
      
      forceSample = false;
    }
  }
}


// [[Rcpp::export]]
Rcpp::List BNB_Algo5(Rcpp::IntegerMatrix restrictTable,
		     int numDrivers,
		     int numGenes,
		     Rcpp::CharacterVector typeCBN_,
		     double s, 
		     double death,
		     double mu,
		     double initSize,
		     double sampleEvery,
		     double detectionSize,
		     double finalTime,
		     int initSp,
		     int initIt,
		     int seed,
		     int verbosity,
		     int speciesFS,
		     double ratioForce,
		     Rcpp::CharacterVector typeFitness_,
		     int maxram,
		     int mutationPropGrowth,
		     int initMutant,
		     double maxWallTime,
		     double keepEvery,
		     double sh,
		     double K,
		     int detectionDrivers,
		     bool onlyCancer,
		     bool errorHitWallTime,
		     int maxNumTries,
		     bool errorHitMaxTries,
		     double minDetectDrvCloneSz,
		     double extraTime
	       ) {
  //BEGIN_RCPP
  // using namespace Rcpp;
  precissionLoss();
  const std::string typeFitness = Rcpp::as<std::string>(typeFitness_); // no need to do [0]
  const std::string typeCBN = Rcpp::as<std::string>(typeCBN_); // no need to do [0]
  // const double genTime = 4.0; // should be a parameter. For Bozic only.
  
  // const IntegerMatrix restrictTable(restrictTable_);
  // const int numDrivers = as<int>(numDrivers_);
  // const int numGenes = as<int>(numGenes_);
  // const std::string typeCBN = as<std::string>(typeCBN_);
  // const std::string typeFitness = as<std::string>(typeFitness_);
  // // birth and death are irrelevant with Bozic
  // const double birthRate = as<double>(birthRate_);
  // const double death = as<double>(death_);
  // const double s = as<double>(s_);
  // const double mu = as<double>(mu_);
  // const double initSize = as<double>(initSize_);
  // const double sampleEvery = as<double>(sampleEvery_);
  // const double detectionSize = as<double>(detectionSize_);
  // const double finalTime = as<double>(finalTime_);
  // const int initSp = as<int>(initSize_species_);
  // const int initIt = as<int>(initSize_iter_); // FIXME: this is a misnomer
  // const int verbosity = as<int>(verbose_);
  // // const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  // double ratioForce = as<double>(ratioForce_); // If a single species this times
  // // detectionSize, force a sampling to prevent going too far.
  // int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // // species
  // const int seed = as<int>(seed_gsl_);
  // const long maxram = as<int>(maxram_);
  // const int mutationPropGrowth = as<int>(mutationPropGrowth_);
  // const int initMutant = as<int>(initMutant_);
  // const double maxWallTime = as<double>(maxWallTime_);
  // const double keepEvery = as<double>(keepEvery_);

  
  // const double alpha = as<double>(alpha_);
  // const double sh = as<double>(sh_); // coeff for fitness
  // // if a driver without dependencies. Like in Datta et al., 2013.
  // const double K = as<double>(K_); //for McFarland
  // //const double endTimeEvery = as<double>(endTimeEvery_); 
  // const int detectionDrivers = as<int>(detectionDrivers_); 
  
  // // const bool errorFinalTime = as<bool>(errorFinalTime_);
  // const bool errorHitWallTime = as<bool>(errorHitWallTime_);
  // const bool onlyCancer = as<bool>(onlyCancer_);
  // const int maxNumTries = as<int>(maxNumTries_);
  // const bool errorHitMaxTries = as<bool>(errorHitMaxTries_);
  // const double minDetectDrvCloneSz = as<double>(minDetectDrvCloneSz_);
  // const double extraTime = as<double>(extraTime_);
  
  // C++11 random number
  std::mt19937 ran_generator(seed);

  // some checks. Do this systematically
  // FIXME: do only if mcfarland!
  if(K < 1 )
    throw std::range_error("K < 1.");
  // verify we are OK with usigned long long
  if( !(static_cast<double>(std::numeric_limits<unsigned long long>::max()) 
  	>= pow(2, 64)) )
    throw std::range_error("The size of unsigned long long is too short.");
  if(numGenes > 64)  
    throw std::range_error("This version only accepts up to 64 genes. Caught in R");

  bool runAgain = true;
  bool reachDetection = false;
  //Output
  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!
  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  double totPopSize = 0;
  std::vector<double> sampleTotPopSize;
  std::vector<double> sampleLargestPopSize;
  std::vector<int> sampleMaxNDr; //The number of drivers in the population
  // with the largest number of drivers; and this for each time sample
  std::vector<int> sampleNDrLargestPop; //Number of drivers in population
  // with largest size (at each time sample)
  sampleTotPopSize.reserve(initIt);
  sampleLargestPopSize.reserve(initIt);
  sampleMaxNDr.reserve(initIt);
  sampleNDrLargestPop.reserve(initIt);


  int outNS_i = -1; // the column in the outNS
  // time limits
  // FIXME think later FIXME
  time_t start_time = time(NULL);
  double runningWallTime = 0;
  bool  hittedWallTime = false;
  bool hittedMaxTries = false;
 
  // spParamsP tmpParam; 
  // std::vector<spParamsP> popParams(1);
  // const int sp_per_period = 5000;

  // popParams.reserve(sp_per_period);
  // Genotypes.reserve(sp_per_period);

  // std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset


  // // multimap to hold nextMutationTime
  // std::multimap<double, int> mapTimes;
  // //std::multimap<double, int>::iterator m1pos;


  // // count troublesome tis
  // int ti_dbl_min = 0;
  // int ti_e3 = 0;


  
  // // Beerenwinkel
  // double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  // //McFarland
  // double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  double e1, n_0, n_1; // for McFarland error
  // double tps_0, tps_1; // for McFarland error
  // tps_0 = 0.0;
  // tps_1 = 0.0;
  e1 = 0.0;
  n_0 = 0.0;
  n_1 = 0.0;

  // // For totPopSize_and_fill and bailing out
  // // should be static vars inside funct,
  // // but they keep value over calls in same R session.
  // int lastMaxDr = 0;
  // double done_at = -9;
  // // totalPopSize at time t, at t-1 and the max error.

  // 5.1 Initialize 

  int numRuns = 0;
  bool forceRerun = false;
  
  double currentTime = 0;
  int iter = 0;
  while(runAgain) {

    // Initialize a bunch of things
    
// #ifdef MIN_RATIO_MUTS
//   g_min_birth_mut_ratio = DBL_MAX;
//   g_min_death_mut_ratio = DBL_MAX;
//   g_tmp1 = DBL_MAX;
// #endif

  
  // untilcancer goes here
  

  
  
  //tmpParam is a temporary holder. 
  // init_tmpP(tmpParam);
  // init_tmpP(popParams[0]);

  // lastStoredSample = 0.0;
  // Genotypes[0].reset();
  // popParams[0].popSize = initSize;
  // totPopSize = initSize;

  // tps_0 = totPopSize;
  // e1 = 0.0;
  // tps_1 = totPopSize;



    try {
      // it is CRUCIAL that several entries are zeroed (or -1) at the
      // start of innerBNB now that we do multiple runs if onlyCancer = true.
      innerBNB(
	       numGenes,
	       initSize,
	       K,
	       // alpha,
	       typeCBN,
	       // genTime,
	       typeFitness,
	       mutationPropGrowth,
	       mu,
	       sh,
	       s,
	       death,
	       // birthRate,
	       keepEvery,
	       sampleEvery,		     
	       numDrivers,
	       initMutant,
	       start_time,
	       maxWallTime,
	       finalTime,
	       detectionSize,
	       //endTimeEvery,
	       detectionDrivers,
	       minDetectDrvCloneSz,
	       extraTime,
	       verbosity,
	       totPopSize,
	       e1,
	       n_0,
	       n_1,
	       ratioForce,
	       currentTime,
	       speciesFS,
	       outNS_i, 
	       iter,
	       genot_out,
	       popSizes_out,
	       index_out,
	       time_out,
	       sampleTotPopSize,
	       sampleLargestPopSize,
	       sampleMaxNDr,
	       sampleNDrLargestPop,
	       reachDetection,
	       ran_generator,
	       runningWallTime,
	       hittedWallTime,
	       restrictTable);
      ++numRuns;
      forceRerun = false;
    } catch (rerunExcept &e) {
      Rcpp::Rcout << "\n Exception " << e.what() 
		  << ". Rerunning.";
      forceRerun = true;
    } catch (const std::exception &e) {
      Rcpp::Rcout << "\n Unrecoverable exception: " << e.what()
		  << ". Aborting. \n";
      return
	Rcpp::List::create(Named("other") =
		     Rcpp::List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = e.what()));
    } catch (...) {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting. \n";
      return
	Rcpp::List::create(Named("other") =
		     Rcpp::List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = "Unknown exception"));
    }


    if(hittedWallTime) {
      Rcpp::Rcout << "\n Hitted wall time. Exiting.";
      runAgain = false;
      if(errorHitWallTime) {
	Rcpp::Rcout << "\n Hitting wall time is regarded as an error. \n";
	return
	  Rcpp::List::create(Named("HittedWallTime") = true,
		       Named("HittedMaxTries") = false, // yes, for
							// coherent return
							// objects
		       Named("other") =
			     Rcpp::List::create(Named("UnrecoverExcept") = false));
      }
    } else if(numRuns > maxNumTries) {
      //  hittedMaxTries
      hittedMaxTries = true;
      Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
      runAgain = false;
      if(errorHitMaxTries) {
	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
	return
	  Rcpp::List::create(Named("HittedWallTime") = false,
		       Named("HittedMaxTries") = true,
		       Named("other") =
		       Rcpp::List::create(Named("UnrecoverExcept") = false));
      }
    } else if(forceRerun) {
      runAgain = true;
      forceRerun = false;
    } else {
      if(onlyCancer) {
	runAgain = !reachDetection;
      } else {
	runAgain = false;
      }
    }
#ifdef DEBUGV
      Rcpp::Rcout << "\n reachDetection = " << reachDetection;
      Rcpp::Rcout << "\n forceRerun =  " << forceRerun  << "\n";
      
#endif
    
  } // runAgain loop
  // FIXME: zz
  // untilcancer
  // inner loop ends above
  // The return objects only created if needed

  
  // If we hit wallTime, we can get done without going through
  // totPopSize.... Problem if sampling at end
  // if ( hittedWallTime ) {
  //   // hitted wall time. So we need to sample at the very end.
  // Nope! Just ensure if hittedWallTime you always sample properly!
  // }
    

  
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species

  // here("right after simuls done");

  // FIXME: all this is ugly and could be a single function
  // up to call to IntegerMatrix
  std::set<unsigned long long> uniqueGenotypes;
  std::vector<unsigned long long> genot_out_ullong(genot_out.size());
  genot_out_to_ullong(genot_out_ullong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ullong);
  std::vector<unsigned long long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);
  // IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);
  Rcpp::IntegerMatrix returnGenotypes(numGenes, uniqueGenotypes_vector.size());
  // here("after creating returnGenotypes");
  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);
  // here("after call to create_returnGenotypes_to_vector");

  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)
  // here("after uniqueGenotypes_to_vector");
 

  int outNS_r, outNS_c, create_outNS;
  if( ( (uniqueGenotypes.size() + 1) *  (outNS_i + 1) ) > ( pow(2, 31) - 1 ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > 2^31 - 1. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else if ( 
	     static_cast<long>((uniqueGenotypes.size()+1) * (outNS_i+1)) * 8 > 
	     (maxram * (1024*1024) ) ) {
    Rcpp::Rcout << "\nWARNING: Return outNS object > maxram. Not created.\n";
    outNS_r = 1;
    outNS_c = 1;
    create_outNS = 0;
  } else {
    outNS_r = outNS_i + 1;
    outNS_c = uniqueGenotypes.size() + 1;
    create_outNS = 1;
  }
  Rcpp::NumericMatrix outNS(outNS_r, outNS_c);  
  if(create_outNS) {
    reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ullong, 
		     popSizes_out, 
		     index_out, time_out);
    
  } else {
    outNS(0, 0) = -99;
  }

  int maxNumDrivers = 0;
  int totalPresentDrivers = 0;
  std::vector<int>countByDriver(numDrivers, 0);
  std::string occurringDrivers;
  

  // here("before count_NumDrivers");
  // IntegerVector totDrivers(returnGenotypes.ncol());
  // count_NumDrivers(maxNumDrivers, returnGenotypes, numDrivers,
  // 		   totDrivers);
  count_NumDrivers(maxNumDrivers, countByDriver,
		   returnGenotypes, numDrivers);

  whichDrivers(totalPresentDrivers, occurringDrivers, countByDriver);

  std::vector<double> sampleLargestPopProp(outNS_i + 1);

  if((outNS_i + 1) != static_cast<int>(sampleLargestPopSize.size()))
    throw std::length_error("outNS_i + 1 != sampleLargestPopSize.size");
  std::transform(sampleLargestPopSize.begin(), sampleLargestPopSize.end(),
		 sampleTotPopSize.begin(),
		 sampleLargestPopProp.begin(),
		 std::divides<double>());

  Rcpp::NumericMatrix perSampleStats(outNS_i + 1, 5);
  fill_SStats(perSampleStats, sampleTotPopSize, sampleLargestPopSize,
	      sampleLargestPopProp, sampleMaxNDr, sampleNDrLargestPop);
  
  // error in mcfarland's
  // if((typeFitness == "mcfarland0") || (typeFitness == "mcfarlandlog"))
  //   e1r = log(e1);
  // if(typeFitness == "mcfarland")
  //   e1r = (1.0/K) * e1;

  // here("before return");

  // // // debuggin: precompute things
  // DP2(simulsDone);
  // DP2(maxWallTime);
  // DP2(hittedWallTime);
  // DP2(outNS_i);
  // DP2( sampleMaxNDr[outNS_i]);
  // DP2(sampleNDrLargestPop[outNS_i]);
  // DP2(sampleLargestPopSize[outNS_i]);
  // DP2(sampleLargestPopProp[outNS_i]);
  // DP2((runningWallTime > maxWallTime));
  // here("after precomp");
  // here("*******************************************");


  return 
    Rcpp::List::create(Named("pops.by.time") = outNS,
		 Named("NumClones") = uniqueGenotypes.size(), 
		 Named("TotalPopSize") = totPopSize,
		 Named("Genotypes") = returnGenotypes,
		 Named("MaxNumDrivers") = maxNumDrivers,
		 // Named("MaxDrivers_PerSample") = wrap(sampleMaxNDr),
		 // Named("NumDriversLargestPop_PerSample") = sampleNDrLargestPop,
		 // Named("TotPopSize_PerSample") = sampleTotPopSize,
		 // Named("LargestPopSize_PerSample") = sampleLargestPopSize,
		 // Named("PropLargestPopSize_PerSample") = sampleLargestPopProp,
		 Named("MaxDriversLast") = sampleMaxNDr[outNS_i],
		 Named("NumDriversLargestPop") =  sampleNDrLargestPop[outNS_i],
		 Named("LargestClone") = sampleLargestPopSize[outNS_i],
		 Named("PropLargestPopLast") = sampleLargestPopProp[outNS_i],
		 // Named("totDrivers") = totDrivers,
		 Named("FinalTime") = currentTime,
		 Named("NumIter") = iter,
		 //		 Named("outi") = outNS_i + 1, // silly. Use the real number of samples. FIXME
		 Named("HittedWallTime") = hittedWallTime, // (runningWallTime > maxWallTime),
		 Named("HittedMaxTries") = hittedMaxTries,
		 // Named("iRunningWallTime") = runningWallTime,
		 // Named("oRunningWallTime") = difftime(time(NULL), start_time),
		 // Named("ti_dbl_min") = ti_dbl_min,
		 // Named("ti_e3") = ti_e3,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 // FIXME: OccurringDrivers underestimates true occurring
		 // drivers if keepEvery < 0, so we only return the last.
		 Named("OccurringDrivers") = occurringDrivers,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = Rcpp::List::create(Named("attemptsUsed") = numRuns,
					       Named("errorMF") = 
						     returnMFE(e1, // K, 
							 typeFitness),
					       Named("errorMF_size") = e1,
					       Named("errorMF_n_0") = n_0,
#ifdef MIN_RATIO_MUTS
					       Named("minDMratio") =
					       g_min_death_mut_ratio,
					       Named("minBMratio") =
					       g_min_birth_mut_ratio,      
#else
					       Named("minDMratio") = -99,
					       Named("minBMratio") = -99,
#endif
					       Named("errorMF_n_1") = n_1,
					       Named("UnrecoverExcept") = false)
		 );

  // END_RCPP
    
    }

// FIXME: it would have been nice to have OccuringDrivers for last sampling

// Count by driver: count the number of time each driver has appeared. Simply count
//                  over genotypes, so a very crude measure.

// OccurringDrivers: which drivers have ever appeared (might have disappeared)

// TotalPresentDrivers: how many drivers have ever appeared. If a driver appears once, 
//                   at one sampling period, it is counted. How is this different
//                   from MaxNumDrivers? This is the individual drivers.

// MaxNumDrivers:  what is the largest number of drivers in any genotype, ever.
//                 Note that that genotype might have disappeared. How is this different
//                 from TotalPresentDrivers? This is the max drivers in a genotype!

// MaxDriversLast: the largest number of drivers in any of the genotypes
//                 that exist at the last sample.

// NumDriversLargestPopLast: the number of drivers in the population with
//                 largest pop. size at the last sample.

// The above are also parte of the matrix "PerSampleStats"

// PerSampleStats: five columns, with one entry per sample
//                 (and you get the time from pops.by.time, the first column)
//                 - Total Population Size
//                 - Pop size of the population with largest size
//                 - Ratio of the previous two
//                 - the largest number of drivers in any of the genotypes
//                 that exist in that sample (i.e., the per-sample equivalent
//                 of MaxDriversLast or MaxNumDrivers)
//                 - the number of drivers in the population with
//                 largest pop. size (i.e., the per-sample equivalent of
//                 NumDriversLargestPopLast)




// ********************************************************
// ********************************************************

// timeLastUpdate, sampling, etc.  

//timeLastUpdate only used on the right (i.e. its value is used for
// something in):
//  5.4 (call to Algo3, update the mutated)
//  5.6 (call to Algo2, update a previous sp. to which it mutates)
//  5.9 (sampling)

// timeLastUpdate ONLY updated on 5.2, when flagged ones are updated.
// (I add a few other places, to catch errors)

// Any species which are modified (5.4 Algo3 in step 5.7, 5.6 via Algo2,
// or created new) will be flagged

// When we sample, all have been unflagged as we unflag all flagged in 5.2


// My former confusion might have come from trying to sample by
// iteration, not time.

// ********************************************************
// ********************************************************



//  ************** Forced sampling idea ************************

// I was doing it wrong. I was forcing sampling, and jumpint to next
// sampling time, but possibly a minNextMutationTime was smaller!
// Now I do it right, but requires ugly hacks.
// Note initial part of 5.2: tSample is set to currentTime
// This is done after a previous asignment to tSample. I force
// getting mutations time, and proceeding pretending we are 
// doing the usual thing, but then set tSample to the previous
// currentTime, which unless ti = 0, will make us enter sampling.
// An ugly hack anyway. Just a way of forcing to go to sampling.

// Alternatively, we could set tSample to currentTime, and not go
// through updateint nextMutation or finding the minNextMutationTime
// This would be slightly faster. But breaks logic of algorithm and
// forces using a &&!forceSample in the test to not sample.

// Current approach would also allow to decreasing sampling time
// at time progresses, respecting logic of algorithm.

// Forced sampling, however, unlikely to avoid getting stuck except
// as a way of early sampling for large pop. Remember
// than at the sampling, at least two pops will have a size >= 1.
//
// Note: I do not got through Algo2 in 5.7, because I check t == 0.0.
//  But warning if enter Algo2 with t == 0 anyway.



// ALTERNATIVES FOR 5.6
// The following are two alternative approaches
// which do one fewer comparison. Harder to understand, though{
// V1
//   s = 0;
//   k = 0;

//   for( ; ; ){
// 	if( k == numGenes) { //pre-existing species
// 	  break;
// 	}
// 	else {
// 	  if(newGenotype[k] == Genotypes(k, s)) k++;
// 	  else {
// 	    s++;
// 	    if(numSpecies >= s) { //create new species
// 	      break;
// 	    } else {
// 	      k = 0;
// 	    }
// 	  }
// 	}
//   }
// V2
//   s = 0;
//   for(;;) {
// 	k = 0;
// 	while ( newGenotype[k] == Genotypes(k, s) ) {
// 	  k++;
// 	  if(k == numGenes ) {
// 	    // s is same as new genotype
// 	    // Thus: a mutation to a pre-existing species
// 	    // break
// 	  }
// 	}
// 	s++;
// 	if( s == numSpecies ) {
// 	  // create new species
// 	  // break
// 	}
//   }   
// }




// suppose this tree
// 1 -> 2
// 5 -> 4 -> 3
// 6

// grandparente has 1 and 2
//  father has 1, 2, 3: thus gets same fitness as grandpa
// child has 1, 2, 3, 4: has 4 has unment, it gets same fitness as father, and thus as 
//   grandpa, but it should get the fitness of three drivers, since
// 4 meets requirement of 3.
// Now grandchild gets 5. here things get fixed.

// And the same if unment.
// Suppose grandparen gets 1, 2, 3.
// If father now gets 6, it counts as four drivers, but it really only has
// three with the restrictions met.
