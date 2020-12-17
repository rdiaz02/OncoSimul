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


// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
// FIXMEmaybe: include randutils. Should work now but will upset many
// tests with fixed seeds.
#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"
#include "new_restrict.h"
#include <cfloat>
#include <limits>
#include <Rcpp.h>
#include <iostream>
#include <random>
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

// To track if mutation is really much smaller than birth/death Yes:
// global vars. FIXMEmaybe: turn them into non-globals. This is ugly as
// hell.
double g_min_birth_mut_ratio_nr = DBL_MAX;
double g_min_death_mut_ratio_nr = DBL_MAX;
double g_tmp1_nr = DBL_MAX;



// Obtain birth and death for a new genotype.
//     return values that signal no-viability when appropriate
//                    (to avoid adding to tables)
//     parentP used to provide the model-specific details
void nr_fitness(spParamsP& tmpP,
		const spParamsP& parentP,
		const Genotype& ge,
		const fitnessEffectsAll& F,
		const TypeModel typeModel,
		std::vector<Genotype>& Genotypes,
		std::vector<spParamsP>& popParams,
		const double& currentTime) {

  // We want a way to signal immediate non-viability of a clone. For
  // "birth-based" models that happens when any s = -1, as the fitness is
  // 0. By setting birth = 0.0 we ensure this clone does not get added and
  // we never reach into algo2, etc, leading to numerical problems.

  // With Bozic models, which are "death-based", it is different. For
  // bozic1, deaths of around 50 lead to numerical issues.  The general
  // rule is: set those mutations to -inf, so prodDeathFitness returns an
  // inf for death, and that is recognized as "no viability" (anything
  // with death > 99)

  // The ones often used are bozic1, exp, mcfarlandlog

  // For frequency dependence we need to add the candidate genotype
  // to evaluate fitness. At end we pop it out.

  // Could use const and create a copy, but possibly expensive? FIXMEmaybe:SPEED/SAFETY
  // This is FDF: those tables are most likely tiny
  if(F.frequencyDependentFitness){
    popParams.push_back(tmpP);
    Genotypes.push_back(ge);
  }

  if(typeModel == TypeModel::bozic1) {
    tmpP.death = prodDeathFitness(evalGenotypeFitness(ge, F, Genotypes, popParams, currentTime));
    if( tmpP.death > 99) {
      tmpP.birth = 0.0;
    } else {
      tmpP.birth = 1.0;
    }
  } else {
    double fitness = prodFitness(evalGenotypeFitness(ge, F, Genotypes, popParams, currentTime));
    if( fitness <= 0.0) {
      tmpP.absfitness = 0.0;
      tmpP.death = 1.0;
      tmpP.birth = 0.0;
    } else{
      // Set appropriate model-specific defaults and change only as needed
      tmpP.death = parentP.death; // will use 1 for exp, log whatever for McF
      tmpP.absfitness = parentP.absfitness; // was used for old Beerenwinkel model
      tmpP.birth = fitness;
      }
  }
  if(F.frequencyDependentFitness){
    popParams.pop_back();
    Genotypes.pop_back();
  }
}


// is this really any faster than the one below?
inline void new_sp_v(unsigned int& sp,
		     const Genotype& newGenotype,
		     const std::vector<Genotype> Genotypes) {
  sp = 0;
  for(sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] )
      break;
  }
}

// inline unsigned int new_sp(const Genotype& newGenotype,
// 		    const std::vector<Genotype> Genotypes) {
//   for(unsigned int sp = 0; sp < Genotypes.size(); ++sp) {
//     if( newGenotype == Genotypes[sp] ) {
//       return sp;
//     }
//   }
//   return Genotypes.size();
// }


void remove_zero_sp_nr(std::vector<int>& sp_to_remove,
			      std::vector<Genotype>& Genotypes,
			      std::vector<spParamsP>& popParams,
			      std::multimap<double, int>& mapTimes) {
  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype>::iterator Genotypes_begin = Genotypes.begin();
  std::vector<int>::reverse_iterator r = sp_to_remove.rbegin();
  int remove_this;
  while(r != sp_to_remove.rend() ) {
    remove_this = *r;
    mapTimes.erase(popParams[remove_this].pv);
    popParams.erase(popParams_begin + remove_this);
    Genotypes.erase(Genotypes_begin + remove_this);
    ++r;
  }
}


inline void driverCounts(int& maxNumDrivers,
			 int& totalPresentDrivers,
			 std::vector<int>& countByDriver,
			 std::vector<int>& presentDrivers,
			 Rcpp::IntegerMatrix& returnGenotypes,
			 const vector<int>& drv){
  // Fill up the "countByDriver" table, how many genotypes each driver is
  // present.  Return the maximum number of mutated drivers in any
  // genotype, the vector with just the present drivers, and the total
  // number of present drivers.

  // We used to do count_NumDrivers and then whichDrivers
  maxNumDrivers = 0;
  int tmpdr = 0;
  int driver_indx = 0; // the index in the driver table
  for(int j = 0; j < returnGenotypes.ncol(); ++j) {
    tmpdr = 0;
    driver_indx = 0;
    for(int i : drv) {
      tmpdr += returnGenotypes(i - 1, j);
      countByDriver[driver_indx] += returnGenotypes(i - 1, j);
      ++driver_indx;
    }
    if(tmpdr > maxNumDrivers) maxNumDrivers = tmpdr;
  }
  if(returnGenotypes.ncol() > 0) {
    STOPASSERT(driver_indx == static_cast<int>( countByDriver.size()));
  } else {
    STOPASSERT(driver_indx <= static_cast<int>( countByDriver.size()));
  }
  for(size_t i = 0; i < countByDriver.size(); ++i) {
    if(countByDriver[i] > 0) {
      presentDrivers.push_back(i + 1);
      ++totalPresentDrivers;
    }
  }
}



// Determine if stopping conditions reached, fill out output structures,
// return sample summaries for popsize, drivers, compute totPopSize
// Could be broken down into a couple of functions (with better names :) )
void nr_totPopSize_and_fill_out_crude_P(int& outNS_i,
					double& totPopSize,
					double& lastStoredSample,
					std::vector<Genotype>& genot_out,
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
					const std::vector<Genotype>& Genotypes,
					const std::vector<spParamsP>& popParams,
					const double& currentTime,
					const double& keepEvery,
					const double& detectionSize,
					const double& finalTime,
					// const double& endTimeEvery,
					const int& detectionDrivers,
					const int& verbosity,
					const double& minDetectDrvCloneSz,
					const double& extraTime,
					const vector<int>& drv,
					const double& cPDetect,
					const double& PDBaseline,
					const double& checkSizePEvery,
					double& nextCheckSizeP,
					std::mt19937& ran_gen,
					const bool& AND_DrvProbExit,
					const std::vector<std::vector<int> >& fixation_l,
					const double& fixation_tolerance,
					const int& min_successive_fixation,
					const double& fixation_min_size,
					int& num_successive_fixation,
					POM& pom,
					const std::map<int, std::string>& intName,
					const fitness_as_genes& genesInFitness,
					const double& fatalPopSize = 1e15
					) {
  
  // This determines if we are done or not by checking popSize, number of
  // drivers, etc
  
  // static int lastMaxDr = 0; // preserves value across calls to Algo5 from R.
  // so can not use it.
  bool storeThis = false;
  bool checkSizePNow = false;
  totPopSize = 0.0;
  
  // this could all be part of popSize_over_m_dr, with a better name
  int tmp_ndr = 0;
  int max_ndr = 0;
  double popSizeOverDDr = 0.0;

  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
    tmp_ndr = getGenotypeDrivers(Genotypes[i], drv).size();
    if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
    if(tmp_ndr >= detectionDrivers) popSizeOverDDr += popParams[i].popSize;
  }
  lastMaxDr = max_ndr;

  // Until fixation done here. Recall we use an OR operation for exiting
  // below.  Could be added to loop above.
  // And we call allGenesinGenotype also above, inside getGenotypeDrivers.
  // So room for speed ups?
  
  // Since we convert each genotype to a sorted allGenesinGenotype, iterate
  // over that first. Add that pop size if the combination is present in genotype.
  bool fixated = false;
  if(totPopSize > 0) { // Avoid silly things
    if( fixation_l.size() ) {
      std::vector<double> popSize_fixation(fixation_l.size());
      for(size_t i = 0; i < popParams.size(); ++i) {
	std::vector<int> thisg = allGenesinGenotype(Genotypes[i]);
	for(size_t fc = 0; fc != popSize_fixation.size(); ++fc) {
	  // Yes, fixation_l is sorted in R.
	  // if fixation_l[fc] starts with a -9, we are asking
	  // for exact genotype equality
	  if(fixation_l[fc][0] == -9) {
	    // // exact genotype identity?
	    std::vector<int> this_fix(fixation_l[fc].begin() + 1,
				      fixation_l[fc].end());
	    if(thisg == this_fix) {
	      popSize_fixation[fc] = popParams[i].popSize;
	    }
	  } else {
	  // gene combination in fixation element present in genotype?
	    if(std::includes(thisg.begin(), thisg.end(),
			     fixation_l[fc].begin(), fixation_l[fc].end()) ) {
	      popSize_fixation[fc] += popParams[i].popSize;
	    }
	  }
	}
      }
      // Any fixated? But avoid trivial of totPopSize of 0!
      // Now check of > 0 is redundant as we check totPopSize > 0
      // Do we want tolerance around that value?
      double max_popSize_fixation =
	*std::max_element(popSize_fixation.begin(), popSize_fixation.end());
      if( (max_popSize_fixation >= fixation_min_size ) &&
	  (max_popSize_fixation >= (totPopSize * (1 - fixation_tolerance) )) ) {
	++num_successive_fixation;
	if( num_successive_fixation >= min_successive_fixation) fixated = true;
      } else {
	num_successive_fixation = 0;
      }
    }
  }

  if (keepEvery < 0) {
    storeThis = false;
  } else if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }

  
  // Doing this is cheaper than drawing unnecessary runifs.
  // Equality, below, leads to suprises with floating point arith.

  // Operates the same as we do with keepEvery, but here we
  // compute the jump in each accepted sample. And here we use >, not
  // >=
  if(currentTime > nextCheckSizeP) {
    checkSizePNow = true;
    nextCheckSizeP = currentTime + checkSizePEvery;
    // minimal jump can be smaller than checkSizePEvery
  } else {
    checkSizePNow = false;
  }

  // We do not verify that conditions for exiting are also satisfied
  // at the exit time when extraTime > 0. We could do that,
  // checking again for the conditions (or the reasonable conditions, so
  // probably not detectSizeP). For instance, with fixated.

  // For fixated in particular, note that we evaluate fixation always, but
  // we might be exiting when there is no longer fixation.  But the logic
  // with fixation is probably to use as large a min_successive_fixation
  // as desired and no extraTime.

  // Probably would not need to check lastMaxDr and popSizeOverDDr
  // as those should never decrease. Really?? FIXME_CHECK
  
  if(AND_DrvProbExit) {
    // The AND of detectionProb and drivers
    // fixated plays no role here, and cannot be passed from R
    if(extraTime > 0) {
      if(done_at <  0) {
	if( (lastMaxDr >= detectionDrivers) &&
	    (popSizeOverDDr >= minDetectDrvCloneSz) &&
	    checkSizePNow  &&
	    detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen) ) {
	  done_at = currentTime + extraTime;
	}
      } else if (currentTime >= done_at) {
  	simulsDone = true;
  	reachDetection = true;
      }
    } else if( (lastMaxDr >= detectionDrivers) &&
	       (popSizeOverDDr >= minDetectDrvCloneSz) &&
	       checkSizePNow  &&
	       detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen) ) {
      simulsDone = true;
      reachDetection = true; 
    }
  } else {
    // The usual OR mechanism of each option
    if(extraTime > 0) {
      if(done_at <  0) {
	if( (fixated) ||
	    (totPopSize >= detectionSize) ||
	    ( (lastMaxDr >= detectionDrivers) &&
	      (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	    ( checkSizePNow  &&
	      detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen))) {
	  done_at = currentTime + extraTime;
	}
      } else if (currentTime >= done_at) {
	  simulsDone = true;
	  reachDetection = true;
      }
    } else if( (fixated) ||
	       (totPopSize >= detectionSize) ||
	       ( (lastMaxDr >= detectionDrivers) &&
		 (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	       ( checkSizePNow  &&
		 detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen)) ) {
      simulsDone = true;
      reachDetection = true;
    }
  }


  if(totPopSize >= fatalPopSize) {
    Rcpp::Rcout << "\n\totPopSize > " << fatalPopSize
		<<". You are likely to loose precision and run into numerical issues\n";
       }

  if(simulsDone)
    storeThis = true;

  // Code repeated below. But avoid looping twice if storeThis
  if( storeThis ) {
    lastStoredSample = currentTime;
    outNS_i++;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    int largest_clone = -99;

    time_out.push_back(currentTime);

    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);

      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = getGenotypeDrivers(Genotypes[i], drv).size();
	largest_clone = i;
      }
    }
    sampleTotPopSize.push_back(totPopSize);
    sampleLargestPopSize.push_back(l_pop_s);
    sampleMaxNDr.push_back(max_ndr);
    sampleNDrLargestPop.push_back(ndr_lp);

    if(l_pop_s > 0) {
      if (largest_clone < 0)
    	throw std::logic_error("largest_clone < 0");
      addToPOM(pom, Genotypes[largest_clone], intName, genesInFitness);
    } else {
      addToPOM(pom, "_EXTINCTION_");
    }
  }

  if( !std::isfinite(totPopSize) ) {
    throw std::range_error("totPopSize not finite");
  }
  if( std::isnan(totPopSize) ) {
    throw std::range_error("totPopSize is NaN");
  }
  // For POM
  if( !storeThis ) {
    double l_pop_s = 0.0;
    int largest_clone = -99;
    for(size_t i = 0; i < popParams.size(); ++i) {
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	largest_clone = i;
      }
    }
    if(l_pop_s > 0) {
      if (largest_clone < 0)
	throw std::logic_error("largest_clone < 0");
      addToPOM(pom, Genotypes[largest_clone], intName, genesInFitness);
    } else {
      addToPOM(pom, "_EXTINCTION_");
    }
  }
}


inline void nr_reshape_to_outNS(Rcpp::NumericMatrix& outNS,
				const vector<vector<int> >& uniqueGenotV,
				const vector<vector<int> >& genot_out_v,
				const vector<double>& popSizes_out,
				const vector<int>& index_out,
				const vector<double>& time_out){

  // I might want to return the actual drivers in each period
  // and the actual drivers in the population with largest popsize
  // Something like what we do now with whichDrivers
  // and count_NumDrivers?

  vector<vector<int> >::const_iterator fbeg = uniqueGenotV.begin();
  vector<vector<int> >::const_iterator fend = uniqueGenotV.end();

  int column;

  for(size_t i = 0; i < genot_out_v.size(); ++i) {
    column = std::distance(fbeg, lower_bound(fbeg, fend, genot_out_v[i]) );
    outNS(index_out[i], column + 1) =  popSizes_out[i];
  }

  for(size_t j = 0; j < time_out.size(); ++j)
    outNS(j, 0) = time_out[j];
}


Rcpp::NumericMatrix create_outNS(const vector<vector<int> >& uniqueGenotypes,
				 const vector<vector<int> >& genot_out_v,
				 const vector<double>& popSizes_out,
				 const vector<int>& index_out,
				 const vector<double>& time_out,
				 const int outNS_i, const int maxram) {
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
    nr_reshape_to_outNS(outNS, uniqueGenotypes,
			genot_out_v,
			popSizes_out,
			index_out, time_out);

  } else {
    outNS(0, 0) = -99;
  }
  return outNS;
}


vector< vector<int> > uniqueGenot_vector(vector<vector<int> >& genot_out) {
  // From genot_out we want the unique genotypes, but each as a single
  // vector. Convert to the vector, then use a set to give unique sorted
  // vector.
  std::set<std::vector<int> > uniqueGenotypes_nr(genot_out.begin(),
						 genot_out.end());
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr (uniqueGenotypes_nr.begin(),
							    uniqueGenotypes_nr.end());
  return uniqueGenotypes_vector_nr;
}



std::vector<std::vector<int> > genot_to_vectorg(const std::vector<Genotype>& go) {
  std::vector<std::vector<int> > go_l;
  std::transform(go.begin(), go.end(), back_inserter(go_l), genotypeSingleVector);
  return go_l;
}



std::string driversToNameString(const std::vector<int>& presentDrivers,
			    const std::map<int, std::string>& intName) {
  std::string strDrivers;
  std::string comma = "";
  for(auto const &g : presentDrivers) {
    strDrivers += (comma + intName.at(g));
    comma = ", ";
  }
  return strDrivers;
}


std::string genotypeToNameString(const std::vector<int>& genotypeV,
				 const fitness_as_genes& fg,
				 const std::map<int, std::string>& intName) {

  // The genotype vectors are returned as a string of names. Similar to
  // the Int version, but we map here to names.

  // As the fitness is stored in terms of modules, not genes, we need to
  // check if a _gene_ is in the order part or not by mapping back to
  // modules. That is the fitness_as_genes argument.

  std::string strGenotype;

  std::vector<int> order_int;
  std::vector<int> rest_int;

  for(auto const &g : genotypeV) {
    if( binary_search(fg.orderG.begin(), fg.orderG.end(), g)) {
      order_int.push_back(g);
    } else {
      rest_int.push_back(g);
    }
  }

  std::string order_sep = " _ ";
  std::string order_part;
  std::string rest;
  std::string comma = "";

  // FIXMEmaybe: when sure no problems, remove if needed for speed.
  for(auto const &g : order_int) {
    order_part += (comma + intName.at(g));
    comma = " > "; // comma = ", ";
  }
  comma = "";
  for(auto const &g : rest_int) {
    rest += (comma + intName.at(g));
    comma = ", ";
  }
  if(fg.orderG.size()) {
    strGenotype = order_part + order_sep + rest;
  } else {
    strGenotype = rest;
  }
  return strGenotype;
}


std::vector<std::string> genotypesToNameString(const std::vector< vector<int> >& uniqueGenotypesV,
					       const fitness_as_genes fg,
					       const std::map<int, std::string>& intName) {
  std::vector<std::string> gs;
  for(auto const &v: uniqueGenotypesV )
      gs.push_back(genotypeToNameString(v, fg, intName));
  return gs;
}


Rcpp::IntegerMatrix nr_create_returnGenotypes(const int& numGenes,
					      const std::vector< vector<int> >& uniqueGenotypesV){
  // We loose order here. Thus, there might be several identical columns.
  Rcpp::IntegerMatrix returnGenotypes(numGenes, uniqueGenotypesV.size());
  for(size_t i = 0; i < uniqueGenotypesV.size(); ++i) {
    for(int j : uniqueGenotypesV[i]) {
      returnGenotypes(j - 1, i) = 1;
    }
  }
  return returnGenotypes;
}






static void nr_sample_all_pop_P(std::vector<int>& sp_to_remove,
				std::vector<spParamsP>& popParams,
				// next only used with DEBUGV
				const std::vector<Genotype>& Genotypes,
				const double& tSample,
				const int& mutationPropGrowth){
  sp_to_remove.clear();

  for(size_t i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
    DEBUG_nr;

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
      sp_to_remove.push_back(i);
      DEBUG_nr2;
    }
  }
}

void addToPhylog(PhylogName& phylog,
		 const Genotype& parent,
		 const Genotype& child,
		 const double time,
		 const std::map<int, std::string>& intName,
		 const fitness_as_genes& fg,
		 const double pop_size_child) {
  phylog.time.push_back(time);
  phylog.parent.push_back(genotypeToNameString(genotypeSingleVector(parent),
					       fg, intName));
  phylog.child.push_back(genotypeToNameString(genotypeSingleVector(child),
					      fg, intName));
  phylog.pop_size_child.push_back(pop_size_child);
}



void addToLOD(std::map<std::string, std::string>& lod,
	      const Genotype& parent,
	      const Genotype& child,
	      const std::map<int, std::string>& intName,
	      const fitness_as_genes& fg) {
  // Only called when the child has pop size of 0
  // so true LOD
  // Use a map for LOD, and overwrite the parent:
  // we only add when the size of the child is 0
  // The key of the map is the child.

  // We might want to store the time? Not really clear even if that
  // makes sense. We would be storing the last time the child (which had 0
  // size at that time) arose from the parent.
  // A simple kludge is to have two maps, the second with child and time.
  // Or do it properly as map<int, genot_time_struct>
  // genot_time_struct {string parent; double time}

  std::string parent_str = genotypeToNameString(genotypeSingleVector(parent),
						fg, intName);
  std::string child_str = genotypeToNameString(genotypeSingleVector(child),
					       fg, intName);
  lod[child_str] = parent_str;
}


void addToPOM(POM& pom,
	      const Genotype& genotype,
	      const std::map<int, std::string>& intName,
	      const fitness_as_genes& fg) {

  if (pom.genotypes.empty()) {
    std::string g = genotypeToNameString(genotypeSingleVector(genotype),
				       fg, intName);
    pom.genotypesString.push_back(g);
    pom.genotypes.push_back(genotype);
  } else if ( !(pom.genotypes.back() == genotype) ) {
    // Insert only if different from previous
    std::string g = genotypeToNameString(genotypeSingleVector(genotype),
				       fg, intName);
    pom.genotypesString.push_back(g);
    pom.genotypes.push_back(genotype);
  }
}

// to explicitly signal extinction
void addToPOM(POM& pom,
	      const std::string string) {
  pom.genotypesString.push_back(string);
}





static void nr_innerBNB (const fitnessEffectsAll& fitnessEffects,
			 const std::vector<double>& initSize,
			 const double& K,
			 const TypeModel typeModel,
			 const int& mutationPropGrowth,
			 const std::vector<double>& mu,
			 const double& death,
			 const double& keepEvery,
			 const double& sampleEvery,
			 const std::vector<std::vector<int> >& initMutant,
			 const time_t& start_time,
			 const double& maxWallTime,
			 const double& finalTime,
			 const double& detectionSize,
			 const int& detectionDrivers,
			 const double& minDetectDrvCloneSz,
			 const double& extraTime,
			 const int& verbosity,
			 double& totPopSize,
			 double& em1,
			 double& em1sc,
			 // double& n_1,
			 // double& en1,
			 double& ratioForce,
			 double& currentTime,
			 int& speciesFS,
			 int& outNS_i,
			 int& iter,
			 std::vector<Genotype>& genot_out,
			 std::vector<double>& popSizes_out,
			 std::vector<int>& index_out,
			 std::vector<double>& time_out,
			 std::vector<double>& sampleTotPopSize,
			 std::vector<double>& sampleLargestPopSize,
			 std::vector<int>& sampleMaxNDr,
			 std::vector<int>& sampleNDrLargestPop,
			 bool& reachDetection,
			 std::mt19937& ran_gen,
			 // randutils::mt19937_rng& ran_gen,
			 double& runningWallTime,
			 bool& hittedWallTime,
			 const std::map<int, std::string>& intName,
			 const fitness_as_genes& genesInFitness,
			 PhylogName& phylog,
			 bool keepPhylog,
			 const fitnessEffectsAll& muEF,
			 const std::vector<int>& full2mutator,
			 const double& cPDetect,
			 const double& PDBaseline,
			 const double& checkSizePEvery,
			 const bool& AND_DrvProbExit,
			 const std::vector< std::vector<int> >& fixation_l,
			 const double& fixation_tolerance,
			 const int& min_successive_fixation,
			 const double& fixation_min_size,
			 int& ti_dbl_min,
			 int& ti_e3,
			 std::map<std::string, std::string>& lod,
			 // LOD& lod,
			 POM& pom) {

  double nextCheckSizeP = checkSizePEvery;
  const int numGenes = fitnessEffects.genomeSize;

  //smallest mutation rate, such as when all loci mutated
  double dummyMutationRate = setDummyMutationRate(mu, verbosity);
  
  genot_out.clear();

  phylog = PhylogName();
  lod.clear();
  pom = POM();

  popSizes_out.clear();
  index_out.clear();
  time_out.clear();
  totPopSize = 0.0;
  sampleTotPopSize.clear();
  currentTime = 0.0;
  iter = 0;

  //Yes, nothing yet. 
  outNS_i = -1; 

  sampleTotPopSize.clear();
  sampleLargestPopSize.clear();
  sampleMaxNDr.clear();
  sampleNDrLargestPop.clear();
 

  bool forceSample = false;
  bool simulsDone = false;
  double lastStoredSample = 0.0;


  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double timeNextPopSample;
  double tSample;

  std::vector<int> newMutations;
  int nextMutant;
  unsigned int numSpecies = 0;
  int numMutablePosParent = 0;


  unsigned int sp = 0;
 
  int iterInterrupt = 50000; //how large should we make this?

  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;

  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);

  // those to update
  int to_update = 1; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  const int sp_per_period = 5000;

 
  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

  // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;

  ti_dbl_min = 0;
  ti_e3 = 0;

 
  //McFarland
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();
  // for McFarland error
  em1 = 0.0;
  em1sc = 0.0;
  
  int lastMaxDr = 0;
  double done_at = -9;
  
  int num_successive_fixation = 0; // none so far
  

  g_min_birth_mut_ratio_nr = DBL_MAX;
  g_min_death_mut_ratio_nr = DBL_MAX;
  g_tmp1_nr = DBL_MAX;

  // Where we will store the genotypes and pops
  std::vector<Genotype> Genotypes(0);
  std::vector<spParamsP> popParams(0);
  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);
  
  // Placeholders for genotype and params created at mutation
  Genotype newGenotype = wtGenotype();
  spParamsP tmpParam;
  init_tmpP(tmpParam);

  // Initialize population from initMutants and initSize
  initPops(numSpecies, totPopSize,  outNS_i, lastStoredSample,  Genotypes,
	   popParams, genot_out, popSizes_out, index_out, time_out,
	   sampleTotPopSize, sampleLargestPopSize, sampleMaxNDr,
	   sampleNDrLargestPop, pom, ran_gen, initMutant, initSize,
	   fitnessEffects, mu, muEF, full2mutator, intName, genesInFitness,
	   dummyMutationRate, K, death, currentTime, keepEvery, mutationPropGrowth,
	   typeModel,  verbosity);

  
  // For McFL error. 
  double totPopSize_previous = totPopSize;
  double DA_previous = log1p(totPopSize_previous/K);
  
  timeNextPopSample = currentTime + sampleEvery;

  while(!simulsDone) {
    runningWallTime = difftime(time(NULL), start_time);
    if( runningWallTime > maxWallTime ) {
      hittedWallTime = true;
      forceSample = true;
      simulsDone = true;
    }

    iter++;

    // Capture use interruptions periodically
    if( !(iter % iterInterrupt)) Rcpp::checkUserInterrupt();

    //  Step  5.2 in algorithm
    message1(verbosity, "Looping through 5.2", iter, currentTime, numSpecies,
	     totPopSize, timeNextPopSample, minNextMutationTime);
    
    tSample = std::min(timeNextPopSample, finalTime);

    DEBUGfs;
    
    if(iter == 1) {
      // DO NOT move this block above (e.g., right after
      // initializing initPop): this depends on tSample
      // which depends on timeNextPopSample, which is changed inside
      // this loop.  You could move the first definition
      // of timeNextPopSample and tSample before the init
      // block but what is the point?

      // This block is kept outside of initPops as this is very
      // specific to Mather's BNB algorithm
      for(size_t m = 0; m != popParams.size(); ++m ) {
	// Do I need this pv? Probably not if I called mapTimes_updateP
	// differently (timeLastUpdate < -1) as .w. it removes the entry.
	// FIXMEmaybe: change the call?
	popParams[m].pv = mapTimes.insert(std::make_pair(-999, m));
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[m],
					   currentTime,
					   tSample,
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, m, tmpdouble1);
	popParams[m].timeLastUpdate = currentTime;
      }
    } else { // any other iter
      if(to_update == 1) {
	// we did not sample or mutate to a different species in previous period
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample,
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

	DEBUG_detect_duplicates(tmpdouble1, u_1);
	DEBUG_52(tmpdouble1, u_1, "update one");
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

	DEBUG_detect_duplicates(tmpdouble1, u_1);
	DEBUG_detect_duplicates(tmpdouble2, u_2);

	DEBUG_52(tmpdouble2, u_2, "update two, u_1");
	DEBUG_52(tmpdouble2, u_2, "update two, u_2");
	
      } else { // we sampled, so update all: i.e. to_update == 3
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;

	  DEBUG_detect_duplicates(tmpdouble1, i);
	  DEBUG_52(tmpdouble1, i, "ti_nextTime, update all");
	}
      }
    }
    if(forceSample) {
      // A VERY ugly hack. Resetting tSample to jump to sampling.
      tSample = currentTime;
      // Need this, o.w. would skip a sampling.
      timeNextPopSample = currentTime;
    }


    // Step  5.3 of the algorithm and do we sample? 
    // Find minimum to know if we need to sample the whole pop
    // We also obtain the nextMutant
    getMinNextMutationTime4(nextMutant, minNextMutationTime,
			    mapTimes);

    message1(verbosity, "after getMinNextMutationTime4",
	     iter, currentTime, popParams.size(),
	     totPopSize, timeNextPopSample, minNextMutationTime);

    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // Step 5.3 in algorithm
      currentTime = minNextMutationTime;
      // Step 5.4. in algorithm
      mutantTimeSinceLastUpdate =
	currentTime -	popParams[nextMutant].timeLastUpdate;

      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);

      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
	DEBUGfs2;
      }
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet
      // know total number of different species
      // This is a protection against things going wild. Should
      // not happen in regular usage.
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2;
	DEBUGfsnl;
      }

      if(popParams[nextMutant].numMutablePos != 0) {
	// This is the usual case. The alternative is the null or dummy
	// mutation --below.
	// Step 5.5 of algorithm

	newMutations.clear();
	obtainMutations(Genotypes[nextMutant],
			fitnessEffects,
			numMutablePosParent,
			newMutations,
			ran_gen,
			mu);

	// Step 5.6 of algorithm
	newGenotype = createNewGenotype(Genotypes[nextMutant],
					newMutations,
					fitnessEffects,
					ran_gen,
					true);
	new_sp_v(sp, newGenotype, Genotypes);

	if(sp == numSpecies) {// New species
	  ++numSpecies;
	  init_tmpP(tmpParam);
	  messageNewSpecies(verbosity, iter, numSpecies, nextMutant);

	  DEBUG_1456;
	  tmpParam.popSize = 1;

	  nr_fitness(tmpParam, popParams[nextMutant],
		     newGenotype,
		     fitnessEffects,
		     typeModel, Genotypes, popParams, currentTime);
	  if(tmpParam.birth > 0.0) {

	    tmpParam.numMutablePos = numMutablePosParent - 1;
	    tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
						    fitnessEffects,
						    mutationPropGrowth, full2mutator,
						    muEF, Genotypes, popParams, currentTime,
						    dummyMutationRate);

	    if (tmpParam.mutation > 1 ) Rcpp::Rcout << "WARNING: mutation > 1\n";
	    if ((numMutablePosParent == 1) && (verbosity >= 1)) {
	      Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	    }
	    W_f_st(tmpParam);
	    R_f_st(tmpParam);
	    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.
	    // as this is a new species
	    popParams.push_back(tmpParam);
	    Genotypes.push_back(newGenotype);
	    to_update = 2;
	    g_tmp1_nr = tmpParam.birth/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;

	    g_tmp1_nr = tmpParam.death/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;


	    // LOD: here first call to addToPhylog, with popSize popParams[sp].popSize
	    // and it is 0
	    if(keepPhylog) addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
				       intName, genesInFitness, 0);
	    addToLOD(lod, Genotypes[nextMutant], newGenotype, intName, genesInFitness);

	  } else {// fitness is 0, so we do not add it
	    --sp;
	    --numSpecies;
	    to_update = 1;
	  }
	  vvmessageNewSpecies(verbosity, sp, newGenotype, Genotypes[nextMutant],
			      tmpParam, popParams[nextMutant]);
	} else {	// A mutation to pre-existing species

	  // What we do here is step 6 of Algorithm 5, in the "Otherwise",
	  // in p. 5 of suppl mat. We will update both, and only these
	  // two.
	  to_update = 2;
	  DEBUG_1536;

	  // Could the if can be removed??
	  // Possibly. But note that the popParams[sp].popSize can be >
	  // 0, but when updated via Algo2 and added to 1.0 we can end
	  // in 1. Why? Because Algo2 can return a 0. The species
	  // "exist" in the sense that it had non-zero pop size when we
	  // last sampled/updated it.

	  if(popParams[sp].popSize > 0.0) {
	    popParams[sp].popSize = 1.0 +
	      Algo2_st(popParams[sp], currentTime, mutationPropGrowth);
	    if(verbosity >= 2) Rcpp::Rcout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  } else {
	    throw std::range_error("\n popSize == 0 but existing? \n");
	  }
	  // here one of the calls to addToPhylog, with popSize popParams[sp].popSize
	  if(keepPhylog)
	    addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
			intName, genesInFitness, popParams[sp].popSize);

	}
	// Step 5.6 of algorithm
	// u_2 irrelevant if to_update = 1;
	u_1 = nextMutant;
	u_2 = static_cast<int>(sp);
      } else { // the null or dummy mutation case
	// We increase size by 1, as we already called Algo3. And then
	// update the ti.
	++popParams[nextMutant].popSize;
	to_update = 1;
	u_1 = nextMutant;
	u_2 = -99;
	if(verbosity >= 1) Rcpp::Rcout << "Note: updating in null mutation\n";
      }
    } else { //       *********** We are sampling **********
      to_update = 3; //short_update = false;
      messageSampling(verbosity, tSample, finalTime, popParams);

      currentTime = tSample;
      nr_sample_all_pop_P(sp_to_remove,
			  popParams, Genotypes, tSample,
			  mutationPropGrowth);
      timeNextPopSample += sampleEvery;
      // When we call nr_totPopSize ... species that existed between
      // their creation and sampling time are never reflected. That is OK.
      // This is on purpose, but if you track the phylogeny, you might see
      // in the phylogeny things that never get reflected in the pops.by.time
      // object.
      if(sp_to_remove.size())
	remove_zero_sp_nr(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();

      // Check stopping conditions and fill up output structures
      nr_totPopSize_and_fill_out_crude_P(outNS_i, totPopSize,
					 lastStoredSample,
					 genot_out,
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
					 keepEvery,
					 detectionSize,
					 finalTime,
					 detectionDrivers,
					 verbosity,
					 minDetectDrvCloneSz,
					 extraTime,
					 fitnessEffects.drv,
					 cPDetect,
					 PDBaseline,
					 checkSizePEvery,
					 nextCheckSizeP,
					 ran_gen,
					 AND_DrvProbExit,
					 fixation_l,
					 fixation_tolerance,
					 min_successive_fixation,
					 fixation_min_size,
					 num_successive_fixation,
					 pom, intName,
					 genesInFitness); //keepEvery is for thinning

      messagePostSampling(verbosity, popParams, totPopSize);

      computeMcFarlandError_new(em1, em1sc, totPopSize_previous, DA_previous,
				typeModel, totPopSize, K);

      if(simulsDone) break; //skip last updateRates

      updateBirthDeathRates(popParams, Genotypes, fitnessEffects, adjust_fitness_MF,
			    K, totPopSize, currentTime, typeModel);
  
      // could go inside sample_all_pop but here we are sure death, etc, are current
      // But I catch them when they are created. Is this really needed?
      for(size_t i = 0; i < popParams.size(); i++) {
	g_tmp1_nr = popParams[i].birth/popParams[i].mutation;
	if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;

	g_tmp1_nr = popParams[i].death/popParams[i].mutation;
	if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;
      }

      forceSample = false;
    }
  }
}



// [[Rcpp::export]]
Rcpp::List nr_BNB_Algo5(Rcpp::List rFE,
			Rcpp::NumericVector mu_,
			double death,
			Rcpp::NumericVector initSize_,
			double sampleEvery,
			double detectionSize,
			double finalTime,
			int initSp,
			int initIt,
			double seed,
			int verbosity,
			int speciesFS,
			double ratioForce,
			Rcpp::CharacterVector typeFitness_,
			int maxram,
			int mutationPropGrowth,
			Rcpp::List initMutant_,
			double maxWallTime,
			double keepEvery,
			double K,
			int detectionDrivers,
			bool onlyCancer,
			bool errorHitWallTime,
			int maxNumTries,
			bool errorHitMaxTries,
			double minDetectDrvCloneSz,
			double extraTime,
			bool keepPhylog,
			Rcpp::List MMUEF,
			Rcpp::IntegerVector full2mutator_,
			double n2,
			double p2,
			double PDBaseline,
			double cPDetect_i,
			double checkSizePEvery,
			bool AND_DrvProbExit,
			Rcpp::List fixation_i) {


  precissionLoss();
  const std::vector<double> mu = Rcpp::as<std::vector<double> >(mu_);
  const std::vector < std::vector<int> >
    initMutant =  list_to_vector_of_int_vectors(initMutant_, false);
  
  const std::vector<double> initSize = Rcpp::as<std::vector<double> >(initSize_);
  
  const TypeModel typeModel = stringToModel(Rcpp::as<std::string>(typeFitness_));

  // A simple, vector-indexed way to map from numeric ids in full to
  // numeric ids in mutator. Recall all genes start with 1. So full2mutator[i-1];
  const std::vector<int> full2mutator = Rcpp::as<std::vector<int> >(full2mutator_);
  // A consistency check

  // Code for using randutils
  // randutils::mt19937_rng ran_gen;
  // if(seed == 0)
  //   ran_gen.seed();
  // else {
  //   ran_gen.seed(static_cast<unsigned int>(seed));
  //   // The next does not solve the differences between clang and gcc. So
  //   // keep it simple.
  //   // std::seed_seq s1{static_cast<unsigned int>(seed)};
  //   // ran_gen.seed(s1);
  // }

  unsigned int rseed = static_cast<unsigned int>(seed);
  if(seed == 0) {
    rseed = std::random_device{}();
  }
  std::mt19937 ran_gen(rseed);

  double cPDetect = cPDetect_i;
  if( (n2 > 0) && (p2 > 0) ) {
    if (PDBaseline <= 0) throw std::range_error("PDBaseline <= 0");
    cPDetect = set_cPDetect(n2, p2, PDBaseline);
    if(verbosity >= 1)
      Rcpp::Rcout << "  cPDetect set at " << cPDetect << "\n";
  }

  if( (K < 1 ) && ( (typeModel ==   TypeModel::mcfarlandlog) ||
		    (typeModel ==   TypeModel::mcfarlandlog) ))
    throw std::range_error("K < 1.");
  
  fitnessEffectsAll fitnessEffects =  convertFitnessEffects(rFE);
  //Used at least twice
  std::map<int, std::string> intName = mapGenesIntToNames(fitnessEffects);
  fitness_as_genes genesInFitness = fitnessAsGenes(fitnessEffects);
  PhylogName phylog;
  // LOD lod;
  std::map<std::string, std::string> lod;
  POM pom;

  // Mutator effects
	fitnessEffectsAll muEF;
  if( (full2mutator.size() != 0) ){
		muEF = convertFitnessEffects(MMUEF);
	} else {
		muEF = nullFitnessEffects();
		muEF.frequencyDependentFitness = fitnessEffects.frequencyDependentFitness;
			}
  // Paranoia. We should never end up here.
  if( (full2mutator.size() != 0) && (muEF.genomeSize == 0))
    throw std::logic_error("full2mutator > 0 with mutatorEffects.genomesize 0");
  if( (full2mutator.size() == 0) && (muEF.genomeSize != 0)) {
    throw std::logic_error("full2mutator 0 with mutatorEffects.genomesize != 0");
  }

  // Fixation stuff: run until some genotype combinations fixed
  double fixation_tolerance = -9;
  int min_successive_fixation = 100;
  double fixation_min_size = 0.0;
  std::vector < std::vector<int> > fixation_l;

  if( fixation_i.size() != 0 ) {
    Rcpp::List fggl = fixation_i["fixation_list"] ;
    fixation_l = list_to_vector_of_int_vectors(fggl, true); 
    fixation_tolerance = Rcpp::as<double>(fixation_i["fixation_tolerance"]);
    min_successive_fixation = Rcpp::as<int>(fixation_i["min_successive_fixation"]);
    fixation_min_size = Rcpp::as<double>(fixation_i["fixation_min_size"]);
  } else { 
     fixation_l.resize(0); // explicit
  }


  bool runAgain = true;
  bool reachDetection = false;
  //Output
  std::vector<Genotype> genot_out;
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
  std::vector<int> sampleMaxNDr; //The largest number of drivers in any
				 //genotype or clone at each  time sample
  std::vector<int> sampleNDrLargestPop; //Number of drivers in the clone
					// with largest size (at each time
					// sample)
  sampleTotPopSize.reserve(initIt);
  sampleLargestPopSize.reserve(initIt);
  sampleMaxNDr.reserve(initIt);
  sampleNDrLargestPop.reserve(initIt);

  int outNS_i = -1; // the column in the outNS
  time_t start_time = time(NULL);
  double runningWallTime = 0;
  bool  hittedWallTime = false;
  bool hittedMaxTries = false;

  double em1, em1sc; // new computation of McFarland error
  em1 = 0.0;
  em1sc = 0.0;


  // 5.1 Initialize some vars (NOT the population)

  int numRuns = 0;
  int numRecoverExcept = 0;
  bool forceRerun = false;

  double currentTime = 0;
  int iter = 0;

  int ti_dbl_min = 0;
  int ti_e3 = 0;

  int accum_ti_dbl_min = 0;
  int accum_ti_e3 = 0;
  
  while(runAgain) {

    if(numRuns >= maxNumTries) {
      // We want this here to avoid an extra run and confusing output
      hittedMaxTries = true;
      Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
      runAgain = false;
      if(errorHitMaxTries) {
	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
	return
	  List::create(Named("HittedWallTime") = false,
		       Named("HittedMaxTries") = true,
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
      }
      break;
    }

    try {
      Rcpp::checkUserInterrupt();
      // it is CRUCIAL that several entries are zeroed (or -1) at the
      // start of innerBNB when we do multiple runs if onlyCancer = true.
      nr_innerBNB(fitnessEffects,
		  initSize,
		  K,
		  typeModel,
		  mutationPropGrowth,
		  mu,
		  death,
		  keepEvery,
		  sampleEvery,
		  initMutant,
		  start_time,
		  maxWallTime,
		  finalTime,
		  detectionSize,
		  detectionDrivers,
		  minDetectDrvCloneSz,
		  extraTime,
		  verbosity,
		  totPopSize,
		  em1,
		  em1sc,
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
		  ran_gen,
		  runningWallTime,
		  hittedWallTime,
		  intName,
		  genesInFitness,
		  phylog,
		  keepPhylog,
		  muEF,
		  full2mutator,
		  cPDetect,
		  PDBaseline,
		  checkSizePEvery,
		  AND_DrvProbExit,
		  fixation_l,
		  fixation_tolerance,
		  min_successive_fixation,
		  fixation_min_size,
		  ti_dbl_min,
		  ti_e3,
		  lod,
		  pom);
      ++numRuns;
      forceRerun = false;
      accum_ti_dbl_min += ti_dbl_min;
      accum_ti_e3 += ti_e3;
    } catch (rerunExcept &e) {
      Rcpp::Rcout << "\n Recoverable exception " << e.what()
		  << ". Rerunning.";
      forceRerun = true;
      ++numRecoverExcept;
      ++numRuns; // exception should count here!
      accum_ti_dbl_min += ti_dbl_min;
      accum_ti_e3 += ti_e3;
    } catch (const std::exception &e) {
      Rcpp::Rcout << "\n Unrecoverable exception: " << e.what()
		  << ". Aborting. \n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = e.what()));
    } catch (...) {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting."
		  << "(User interrupts also generate this).\n";
      return
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = "Unknown exception"));
    }
    if(hittedWallTime) {
      Rcpp::Rcout << "\n Hitted wall time. Exiting.";
      runAgain = false;
      if(errorHitWallTime) {
	Rcpp::Rcout << "\n Hitting wall time is regarded as an error. \n";
	return
	  List::create(Named("HittedWallTime") = true,
		       Named("HittedMaxTries") = false, // yes, for
							// coherent return
							// objects
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
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
    DEBUG_rrr;
  } // runAgain loop


  std::vector<std::vector<int> > genot_out_v = genot_to_vectorg(genot_out);
  std::vector<std::vector<int> > uniqueGenotypes_vector_nr  =
    uniqueGenot_vector(genot_out_v);
  IntegerMatrix returnGenotypes =
    nr_create_returnGenotypes(fitnessEffects.genomeSize,
  			      uniqueGenotypes_vector_nr);
  Rcpp::NumericMatrix outNS = create_outNS(uniqueGenotypes_vector_nr,
  					   genot_out_v,
  					   popSizes_out,
  					   index_out, time_out,
  					   outNS_i, maxram);

  int maxNumDrivers = 0;
  int totalPresentDrivers = 0;
  std::vector<int> countByDriver(fitnessEffects.drv.size(), 0);
  std::vector<int> presentDrivers;
  driverCounts(maxNumDrivers, totalPresentDrivers,
	       countByDriver, presentDrivers,
	       returnGenotypes, fitnessEffects.drv);


  std::vector<std::string> genotypesAsStrings =
    genotypesToNameString(uniqueGenotypes_vector_nr, genesInFitness, intName);
  std::string driversAsString =
    driversToNameString(presentDrivers, intName);


  std::vector<double> sampleLargestPopProp(outNS_i + 1);
  if((outNS_i + 1) != static_cast<int>(sampleLargestPopSize.size()))
    throw std::length_error("outNS_i + 1 != sampleLargestPopSize.size");
  std::transform(sampleLargestPopSize.begin(), sampleLargestPopSize.end(),
  		 sampleTotPopSize.begin(),
  		 sampleLargestPopProp.begin(),
  		 std::divides<double>());
  NumericMatrix perSampleStats(outNS_i + 1, 5);
  fill_SStats(perSampleStats, sampleTotPopSize, sampleLargestPopSize,
  	      sampleLargestPopProp, sampleMaxNDr, sampleNDrLargestPop);

  // create the lod return pieces. Move to a function later
  std::vector<std::string> lod_parent;
  std::vector<std::string> lod_child;
  for (const auto &l : lod) {
    lod_child.push_back(l.first);
    lod_parent.push_back(l.second);
  }

  return
    List::create(Named("pops.by.time") = outNS,
		 Named("NumClones") = uniqueGenotypes_vector_nr.size(),
		 Named("TotalPopSize") = totPopSize,
		 Named("Genotypes") = returnGenotypes,
		 Named("GenotypesWDistinctOrderEff") = Rcpp::wrap(uniqueGenotypes_vector_nr),
		 Named("GenotypesLabels") = Rcpp::wrap(genotypesAsStrings),
		 Named("MaxNumDrivers") = maxNumDrivers,
		 Named("MaxDriversLast") = sampleMaxNDr[outNS_i],
		 Named("NumDriversLargestPop") =  sampleNDrLargestPop[outNS_i],
		 Named("LargestClone") = sampleLargestPopSize[outNS_i],
		 Named("PropLargestPopLast") = sampleLargestPopProp[outNS_i],
		 Named("FinalTime") = currentTime,
		 Named("NumIter") = iter,
		 Named("HittedWallTime") = hittedWallTime,
		 Named("HittedMaxTries") = hittedMaxTries,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 Named("OccurringDrivers") = driversAsString,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = List::create(Named("attemptsUsed") = numRuns,
					       Named("errorMF") =
					       returnMFE_new(em1sc, typeModel),
					       Named("errorMF_size") =
					       returnMFE_new(em1, typeModel), 
					       Named("minDMratio") =
					       g_min_death_mut_ratio_nr,
					       Named("minBMratio") =
					       g_min_birth_mut_ratio_nr,
					       Named("PhylogDF") =  DataFrame::create(
										      Named("parent") = phylog.parent,
										      Named("child") = phylog.child,
										      Named("time") = phylog.time,
										      Named("pop_size_child") = phylog.pop_size_child
										      ),
					       Named("UnrecoverExcept") = false,
					       Named("numRecoverExcept") = numRecoverExcept,
					       Named("accum_ti_dbl_min") = accum_ti_dbl_min,
					       Named("accum_ti_e3") = accum_ti_e3,
					       Named("LOD_DF") = DataFrame::create(
										   Named("parent") = lod_parent, 
										   Named("child") = lod_child 
										   ),
					       Named("POM") = Rcpp::wrap(pom.genotypesString)
					       )
		 );
}


// Creating return object:

// The 0, 1 representation is how most of the work is done in R: do I want
// to change that?

// Order: beware of two things: order is important for the "true"
// genotypes, but is not immediately observable. So for 0,1
// representation, not needed or used. Thus, maybe I want two
// representations.

// Yes, the full Genotye structure is only used when assigning fitness. So
// could we use a collapsed one with: order + rest? Nope, as whenever I'd
// create a child from a genotype, I'd need like deconvolve, and go back
// to the three piece structure. This seems much more expensive than the
// overloaded == and the usage of the overloaded < (this is only used at
// the end, when producing the output objects)



// When you need to alter things that affect a genotype periodically, like
// the FDF, or the McFarland rates, etc, you will want to do it here:
// when we call updateRatesMcFarlandLog and similar, e.g.,
// In BNB_nr.cpp, in nr_innerBNB function
// when we are sampling.
// Probably also where nr_fitness is called
// And generally, in general, were we pass currentTime



// Old notes about models
// Exp and McFarland and McFarlandlog are also like Datta et al., 2013
// An additional driver gene mutation increases a cell’s fitness by a
// factor of (1+sd), whereas an additional housekeeper gene mutation
// decreases fitness by a factor of (1-sh) and the effect of multiple
// mutations is multiplicative
