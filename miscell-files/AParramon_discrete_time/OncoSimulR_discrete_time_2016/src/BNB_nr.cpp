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



// #include "OncoSimul.h"
// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include "debug_common.h"
#include "common_classes.h"
#include "bnb_common.h"
#include "new_restrict.h"
#include <Rcpp.h>
#include <limits>
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

// To track if mutation is really much smaller than birth/death
#define MIN_RATIO_MUTS_NR
#ifdef MIN_RATIO_MUTS_NR
// There is really no need for these to be globals?
// Unless I wanted to use them inside some function. So leave as globals.
double g_min_birth_mut_ratio_nr = DBL_MAX;
double g_min_death_mut_ratio_nr = DBL_MAX;
double g_tmp1_nr = DBL_MAX;
#endif


void nr_fitness(spParamsP& tmpP,
		const spParamsP& parentP,
		const Genotype& ge,
		const fitnessEffectsAll& F,
		const TypeModel typeModel) {
		       // const double& genTime,
		       // const double& adjust_fitness_B,
		       // const double& adjust_fitness_MF) {

  // We want a way to signal immediate non-viability of a clone. For
  // "birth-based" models that happens when any s = -1, as the fitness is
  // 0. By setting birth = 0.0 we ensure this clone does not get added and
  // we never reach into algo2, etc, leading to numerical problems.

  // With Bozic models, which are "death-based", it is different. For
  // bozic2, birth is bounded, so any death > 2 would lead to birth <
  // 0. For bozic1, deaths of around 50 lead to numerical issues.  The
  // general rule is: set those mutations to -inf, so prodDeathFitness
  // returns an inf for death, and that is recognized as "no
  // viability" (anything with death > 99)

  // The ones often used are bozic1, exp, mcfarlandlog

  if(typeModel == TypeModel::bozic1) {
    tmpP.death = prodDeathFitness(evalGenotypeFitness(ge, F));
    if( tmpP.death > 99) {
      tmpP.birth = 0.0; 
    } else {
      tmpP.birth = 1.0;
    }
  // } else if (typeModel == TypeModel::bozic2) {
  //   double pp = prodDeathFitness(evalGenotypeFitness(ge, F));
  //   tmpP.birth = std::max(0.0, (1.0/genTime) * (1.0 - 0.5 * pp ));
  //   tmpP.death = (0.5/genTime) * pp;
  } else {
    double fitness = prodFitness(evalGenotypeFitness(ge, F));
    if( fitness <= 0.0) {
      tmpP.absfitness = 0.0;
      tmpP.death = 1.0;
      tmpP.birth = 0.0; 
    } else{
      // Set appropriate defaults and change only as needed
      tmpP.death = parentP.death;
      tmpP.absfitness = parentP.absfitness;
      tmpP.birth = fitness;
      // exp, mcfarland, and mcfarlandlog as above. Next are the two exceptions.
      // if(typeModel == TypeModel::beerenwinkel) {
      // 	tmpP.absfitness = fitness; 
      // 	tmpP.birth = adjust_fitness_B * tmpP.absfitness;
      // } else if(typeModel == TypeModel::mcfarland0) {
      // 	tmpP.absfitness = fitness;
      // 	tmpP.birth = adjust_fitness_MF * tmpP.absfitness;
      // }
    }
  }
  // Exp and McFarland and McFarlandlog are also like Datta et al., 2013
  // An additional driver gene mutation increases a cell’s fitness by a
  // factor of (1+sd), whereas an additional housekeeper gene mutation
  // decreases fitness by a factor of (1-sh) and the effect of multiple
  // mutations is multiplicative
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

inline unsigned int new_sp(const Genotype& newGenotype,
		    const std::vector<Genotype> Genotypes) {
  for(unsigned int sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] ) {
      return sp;
    }
  }
  return Genotypes.size();
}

void driverCounts(int& maxNumDrivers,
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

// FIXME: why not keep the number of present drivers in the genotype? We
// call often the getGenotypeDrivers(ge, drv).size()


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



void nr_totPopSize_and_fill_out_crude_P(int& outNS_i,
					double& totPopSize, 
					double& lastStoredSample,
					std::vector<Genotype>& genot_out,
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
					const double& fatalPopSize = 1e15) {
  // Fill out, but also compute totPopSize
  // and return sample summaries for popsize, drivers.

  
  

  
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

  
  if (keepEvery < 0) {
    storeThis = false;
  } else if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }

  // FIXME
  // this is the usual exit condition
  // (totPopSize >= detectionSize) ||
  // 	  ( (lastMaxDr >= detectionDrivers) &&
  // 	    (popSizeOverDDr >= minDetectDrvCloneSz)

  // Now add the prob. of exiting.

  // Doing this is cheaper than drawing unnecessary runifs.
  // Equality, below, leads to suprises with floating point arith.

  // Operates the same as we do with keepEvery, but here we
  // compute the jump in each accepted sample. And here we use >, not
  // >=
  if(currentTime > nextCheckSizeP) {
    checkSizePNow = true;
    nextCheckSizeP = currentTime + checkSizePEvery;
    // Nope; minimal jump can be smaller than checkSizePEvery
    // nextCheckSizeP += checkSizePEvery;
  } else {
    checkSizePNow = false;
  }
  
  if(extraTime > 0) {
    if(done_at <  0) {
      if( (totPopSize >= detectionSize) ||
	  ( (lastMaxDr >= detectionDrivers) &&
	    (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	  ( checkSizePNow &&
	    detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen))) {
	done_at = currentTime + extraTime;
      }
    } else if (currentTime >= done_at) {
	simulsDone = true;
	reachDetection = true; 
      }
  } else if( (totPopSize >= detectionSize) ||
	     ( (lastMaxDr >= detectionDrivers) &&
	       (popSizeOverDDr >= minDetectDrvCloneSz) ) ||
	     ( checkSizePNow &&
	       detectedSizeP(totPopSize, cPDetect, PDBaseline, ran_gen)) ) {
    simulsDone = true;
    reachDetection = true; 
  }
  
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
	ndr_lp = getGenotypeDrivers(Genotypes[i], drv).size();
      }
    }
    sampleTotPopSize.push_back(totPopSize);
    sampleLargestPopSize.push_back(l_pop_s);
    sampleMaxNDr.push_back(max_ndr);
    sampleNDrLargestPop.push_back(ndr_lp);
  } 
  
  if( !std::isfinite(totPopSize) ) {
    throw std::range_error("totPopSize not finite");
  }
  if( std::isnan(totPopSize) ) {
    throw std::range_error("totPopSize is NaN");
  }
  

}

// FIXME: I might want to return the actual drivers in each period
// and the actual drivers in the population with largest popsize
// Something like what we do now with whichDrivers
// and count_NumDrivers



inline void nr_reshape_to_outNS(Rcpp::NumericMatrix& outNS,
				const vector<vector<int> >& uniqueGenotV,
				const vector<vector<int> >& genot_out_v,
				const vector<double>& popSizes_out,
				const vector<int>& index_out,
				const vector<double>& time_out){
  
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



// FIXME: when creating the 0/1, collapse those that are the same


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

// No longer used.
// std::string genotypeToIntString(const std::vector<int>& genotypeV,
// 				   const fitness_as_genes& fg) {
  
//   // The genotype vectors are returned as a string of ints.
  
//   std::string strGenotype;

//   std::vector<int> order_int;
//   std::vector<int> rest_int;

//   for(auto const &g : genotypeV) {
//     if( binary_search(fg.orderG.begin(), fg.orderG.end(), g)) {
//       order_int.push_back(g);
//     } else {
//       rest_int.push_back(g);
//     }
//   }

//   std::string order_sep = "_";
//   std::string order_part;
//   std::string rest;
//   std::string comma = "";

  
//   for(auto const &g : order_int) {
// #ifdef _WIN32  
//      order_part += (comma + SSTR(g));
// #endif
// #ifndef _WIN32
//     order_part += (comma + std::to_string(g));
// #endif
//     comma = ", ";
//   }
//   comma = "";
//   for(auto const &g : rest_int) {
// #ifdef _WIN32  
//      rest += (comma + SSTR(g));
// #endif
// #ifndef _WIN32
//     rest += (comma + std::to_string(g));
// #endif
//     comma = ", ";
//   }
//   if(fg.orderG.size()) {
//     strGenotype = order_part + order_sep + rest;
//   } else {
//     strGenotype = rest;
//   }
//   return strGenotype;
// }


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
  // FIXME: when sure no problems, remove at if needed for speed.
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
					       // const fitnessEffectsAll& F,
					       const std::map<int, std::string>& intName) {
  //fitness_as_genes fg = fitnessAsGenes(F); // I use this before;
  std::vector<std::string> gs;
  for(auto const &v: uniqueGenotypesV )
      gs.push_back(genotypeToNameString(v, fg, intName));
  return gs;
}


// std::vector<std::string> genotypesToString(const std::vector< vector<int> >& uniqueGenotypesV,
// 					   const fitnessEffectsAll& F,
// 					   bool names = true) {
//   fitness_as_genes fg = fitnessAsGenes(F);
//   std::vector<std::string> gs;

//   if(names) {
//     std::map<int, std::string> intName = mapGenesIntToNames(F);
//     for(auto const &v: uniqueGenotypesV )
//       gs.push_back(genotypeToNameString(v, fg, intName));
//   } else {
//       for(auto const &v: uniqueGenotypesV )
// 	gs.push_back(genotypeToIntString(v, fg));
//   }
  
//   // exercise: do it with lambdas
//   // std::transform(uniqueGenotypesV.begin(), uniqueGenotypesV.end(),
//   // 		 back_inserter(gs), vectorGenotypeToString);
//   return gs;
// }

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
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
		<< "\n      Genotype = ";
    print_Genotype(Genotypes[i]); //genotypeSingleVector(Genotypes[i])
    //	      << "\n      sp_id = " << genotypeSingleVector(Genotypes[i]) // sp_id[i]  
    Rcpp::Rcout << "\n      pre-update popSize = " 
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
      sp_to_remove.push_back(i);

#ifdef DEBUGV
      Rcpp::Rcout << "\n\n     Removing species i = " << i 
		  << " with genotype = ";
      print_Genotype(Genotypes[i]); //genotypeSingleVector(Genotypes[i]);
#endif
    } 
#ifdef DEBUGV
    Rcpp::Rcout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
  }
}


void addToPhylog(PhylogName& phylog,
		 const Genotype& parent,
		 const Genotype& child,
		 double time,
		 const std::map<int, std::string>& intName,
		 const fitness_as_genes& fg) {
  phylog.time.push_back(time);
  phylog.parent.push_back(genotypeToNameString(genotypeSingleVector(parent),
					       fg, intName));
  phylog.child.push_back(genotypeToNameString(genotypeSingleVector(child),
					      fg, intName));
}



static void nr_innerBNB(const fitnessEffectsAll& fitnessEffects,
			const double& initSize,
			const double& K,
			// const double& alpha,
			// const double& genTime,
			const TypeModel typeModel,
			const int& mutationPropGrowth,
			const std::vector<double>& mu,
			// const double& mu,
			const double& death,
			const double& keepEvery,
			const double& sampleEvery,		     
			const std::vector<int>& initMutant,
			const time_t& start_time,
			const double& maxWallTime,
			const double& finalTime,
			const double& detectionSize,
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
			const double& checkSizePEvery) {
  
  double nextCheckSizeP = checkSizePEvery;
  const int numGenes = fitnessEffects.genomeSize;

  double mymindummy = 1.0e-11; //1e-10
  double targetmindummy = 1.0e-10; //1e-9
  double minmu = *std::min_element(mu.begin(), mu.end());
  // Very small, but no less than mymindummy, for numerical issues.
  // We can probably go down to 1e-13. 1e-16 is not good as we get lots
  // of pE.f not finite. 1e-15 is probably too close, and even if no pE.f
  // we can get strange behaviors.
  double dummyMutationRate = std::max(std::min(minmu/1.0e4, targetmindummy),
				      mymindummy);
  // This should very rarely happen:
  if(minmu <= 1e-9 ) {
    double newdd = minmu/100.0;
    Rcpp::Rcout << "WARNING: the smallest mutation rate is "
		<< "<= " << mymindummy << ". That is a really small value"
		<< "(per-base mutation rate in the human genome is"
		<< " ~ 1e-11 to 1e-9). "
		<< "Setting dummyMutationRate to your min/100 = "
		<< newdd
		<< ". There can be numerical problems later.\n";
    dummyMutationRate = newdd;
  }
  // double dummyMutationRate = 1e-10;
  // ALWAYS initialize this here, or reinit or rezero
  genot_out.clear();

  phylog = PhylogName();
  
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

  std::vector<int> newMutations;
  int nextMutant;
  unsigned int numSpecies = 0;
  int numMutablePosParent = 0;


  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 100; 
  //int timeL = 1000;

  int iterInterrupt = 50000; //how large should we make this?
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);

  // those to update
  int to_update = 1; //1: one species; 2: 2 species; 3: all.
  int u_1 = -99;
  int u_2 = -99;

  Genotype newGenotype;
  std::vector<Genotype> Genotypes(1);
  Genotypes[0] = wtGenotype(); //Not needed, but be explicit.
  
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

#ifdef MIN_RATIO_MUTS_NR
  g_min_birth_mut_ratio_nr = DBL_MAX;
  g_min_death_mut_ratio_nr = DBL_MAX;
  g_tmp1_nr = DBL_MAX;
#endif

      // // FIXME debug
      // Rcpp::Rcout << " popSize[0]  at 1b ";
      // print_spP(popParams[0]);
      // // end debug

    // This long block, from here to X1, is ugly and a mess!
  // This is what takes longer to figure out whenever I change
  // anything. FIXME!!
  if(initMutant.size() > 0) {
    Genotypes[0] = createNewGenotype(wtGenotype(),
				     initMutant,
				     fitnessEffects,
				     ran_gen,
				     false);
    int numGenesInitMut = Genotypes[0].orderEff.size() +
      Genotypes[0].epistRtEff.size() + Genotypes[0].rest.size();
    int numGenesGenotype = fitnessEffects.allGenes.size();
    popParams[0].numMutablePos = numGenesGenotype - numGenesInitMut;
    // Next two are unreachable since caught in R.
    // But just in case, since it would lead to seg fault.
    if(popParams[0].numMutablePos < 0)
      throw std::invalid_argument("initMutant's genotype has more genes than are possible.");
    if(popParams[0].numMutablePos == 0)
      throw std::invalid_argument("initMutant has no mutable positions: genotype with all genes mutated.");
    // popParams[0].numMutablePos = numGenes - 1;
    // From obtainMutations, but initMutant an int vector. But cumbersome.
    // std::vector<int> sortedg = convertGenotypeFromInts(initMutant);
    // sort(sortedg.begin(), sortedg.end());
    // std::vector<int> nonmutated;
    // set_difference(fitnessEffects.allGenes.begin(), fitnessEffects.allGenes.end(),
    // 		   sortedg.begin(), sortedg.end(),
    // 		   back_inserter(nonmutated));
    // popParams[0].numMutablePos = nonmutated.size();



    // Commenting out the unused models!
    // if(typeModel == TypeModel::beerenwinkel) {
      
    //   popParams[0].death = 1.0; //note same is in McFarland.
    //   // But makes sense here; adjustment in beerenwinkel is via fitness
      
    //   // initialize to prevent birth/mutation warning with Beerenwinkel
    //   // when no mutator. O.w., the defaults
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
    // 			currentTime, alpha, initSize, 
    // 			mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland0) {
    //   // death equal to birth of a non-mutant.
    //   popParams[0].death = log1p(totPopSize/K); // log(2.0), except rare cases
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
    // 			    totPopSize,
    // 			    mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland) {
    //   popParams[0].death = totPopSize/K;
    //   popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
    // 								fitnessEffects));
    // } else       if(typeModel == TypeModel::mcfarlandlog) {

    if(typeModel == TypeModel::mcfarlandlog) {
      popParams[0].death = log1p(totPopSize/K);
      popParams[0].birth = prodFitness(evalGenotypeFitness(Genotypes[0],
								fitnessEffects));
    } else if(typeModel == TypeModel::bozic1) {
      tmpParam.birth =  1.0;
      tmpParam.death = -99.9;
    // } else if (typeModel == TypeModel::bozic2) {
    //   tmpParam.birth =  -99;
    //   tmpParam.death = -99;
    } else if (typeModel == TypeModel::exp) {
      tmpParam.birth =  -99;
      tmpParam.death = death;
    } else {
      // caught in R, so unreachable here
      throw std::invalid_argument("this ain't a valid typeModel");
    } 
    // if( (typeModel != TypeModel::beerenwinkel) && (typeModel != TypeModel::mcfarland0) 
    // 	&& (typeModel != TypeModel::mcfarland) && (typeModel != TypeModel::mcfarlandlog)) // wouldn't matter
    //   nr_fitness(popParams[0], tmpParam,
    // 		 Genotypes[0],
    // 		 fitnessEffects,
    // 		 typeModel, genTime,
    // 		 adjust_fitness_B, adjust_fitness_MF);
    if( (typeModel != TypeModel::mcfarlandlog)) // wouldn't matter
      nr_fitness(popParams[0], tmpParam,
		 Genotypes[0],
		 fitnessEffects,
		 typeModel);
    // , genTime);
    //		 adjust_fitness_B, adjust_fitness_MF);
    // we pass as the parent the tmpParam; it better initialize
    // everything right, or that will blow. Reset to init
    init_tmpP(tmpParam);
  } else { //no initMutant
    popParams[0].numMutablePos = numGenes;
    // if(typeModel == TypeModel::beerenwinkel) {
    //   popParams[0].death = 1.0;
    //   // initialize to prevent birth/mutation warning with Beerenwinkel
    //   // when no mutator. O.w., the defaults
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesBeeren(popParams, adjust_fitness_B, initSize,
    // 			currentTime, alpha, initSize, 
    // 			mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland0) {
    //   popParams[0].death = log1p(totPopSize/K);
    //   if(!mutationPropGrowth)
    // 	popParams[0].mutation = mu * popParams[0].numMutablePos;
    //   popParams[0].absfitness = 1.0;
    //   updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
    // 			    totPopSize,
    // 			    mutationPropGrowth, mu);
    // } else if(typeModel == TypeModel::mcfarland) {
    //   popParams[0].birth = 1.0;
    //   popParams[0].death = totPopSize/K;
    //   // no need to call updateRates
    // } else if(typeModel == TypeModel::mcfarlandlog) {
    if(typeModel == TypeModel::mcfarlandlog) {
      popParams[0].birth = 1.0;
      popParams[0].death = log1p(totPopSize/K);
      // no need to call updateRates
    } else if(typeModel == TypeModel::bozic1) {
       popParams[0].birth = 1.0;
       popParams[0].death = 1.0;
    // } else if (typeModel == TypeModel::bozic2) {
    //   popParams[0].birth = 0.5/genTime;
    //   popParams[0].death = 0.5/genTime;
    } else if (typeModel == TypeModel::exp) {
      popParams[0].birth = 1.0;
      popParams[0].death = death;
    } else {
      throw std::invalid_argument("this ain't a valid typeModel");
    }
  }


  

  
  // // these lines (up to, and including, R_F_st)
  // // not needed with mcfarland0 or beerenwinkel
  // if(mutationPropGrowth)
  //   popParams[0].mutation = mu * popParams[0].birth * popParams[0].numMutablePos;
  // else
  //   popParams[0].mutation = mu * popParams[0].numMutablePos;

  popParams[0].mutation = mutationFromScratch(mu, popParams[0], Genotypes[0],
					      fitnessEffects, mutationPropGrowth,
					      full2mutator, muEF);  
  W_f_st(popParams[0]);
  R_f_st(popParams[0]);

  
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
    sampleMaxNDr.push_back(getGenotypeDrivers(Genotypes[0],
					      fitnessEffects.drv).size());
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




  
  while(!simulsDone) {


    // Check how we are doing with time as first thing.
    runningWallTime = difftime(time(NULL), start_time);
    if( runningWallTime > maxWallTime ) {
      hittedWallTime = true;
      forceSample = true;
      simulsDone = true;
    }
    
    iter++;
    
    if( !(iter % iterInterrupt))
      Rcpp::checkUserInterrupt();
    
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
      Rcpp::Rcout <<"\n\n\n*** Looping through 5.2. Iter = " << iter
		  << ".  Current time " << currentTime <<	" \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    Rcpp::Rcout << " DEBUGV\n";
    Rcpp::Rcout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    if(iter == 1) { // handle special case of first iter
    tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample, 
					 ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      //popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(to_update == 1) {
	// we did not sample or mutate to a different species in previous period
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, 
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

#ifdef DEBUGV
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update one ******\n For to_update = \n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_1]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_1].popSize 
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
	Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update two ******\n " 
		  << "     tSample  = " << tSample
	    
		  << "\n\n**   Species  = " << u_1 
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_1]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_1].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble1
		  << " \n     species R " << popParams[u_1].R
		  << " \n     species W " << popParams[u_1].W
		  << " \n     species death " << popParams[u_1].death
		  << " \n     species birth " << popParams[u_1].birth


		  << "\n\n**     Species  = " << u_2 
		    << "\n       genotype =  ";
	print_Genotype(Genotypes[u_2]);
	Rcpp::Rcout << "\n       popSize = " << popParams[u_2].popSize 
		  << "\n       currentTime = " << currentTime 
		  << "\n       popParams[i].nextMutationTime = " 
		  << tmpdouble2
		  << " \n     species R " << popParams[u_2].R
		  << " \n     species W " << popParams[u_2].W
		  << " \n     species death " << popParams[u_2].death
		  << " \n     species birth " << popParams[u_2].birth;
#endif

      } else { // we sampled, so update all: i.e. to_update == 3
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;
	  
#ifdef DEBUGV
	  Rcpp::Rcout << "\n\n     ********* 5.2: call to ti_nextTime, update all ******\n " 
		    << "     Species  = " << i 
		      << "\n       genotype =  ";
	  print_Genotype(Genotypes[i]);
	  Rcpp::Rcout << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
		      << "\n       popParams[i].nextMutationTime = " 
		      << tmpdouble1
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
      // know total number of different species
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

      if(popParams[nextMutant].numMutablePos != 0) {
	// this is the usual case. The alternative is the dummy or null mutation

      
	// ************   5.5   ***************

	newMutations.clear();
	// FIXME: nonmutated also returned here
	obtainMutations(Genotypes[nextMutant],
			fitnessEffects,
			numMutablePosParent,
			newMutations,
			ran_gen,
			mu);
	// nr_change
	// getMutatedPos_bitset(mutatedPos, numMutablePosParent, // r,
	// 		     ran_gen,
	// 		     mutablePos,
	// 		     Genotypes[nextMutant], 
	// 		     numGenes);
      
	// ************   5.6   ***************
	newGenotype = createNewGenotype(Genotypes[nextMutant],
					newMutations,
					fitnessEffects,
					ran_gen,
					true);
	// nr_change
	// newGenotype = Genotypes[nextMutant];
	// newGenotype.set(mutatedPos);
	// newGenotype[mutatedPos] = 1;

	// FIXME
	// any speed diff between a) and b)?
	// a)
	new_sp_v(sp, newGenotype, Genotypes);
	// b)
	// sp = 0;
	// sp = new_sp(newGenotype, Genotypes);
	
	// nr_change
	// new_sp_bitset(sp, newGenotype, Genotypes);

	if(sp == numSpecies) {// New species
	  ++numSpecies;
	  init_tmpP(tmpParam);

	  if(verbosity >= 2) {
	    Rcpp::Rcout <<"\n     Creating new species   " << (numSpecies - 1)
			<< "         from species "  <<   nextMutant;
	  }
	
	  tmpParam.popSize = 1;

	  nr_fitness(tmpParam, popParams[nextMutant],
		     newGenotype,
		     fitnessEffects,
		     typeModel);// , genTime,
		     // adjust_fitness_B, adjust_fitness_MF);

	  if(tmpParam.birth > 0.0) {
	    // if(keepMutationTimes)
	    //   update_mutation_freqs(newMutation, currentTime, mutation_freq_at);
	    //FIXME: phylog
	    if(keepPhylog)
	      addToPhylog(phylog, Genotypes[nextMutant], newGenotype, currentTime,
			  intName, genesInFitness);
	    
	    tmpParam.numMutablePos = numMutablePosParent - 1;
	    tmpParam.mutation = mutationFromScratch(mu, tmpParam, newGenotype,
					       fitnessEffects,
					       mutationPropGrowth, full2mutator,
						    muEF);
	    // tmpParam.mutation = mutationFromParent(mu, tmpParam, popParams[nextMutant],
	    // 					   newMutations, mutationPropGrowth,
	    // 					   newGenotype, full2mutator,
	    // 					   muEF);

	    
	    //tmpParam.mutation = mu * (numMutablePosParent - 1);
	    if (tmpParam.mutation > 1 )
	      Rcpp::Rcout << "WARNING: mutation > 1\n";
	    if (numMutablePosParent == 1) {
	      if(verbosity >= 1)
		Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	      // FIXME:varmutrate: give the value of dummy here.
	      tmpParam.mutation = dummyMutationRate; // dummy mutation here. Set some mu.
	    }
	    W_f_st(tmpParam);
	    R_f_st(tmpParam);
	    tmpParam.timeLastUpdate = -99999.99999; //mapTimes_updateP does what it should.
	    // as this is a new species
	    popParams.push_back(tmpParam);
	    Genotypes.push_back(newGenotype);
	    to_update = 2;
#ifdef MIN_RATIO_MUTS_NR
	    g_tmp1_nr = tmpParam.birth/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;
	  
	    g_tmp1_nr = tmpParam.death/tmpParam.mutation;
	    if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;	
#endif	  
	  } else {// fitness is 0, so we do not add it
	    --sp;
	    --numSpecies;
	    to_update = 1;
	  }
	  // #ifdef DEBUGV	
	  if(verbosity >= 3) {
	    Rcpp::Rcout << " \n\n\n Looking at NEW species " << sp << " at creation";
	    Rcpp::Rcout << "\n New Genotype :";
	    print_Genotype(newGenotype);
	    Rcpp::Rcout << "\n Parent Genotype :";
	    print_Genotype(Genotypes[nextMutant]);
	    // Rcpp::Rcout << "\n Genotype = " << genotypeSingleVector(newGenotype); //Genotypes[sp];
	    //Genotypes[sp].to_ullong();
	    Rcpp::Rcout << "\n birth of sp = " << tmpParam.birth;
	    Rcpp::Rcout << "\n death of sp = " << tmpParam.death;
	    // Rcpp::Rcout << "\n s = " << s;
	    Rcpp::Rcout << "\n parent birth = " << popParams[nextMutant].birth;
	    Rcpp::Rcout << "\n parent death = " << popParams[nextMutant].death;
	    // Rcpp::Rcout << "\n parent Genotype = " << genotypeSingleVector(Genotypes[nextMutant]);
	    Rcpp::Rcout << "\n\n popParams parent: \n";
	    print_spP(popParams[nextMutant]);
	    Rcpp::Rcout << "\n\npopParams child: \n";
	    print_spP(tmpParam);
	    }
	  // #endif
	} else {	// A mutation to pre-existing species

	  // What we do here is step 6 of Algorithm 5, in the "Otherwise",
	  // in p. 5 of suppl mat. We will update both, and only these
	  // two.
	  to_update = 2; 

#ifdef DEBUGW
	  if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0) {
	    DP2(currentTime);
	    DP2(sp);
	    DP2(popParams[sp].timeLastUpdate);
	    print_spP(popParams[sp]);
	    throw std::out_of_range("currentTime - timeLastUpdate out of range. Serious bug!");
	  }
#endif
	  // if(verbosity >= 2) {
#ifdef DEBUGV
	    Rcpp::Rcout <<"\n     Mutated to existing species " << sp 
			<< " (Genotype = ";
	    print_Genotype(Genotypes[sp]); 
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
	    Rcpp::Rcout << ")"
			<< "\n from species "  <<   nextMutant
			<< " (Genotypes = ";
	    print_Genotype(Genotypes[nextMutant]); 
	      // << "; sp_id = " << Genotypes[sp].to_ullong()
	    Rcpp::Rcout	<< ")";
	    // }
#endif
	  // FIXME00: the if can be removed??
	    // Possibly. But note that the popParams[sp].popSize can be >
	    // 0, but when updated via Algo2 and added to 1.0 we can end
	    // in 1. Why? Because Algo2 can return a 0. The species
	    // "exist" in the sense that it had non-zero pop size when we
	    // last sampled/updated it.
	    
	    // What we do here is step 6 of Algorithm 5, in the
	    // "Otherwise", in p. 5 of suppl mat.

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
	// We increase size by 1, as we already called Algo3. And then
	// update the ti.
	++popParams[nextMutant].popSize;
	to_update = 1;
	u_1 = nextMutant;
	u_2 = -99;
	if(verbosity >= 1)
	  Rcpp::Rcout << "Note: updating in null mutation\n";
      }
    } else { //       *********** We are sampling **********
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
      
      nr_totPopSize_and_fill_out_crude_P(outNS_i, totPopSize, 
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
					 keepEvery,
					 detectionSize,
					 finalTime,
					 //endTimeEvery,
					 detectionDrivers,
					 verbosity,
					 minDetectDrvCloneSz,
					 extraTime,
					 fitnessEffects.drv,
					 cPDetect,
					 PDBaseline,
					 checkSizePEvery,
					 nextCheckSizeP,
					 ran_gen); //keepEvery is for thinning
      if(verbosity >= 3) {
	Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() 
		  << "\n totPopSize after sampling " << totPopSize << "\n";
      }
      
      computeMcFarlandError(e1, n_0, n_1, tps_0, tps_1, 
			    typeModel, totPopSize, K); //, initSize);

      if(simulsDone)
	break; //skip last updateRates

      // if( (typeModel == TypeModel::beerenwinkel) ) {
      // 	updateRatesBeeren(popParams, adjust_fitness_B,
      // 			  initSize, currentTime, alpha, totPopSize,
      // 			  mutationPropGrowth, mu);
      // } else if( (typeModel == TypeModel::mcfarland0) ) {
      // 	updateRatesMcFarland0(popParams, adjust_fitness_MF,
      // 			     K, totPopSize,
      // 			     mutationPropGrowth, mu);
      // } else if( (typeModel == TypeModel::mcfarland) ) {
      // 	updateRatesMcFarland(popParams, adjust_fitness_MF,
      // 			     K, totPopSize);
      // } else if( (typeModel == TypeModel::mcfarlandlog) ) {
      if( (typeModel == TypeModel::mcfarlandlog) ) {
	updateRatesMcFarlandLog(popParams, adjust_fitness_MF,
			     K, totPopSize);
      }
      
#ifdef MIN_RATIO_MUTS_NR
      // could go inside sample_all_pop but here we are sure death, etc, current
      // But I catch them when they are created. Is this really needed?
      for(size_t i = 0; i < popParams.size(); i++) {
	g_tmp1_nr = popParams[i].birth/popParams[i].mutation;
	if(g_tmp1_nr < g_min_birth_mut_ratio_nr) g_min_birth_mut_ratio_nr = g_tmp1_nr;
	
	g_tmp1_nr = popParams[i].death/popParams[i].mutation;
	if(g_tmp1_nr < g_min_death_mut_ratio_nr) g_min_death_mut_ratio_nr = g_tmp1_nr;
      }
#endif
      
      forceSample = false;
    }
  }
}



// [[Rcpp::export]]
Rcpp::List nr_BNB_Algo5(Rcpp::List rFE,
			Rcpp::NumericVector mu_,
			double death,
			double initSize,
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
			Rcpp::IntegerVector initMutant_, 
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
			double checkSizePEvery) {
  // double cPDetect){
  // double n2,
  // double p2,
  // double PDBaseline) {

  
  precissionLoss();
  const std::vector<double> mu = Rcpp::as<std::vector<double> >(mu_);
  const std::vector<int> initMutant = Rcpp::as<std::vector<int> >(initMutant_);
  const TypeModel typeModel = stringToModel(Rcpp::as<std::string>(typeFitness_));

  // A simple, vector-indexed way to map from numeric ids in full to
  // numeric ids in mutator. Recall all genes start with 1. So full2mutator[i-1];
  const std::vector<int> full2mutator = Rcpp::as<std::vector<int> >(full2mutator_);
  // A consistency check

  // const double genTime = 4.0; // should be a parameter. For Bozic only.

  //If seed is -9, then use automatic seed.


  // Code when using randutils
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
    cPDetect = set_cPDetect(n2, p2, PDBaseline);
    if(verbosity >= 1)
      Rcpp::Rcout << "  cPDetect set at " << cPDetect << "\n";
  }
  
  if( (K < 1 ) && ( typeModel ==   TypeModel::mcfarlandlog) )
    throw std::range_error("K < 1.");
  fitnessEffectsAll fitnessEffects =  convertFitnessEffects(rFE);
  //Used at least twice
  std::map<int, std::string> intName = mapGenesIntToNames(fitnessEffects);
  fitness_as_genes genesInFitness = fitnessAsGenes(fitnessEffects);
  PhylogName phylog;

  // Mutator effects
  fitnessEffectsAll muEF;
  if( (full2mutator.size() != 0) ) 
    muEF = convertFitnessEffects(MMUEF);
  else
    muEF = nullFitnessEffects();
  // Paranoia. We should never end up here.
  if( (full2mutator.size() != 0) && (muEF.genomeSize == 0))
    throw std::logic_error("full2mutator > 0 with mutatorEffects.genomesize 0");
  if( (full2mutator.size() == 0) && (muEF.genomeSize != 0)) {
    DP2(muEF.genomeSize);
    throw std::logic_error("full2mutator 0 with mutatorEffects.genomesize != 0");
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

    if(numRuns >= maxNumTries) {
      //  hittedMaxTries This we want here to avoid an extra run and
      //  confusing output
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
      // start of innerBNB now that we do multiple runs if onlyCancer = true.

      nr_innerBNB(
		  fitnessEffects,
		  initSize,
	       K,
		  // alpha,
		  // genTime,
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
		  checkSizePEvery);
      ++numRuns;
      forceRerun = false;
    } catch (rerunExcept &e) {
      Rcpp::Rcout << "\n Recoverable exception " << e.what() 
		  << ". Rerunning.";
      forceRerun = true;
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
    // } else if(numRuns > maxNumTries) {
    //   //  hittedMaxTries FIXME this is very, very confusing in limit
    //   // cases.  suppose maxNumTries = 1. We will run two times, and the
    //   // second might have reached cancer, but we will bail out here, as
    //   // numRuns is actually 2. However, we report the value. And we run
    //   // once more than needed.
    //   hittedMaxTries = true;
    //   Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
    //   runAgain = false;
    //   if(errorHitMaxTries) {
    // 	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
    // 	return
    // 	  List::create(Named("HittedWallTime") = false,
    // 		       Named("HittedMaxTries") = true,
    // 		       Named("other") =
    // 		       List::create(Named("UnrecoverExcept") = false));
    //   }
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
					       returnMFE(e1, //K, 
							 typeModel),
					       Named("errorMF_size") = e1,
					       Named("errorMF_n_0") = n_0,
#ifdef MIN_RATIO_MUTS_NR
					       Named("minDMratio") =
					       g_min_death_mut_ratio_nr,
					       Named("minBMratio") =
					       g_min_birth_mut_ratio_nr,      
#else
					       Named("minDMratio") = -99,
					       Named("minBMratio") = -99,
#endif
					       Named("errorMF_n_1") = n_1,
					       Named("PhylogDF") =  DataFrame::create(
										      Named("parent") = phylog.parent,
										      Named("child") = phylog.child,
										      Named("time") = phylog.time
										      ),
					       Named("UnrecoverExcept") = false)
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
