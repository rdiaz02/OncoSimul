//     Copyright 2013, 2014 Ramon Diaz-Uriarte

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



// FIXME: install and use AMD libraries and recompile R, etc?
// http://blogs.amd.com/developer/2012/04/23/gcc-4-7-is-available-with-support-for-amd-opteron%E2%84%A2-6200-series-and-amd-fx-series-processors/
// Open64: http://devgurus.amd.com/message/1282903#1282903
// With gcc and new flags, I see little or no difference.

// FIXME: openMP?? a few loops, but I call R code for random number
// generator.

#include "OncoSimul.h"
#include <limits>
#include <iostream>
//#include <gsl/gsl_rng.h> // here? in the .h
#include <random>
#include <bitset>
#include <set>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <ctime>
// #include <sys/resource.h> 
#include <sys/time.h> 


// From http://stackoverflow.com/a/5590404
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
       ( std::ostringstream() << std::dec << x ) ).str()
// no longer needed?


#define STOPASSERT(x);

// Will this work under Windows? Probably not. OK, I do not care much.
// void setmemlimit(const long maxram){
//   // try to prevent any overhead if not set
//   if(maxram > 0) {
//     struct rlimit memlimit;
//     long bytes;
    
//     bytes = maxram * (1024*1024);
//     memlimit.rlim_cur = bytes;
//     memlimit.rlim_max = bytes;
//     setrlimit(RLIMIT_AS, &memlimit);
//   }
// }


// #ifdef _WIN32
// void setmemlimit(const long maxram){
// }
// #endif
// #ifndef _WIN32
// // Will this work under Windows? Probably not. OK, I do not care much.
// // No longer used
// void setmemlimit(const long maxram){
//   // try to prevent any overhead if not set
//   if(maxram > 0) {
//     struct rlimit memlimit;
//     long bytes;
    
//     bytes = maxram * (1024*1024);
//     memlimit.rlim_cur = bytes;
//     memlimit.rlim_max = bytes;
//     setrlimit(RLIMIT_AS, &memlimit);
//   }
// }
// #endif


typedef int myT;
typedef std::bitset<64> Genotype64;

struct spParamsP {
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  double timeLastUpdate;
  std::multimap<double, int>::iterator pv;
  double absfitness; 
  int numMutablePos; 
};


static inline void W_f_st(spParamsP& spP){
  spP.W = spP.death + spP.birth + spP.mutation;
}

static inline void R_f_st(spParamsP& spP) {
  spP.R = sqrt( pow( spP.birth - spP.death, 2) + 
		( 2.0 * (spP.birth + spP.death) + 
		  spP.mutation) * spP.mutation );
}

void print_spP(const spParamsP& spP) {
  Rcpp::Rcout <<"\n this is spP\n" 
	    <<"\n popSize = " << spP.popSize
	    <<"\n birth = " << spP.birth
	    <<"\n death = " << spP.death
	    <<"\n W = " << spP.W
	    <<"\n R = " << spP.R
	    <<"\n mutation = " << spP.mutation
	    <<"\n timeLastUpdate = " << spP.timeLastUpdate
	    <<"\n absfitness = " << spP.absfitness
	    <<"\n numMutablePos =" << spP.numMutablePos    
	    <<"\n";
}

static double pM_f_st(const double& t, 
	       const spParamsP& spP){

  long double Ct = cosh(spP.R * t/2.0);
  long double St = sinh(spP.R * t/2.0);
  long double lpM = -99.99;
    
  if( (!std::isfinite(Ct) ) || (!std::isfinite(St)) ) {
    throw std::range_error("pM.f: Ct or St too big");
  }
    

  lpM = (spP.R * Ct + 2.0 * spP.death * St - spP.W * St)/
    (spP.R * Ct - 2.0 * spP.birth * St + spP.W * St);    
  
  double pM = static_cast<double>(lpM);

  if( !std::isfinite(pM) ) {
    print_spP(spP);
    throw std::range_error("pM.f: pM not finite");
  }
  if(pM <= 0.0) {
    print_spP(spP);
    throw std::range_error("pM.f: pM <= 0.0");
  }
  return pM;
}


static inline double pE_f_st(double& pM, const spParamsP& spP){
  double pE = (spP.death * (1.0 - pM ) )/(spP.W - spP.death - spP.birth * pM );
  if( !std::isfinite(pE) ) {
    throw std::range_error("pE.f: pE not finite");
  }
  return pE;
}

static inline double pB_f_st(const double& pE,
			     const spParamsP& spP) {
  return (spP.birth * pE)/spP.death; 
}


static double ti_nextTime_tmax_2_st(const spParamsP& spP,
				    const double& currentTime,
				    const double& tSample,
				    int& ti_dbl_min,
				    int& ti_e3) {
  using namespace Rcpp ;

  double r1;
  double ti;
  double pM;

  // FIXME: should never happen
  if(spP.popSize <= 0.0)
    throw std::range_error("ti: popSize < 0");

  // long double invpop = 1/spP.popSize;
  // long double r;

  double invpop = 1/spP.popSize;
  double r;


  const double epsilon = 10.0;

  if(spP.mutation == 0) { //   spP.W <= -90.0) {
    ti = tSample + 2.0 * epsilon;
    // yes, this is silly but to differentiate from
    // r < pM without further info
    // and to lead to finite value in loop for min.
  } else {
    RNGScope scope;
    r1 = ::Rf_runif(0.0, 1.0);
    r = pow(r1, invpop); 
    pM = pM_f_st(tSample - currentTime, spP);

    if( r < pM) {
      ti = tSample + epsilon;
    } else {
      double tmp2 =  2.0L * spP.mutation;
      double tmp = (spP.birth - spP.death) - spP.mutation;
      double oneminusr = 1.0L - r;
      double numerator =  oneminusr * (tmp + spP.R) + tmp2;
      double denominator = oneminusr * (tmp - spP.R ) + tmp2;
      double invspr = 1.0L/spP.R;
      ti = invspr * log(numerator/denominator);

      if(ti < 0.0) {
	throw std::range_error("ti: eq.11 < 0");
      } 
      if( !std::isfinite(ti) ) {
	throw std::range_error("ti: ti not finite");
      }
      if(ti == 0.0) {
	++ti_dbl_min;
	ti = DBL_MIN;
	// Beware of this!!
	throw std::range_error("ti set to DBL_MIN");
	// Do not exit. Record it. Maybe abort simulation and go to a new one?
	// Rcpp::Rcout << "ti set to DBL_MIN\n";
      }
      if(ti < 0.001) ++ti_e3;
      ti += currentTime;
    } 
  }
  return ti;
}

static double Algo2_st(const spParamsP& spP,
		       const double& ti) {
  using namespace Rcpp ;
  double t = ti - spP.timeLastUpdate;

  if (spP.popSize == 0.0) {
    return 0.0;
  }
  
  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);
  double m; 

  double rnb;
  double retval; 

    if( (1.0 - pe/pm) > 1.0) {
      Rcpp::Rcout << "\n ERROR: Algo 2: (1.0 - pe/pm) > 1.0\n"; 
      throw std::range_error("Algo 2:  1 - pe/pm > 1");
    }

    if( (1.0 - pe/pm) < 0.0 ) {
      throw std::range_error("Algo 2: 1 - pe/pm < 0");
    }

    if( pb > 1.0 ) {
      throw std::range_error("Algo 2: pb > 1 ");
    }

    if( pb < 0.0 ) {
      throw std::range_error("Algo 2: pb < 0");
    }



  if( pe == pm ) {
    // Should never happen. Exact identity??
    Rcpp::Rcout << "\n WARNING: Algo 2: pe == pm \n" ;
    return 0.0;
  }

  
  RNGScope scope;
  m = ::Rf_rbinom(spP.popSize, 1.0 - (pe/pm));
  // this is dangerous. I'd rather throw an exception and bail out soon
  // if(std::isnan(m)) {
  //   // we can get issues with rbinom and odd numbers > 1e15
  //   // see "example-binom-problems.cpp"
  //   // hack this, and issue a warning
  //   Rcpp::Rcout << "\n\nWARNING: Using hack around rbinom NaN problem in Algo2\n";
  //   m = ::Rf_rbinom(spP.popSize + 1, 1.0 - (pe/pm));
  // }
  if(m <= 0.5) { // they are integers, so 0 or 1.
    retval = 0.0;
  } else {
    rnb = ::Rf_rnbinom(m, 1.0 - pb);
    // if(std::isnan(rnb)) {
    //   Rcpp::Rcout << "\n\nWARNING: Using hack around rnbinom NaN problem in Algo2\n";
    //   rnb = ::Rf_rnbinom(m + 1, 1.0 - pb);
    // }
    retval = m + rnb;
  }
 
  if( !std::isfinite(retval) )  {
    throw std::range_error("Algo 2: retval not finite");
  }
  if( std::isnan(retval) )  {
    throw std::range_error("Algo 2: retval is NaN");
  }
  return retval;
}

static double Algo3_st(const spParamsP& spP, const double& t){
  
  using namespace Rcpp ;

  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);

  double m; 
  double retval;
  double rnb;

  if( (1.0 - pe/pm) > 1.0) {
    throw std::range_error("Algo 3:  1 - pe/pm > 1");
  }
  
  if( (1.0 - pe/pm) < 0.0 ) {
    throw std::range_error("Algo 3: 1 - pe/pm < 0");
  }
  
  if( pb > 1.0 ) {
    throw std::range_error("Algo 3: pb > 1 ");
  }
  
  if( pb < 0.0 ) {
    throw std::range_error("Algo 3: pb < 0");
  }
  
  if( pe == pm ) {
    // Should never happen. Exact identity??
    Rcpp::Rcout << "\n WARNING: Algo 3: pm == pe\n"; 
    return 0.0;
  }

  RNGScope scope;
  m = ::Rf_rbinom(spP.popSize - 1.0, 1.0 - (pe/pm));
  // dangerous
  // if(std::isnan(m)) {
  //   // we can get issues with rbinom and odd numbers > 1e15
  //   // see "example-binom-problems.cpp"
  //   // hack this, and issue a warning
  //   Rcpp::Rcout << "\n\nWARNING: Using hack around rbinom NaN problem in Algo3\n";
  //   m = ::Rf_rbinom(spP.popSize, 1.0 - (pe/pm));
  // }
  rnb = ::Rf_rnbinom(m + 2.0, 1.0 - pb);
  // if(std::isnan(rnb)) {
  //   Rcpp::Rcout << "\n\nWARNING: Using hack around rnbinom NaN problem in Algo3\n";
  //   rnb = ::Rf_rnbinom(m + 1.0, 1.0 - pb);
  // }
  retval = m + 1 + rnb;

  if( !std::isfinite(retval) )  {
    throw std::range_error("Algo 3: retval not finite");
  }
  if( !std::isfinite(retval) )  {
    throw std::range_error("Algo 3: retval is NaN");
  }
  return retval;
}

// this is the log of the ratio of death rates
// so the the difference of the successive death rates, as using
// the log version.
static double returnMFE(double& e1,
			const double& K,
			const std::string typeFitness) {
  if(typeFitness == "mcfarlandlog")
    return log(e1);
  else
    return -99;
}

// FIXME But I'd probably want a percent error, compared to the death rate
// something like (log(1+N1/K) - log(1+N2/K))/(log(1+N1/K))
static void computeMcFarlandError(double& e1,
				  double& n_0,
				  double& n_1,
				  double& tps_0,
				  double& tps_1,
				  const std::string typeFitness,
				  const double& totPopSize,
				  const double& K,
				  const double& initSize) {
  // static double tps_0 = initSize;
  // static double tps_1 = 0.0;

  if(typeFitness == "mcfarlandlog" ) {
    double etmp;
    tps_1 = totPopSize;
      if( (tps_0 + 1.0) > tps_1 ) 
	etmp = (K + tps_0 + 1.0)/(K + tps_1);
      else
	etmp = (K + tps_1)/(K + tps_0 + 1);
    if(etmp > e1) {
      e1 = etmp;
      n_0 = tps_0;
      n_1 = tps_1;
    }
    tps_0 = tps_1;
  }
}

static void updateRatesMcFarlandLog(std::vector<spParamsP>& popParams,
				    double& adjust_fitness_MF,
				    const double& K,
				    const double& totPopSize){

  adjust_fitness_MF = log(1.0 + totPopSize/K);

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}

static void fitness(spParamsP& tmpP,
		    const spParamsP& parentP,
		    const int& mutatedPos, 
		    Rcpp::IntegerMatrix restrictTable,
		    const std::string& typeCBN,
		    const Genotype64& newGenotype,
		    const double& birthRate, 
		    const double& s,
		    const double& death,
		    const int& numDrivers,
		    const std::string typeFitness,
		    const double& genTime,
		    const double& adjust_fitness_B,
		    const double& sh,
		    const double& adjust_fitness_MF) {

  using namespace Rcpp;
  int numDependencies;
  int sumDriversMet = 0;
  int sumDriversNoMet = 0;
  int sumDependenciesMet = 0;

  tmpP.birth = parentP.birth;
  tmpP.death = parentP.death;
  tmpP.absfitness = parentP.absfitness;

  if(mutatedPos >= numDrivers) { 
    return;
  } else {
    for(int m = 0; m < numDrivers; ++m) {
      if( newGenotype[m] ) { // this m is mutated
	const Rcpp::IntegerMatrix::Column thisRestrict = 
	  restrictTable(_, m);
	numDependencies = thisRestrict[1];
	if(!numDependencies) { 
	  sumDriversMet++;
	}
	else {
	  sumDependenciesMet = 0;
	  for(int i = 2; i < (2 + numDependencies); i++) {
	    sumDependenciesMet += newGenotype[ thisRestrict[i] ];
	  }
	  if( sumDependenciesMet == numDependencies ) {
	    sumDriversMet++;   
	  } else {
	    sumDriversNoMet++;
	  }
	}
      }
    }
  }

  if((sh < 0) && sumDriversNoMet) {
    tmpP.absfitness = 0.0;
    tmpP.death = 1.0;
    tmpP.birth = 0.0; 
  } else {
    if(typeFitness == "bozic1") {
      tmpP.death = pow( 1.0 - s, sumDriversMet) * 
	pow( 1.0 + sh, sumDriversNoMet);
      tmpP.birth = 1.0;
    } else if(typeFitness == "mcfarlandlog") {
      tmpP.birth = pow(1.0 + s, sumDriversMet) / 
	pow( 1.0 + sh, sumDriversNoMet);
    } else if (typeFitness == "exp") { 
      tmpP.birth = pow(1.0 + s, sumDriversMet) * 
	pow( 1.0 - sh, sumDriversNoMet);
    } else {
      throw std::range_error("Eh?!! What are you doing here. Pass a recognized typeFitness.");

    } 
  }
}

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
}

static inline void mapTimes_updateP(std::multimap<double, int>& mapTimes,
			     std::vector<spParamsP>& popParams,
			     const int index,
			     const double time) {
  if(popParams[index].timeLastUpdate > -1)
    mapTimes.erase(popParams[index].pv);
  popParams[index].pv = mapTimes.insert(std::make_pair(time, index));
}


static inline void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes) {
  nextMutant = mapTimes.begin()->second;
  minNextMutationTime = mapTimes.begin()->first;
}

static void remove_zero_sp_v7(std::vector<int>& sp_to_remove,
			      std::vector<Genotype64>& Genotypes,
			      std::vector<spParamsP>& popParams,
			      std::multimap<double, int>& mapTimes) {
  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();
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
					    int& lastMaxDr,
					    double& done_at,
					    const std::vector<Genotype64>& Genotypes,
					    const std::vector<spParamsP>& popParams, 
					    const double& currentTime,
					    const int& NumDrivers,
					    const double& keepEvery,
					    const double& detectionSize,
					    const double& finalTime,
					    const double& endTimeEvery,
					    const int& finalDrivers,
					    const int& verbosity) {
  bool storeThis = false;
  totPopSize = 0.0;
  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
  }

  if( currentTime >= (lastStoredSample + keepEvery) ) {
    storeThis = true;
  }

  if( (totPopSize <= 0.0) || (currentTime >= finalTime)  ) {
    simulsDone = true;
  }

  if(endTimeEvery > 0) {
    if(done_at <= 0 ) {
      if( (totPopSize >= detectionSize) ||
	   (lastMaxDr >= finalDrivers)  )
	done_at = currentTime + endTimeEvery;
    } else if (currentTime >= done_at) {
      if( (totPopSize >= detectionSize) ||
	  (lastMaxDr >= finalDrivers)  )
	simulsDone = true;
      else
	done_at = -9;
    }
  } else if( (totPopSize >= detectionSize) ||
	     (lastMaxDr >= finalDrivers) )  {	
      simulsDone = true;
  }

  if(simulsDone)
    storeThis = true;

  
  if( storeThis ) {
    lastStoredSample = currentTime;
    outNS_i++;
    int tmp_ndr = 0;
    int max_ndr = 0;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    
    time_out.push_back(currentTime);
    
    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);
      
      tmp_ndr = count_NDrivers(Genotypes[i], NumDrivers);
      if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = tmp_ndr;
      }
      lastMaxDr = max_ndr;
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
  
  if(totPopSize > (4.0 * 1e15)) {
    if(verbosity > 0)
      Rcpp::Rcout << "\nNOTE: popSize > 4e15. Likely loss of precission\n";
  }
}

static inline void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
			       const std::vector<double>& sampleTotPopSize,
			       const std::vector<double>& sampleLargestPopSize,
			       const std::vector<double>& sampleLargestPopProp,
			       const std::vector<int>& sampleMaxNDr,
			       const std::vector<int>& sampleNDrLargestPop){

  for(size_t i = 0; i < sampleTotPopSize.size(); ++i) {
    perSampleStats(i, 0) = sampleTotPopSize[i];
    perSampleStats(i, 1) = sampleLargestPopSize[i];
    perSampleStats(i, 2) = sampleLargestPopProp[i];
    perSampleStats(i, 3) = static_cast<double>(sampleMaxNDr[i]);
    perSampleStats(i, 4) = static_cast<double>(sampleNDrLargestPop[i]);
  }
}

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
  maxNumDrivers = 0;
  int tmpdr = 0;
  
  for(int j = 0; j < returnGenotypes.ncol(); ++j) {
    tmpdr = 0;
    for(int i = 0; i < numDrivers; ++i) {
      tmpdr += returnGenotypes(i, j);
      countByDriver[i] += returnGenotypes(i, j);
    }
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
			     const std::vector<Genotype64>& Genotypes,
			     const double& tSample){

  sp_to_remove.clear();

  for(size_t i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize = 
	Algo2_st(popParams[i], tSample);
    }
    if( popParams[i].popSize <=  0.0 ) {
      sp_to_remove.push_back(i);
    } 
  }
}


static void precissionLoss(){
  // We are storing population sizes as doubles.
  // Should not loose any precission up to 2^53 - 1
  // (e.g., http://stackoverflow.com/a/1848762)
  // but double check if optims break it.

  // Problems are likely to arise sooner, with 4.5e15, because
  // of rbinom. See notes in example-binom-problems.cpp
  // We warn about that if totPopSize > 4e15
  double a, b, c, d;
  int e, f;
  a = pow(2, 52) + 1.0;	
  b = pow(2, 52); // 2^53 a little over 9*1e15
  c = (9.0 * 1e15) + 1.0;
  d = (9.0 * 1e15);

  e = static_cast<int>(a - b);
  f = static_cast<int>(c - d);

  if( a == b) Rcpp::Rcout << "WARNING!!!! \n Precission loss: a == b\n";
  if( !(a > b)) Rcpp::Rcout << "WARNING!!!! \n Precission loss: !(a > b)\n";
  if(c == d) Rcpp::Rcout << "WARNING!!!! \n Precission loss: c == d\n";
  if( !(c > d)) Rcpp::Rcout << "WARNING!!!! \n Precission loss: !(c > d)\n";
  if( e != 1 ) Rcpp::Rcout << "WARNING!!!! \n Precission loss: e != 1\n";
  if( f != 1 ) Rcpp::Rcout << "WARNING!!!! \n Precission loss: f != 1\n";
}

static void init_tmpP(spParamsP& tmpParam) {
  tmpParam.popSize = -std::numeric_limits<double>::infinity();
  tmpParam.birth = -std::numeric_limits<double>::infinity();
  tmpParam.death = -std::numeric_limits<double>::infinity();
  tmpParam.W = -std::numeric_limits<double>::infinity();
  tmpParam.R = -std::numeric_limits<double>::infinity();
  tmpParam.mutation = -std::numeric_limits<double>::infinity();
  tmpParam.timeLastUpdate = std::numeric_limits<double>::infinity();
  tmpParam.absfitness = -std::numeric_limits<double>::infinity();
  tmpParam.numMutablePos = -999999;
}

// [[register]]
SEXP Algorithm5(SEXP restrictTable_,
		SEXP numDrivers_,
		SEXP numGenes_,
		SEXP typeCBN_,
		SEXP birthRate_, 
		SEXP s_, 
		SEXP death_,
		SEXP mu_,
		SEXP initSize_,
		SEXP sampleEvery_,
		SEXP detectionSize_,
		SEXP finalTime_,
		SEXP initSize_species_,
		SEXP initSize_iter_,
		SEXP seed_gsl_,
		SEXP verbose_,
		SEXP speciesFS_,
		SEXP ratioForce_,
		SEXP typeFitness_,
		SEXP maxram_,
		SEXP mutatorGenotype_,
		SEXP initMutant_,
		SEXP maxWallTime_,
		SEXP keepEvery_,
		SEXP alpha_,
		SEXP sh_,
		SEXP K_,
		SEXP endTimeEvery_,
		SEXP finalDrivers_) {
  
  BEGIN_RCPP
    
    using namespace Rcpp;

  precissionLoss();
  
  const IntegerMatrix restrictTable(restrictTable_);
  const int numDrivers = as<int>(numDrivers_);
  const int numGenes = as<int>(numGenes_);
  std::string typeCBN = as<std::string>(typeCBN_); typeCBN = "CBN";
  const std::string typeFitness = as<std::string>(typeFitness_);
  // birth and death are irrelevant with Bozic
  const double birthRate = as<double>(birthRate_);
  const double death = as<double>(death_);
  const double s = as<double>(s_);
  const double mu = as<double>(mu_);
  const double initSize = as<double>(initSize_);
  const double sampleEvery = as<double>(sampleEvery_);
  const double detectionSize = as<double>(detectionSize_);
  const double finalTime = as<double>(finalTime_);
  const int initSp = as<int>(initSize_species_);
  const int initIt = as<int>(initSize_iter_); // FIXME: this is a misnomer
  int verbosity = as<int>(verbose_); // ++verbosity; // will use later. for now, shut up the unused variable warnings
  // const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); 
  int speciesFS = as<int>(speciesFS_);
  const int seed = as<int>(seed_gsl_);
  const long maxram = as<int>(maxram_);
  const int mutatorGenotype = as<int>(mutatorGenotype_);
  int initMutant = as<int>(initMutant_); initMutant++;
  const double maxWallTime = as<double>(maxWallTime_);
  const double keepEvery = as<double>(keepEvery_);
  double alpha = as<double>(alpha_); ++alpha;
  const double sh = as<double>(sh_); 
  const double K = as<double>(K_); 
  const double endTimeEvery = as<double>(endTimeEvery_); 
  const int finalDrivers = as<int>(finalDrivers_); 


  double lastStoredSample;
  const double genTime = 4.0; // should be a parameter. For Bozic only.

  // Nope, do not use as you need to exit R to restore
// #ifndef _WIN32  
//   if(maxram)  setmemlimit(maxram);
// #endif

  time_t start_time = time(NULL);
  double runningWallTime = 0;

  bool forceSample = false;
  bool simulsDone = false;

  double totPopSize = 0;
  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double currentTime = 0.0;
  double timeNextPopSample;
  double tSample;

  int nextMutant;
  unsigned int numSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int outNS_i = 0; 
  unsigned int sp = 0;
  

  int iterL = 1000; iterL++;
  int speciesL = 1000; speciesL++;
  //int timeL = 1000;
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  std::vector<int>sp_to_remove(1);
  sp_to_remove.reserve(10000);
  

  // some checks. Do this systematically
  // FIXME: do only if mcfarland!
  if(K < 1 )
    throw std::range_error("K < 1.");

  // verify we are OK with usigned long long
  if( !(static_cast<double>(std::numeric_limits<unsigned long long>::max()) 
  	>= pow(2, 64)) )
    throw std::range_error("The size of unsigned long long is too short.");

  if(numGenes > 64)  
    throw std::range_error("This version only accepts up to 64 genes.");

  int to_update = 1; 
  int u_1 = -99;
  int u_2 = -99;

  //GSL rng
  // gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  // gsl_rng_set (r, (unsigned long) seed);

  // C++11 random number
  std::mt19937 ran_generator(seed);



  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  //  std::set<Genotype64> uniqueGenotypes;
  std::set<unsigned long long> uniqueGenotypes;
  spParamsP tmpParam; 
  std::vector<spParamsP> popParams(1);
  const int sp_per_period = 5000;

  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

  std::vector<int>mutablePos(numGenes); 

  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; 

  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  std::multimap<double, int> mapTimes;

  int ti_dbl_min = 0;
  int ti_e3 = 0;

  int maxNumDrivers = 0;
  int totalPresentDrivers = 0;
  std::vector<int>countByDriver(numDrivers, 0);
  std::string occurringDrivers;
  
  std::vector<double> sampleTotPopSize;
  std::vector<double> sampleLargestPopSize;
  std::vector<int> sampleMaxNDr; 
  std::vector<int> sampleNDrLargestPop;
  sampleTotPopSize.reserve(initIt);
  sampleLargestPopSize.reserve(initIt);
  sampleMaxNDr.reserve(initIt);
  sampleNDrLargestPop.reserve(initIt);

  double adjust_fitness_B = -std::numeric_limits<double>::infinity();
  double adjust_fitness_MF = -std::numeric_limits<double>::infinity();

  double e1, n_0, n_1, tps_0, tps_1; 
  tps_0 = 0.0;
  tps_1 = 0.0;
  e1 = 0.0;
  n_0 = 0.0;
  n_1 = 0.0;

  int lastMaxDr = 0;
  double done_at = -9;

  init_tmpP(tmpParam);
  init_tmpP(popParams[0]);

  lastStoredSample = 0.0;
  Genotypes[0].reset();
  popParams[0].popSize = initSize;
  totPopSize = initSize;

  tps_0 = totPopSize;
  e1 = 0.0;
  tps_1 = totPopSize;
  
  popParams[0].numMutablePos = numGenes;
  if(typeFitness == "mcfarlandlog") {
    popParams[0].birth = 1.0;
    popParams[0].death = log(1.0 + totPopSize/K);
  } else if(typeFitness == "bozic1") {
    popParams[0].birth = 1.0;
    popParams[0].death = 1.0;
  } else if (typeFitness == "exp") {
    popParams[0].birth = 1.0;
    popParams[0].death = death;
  } else { 
    throw std::range_error("Eh?!! What are you doing here. Pass a recognized typeFitness.");
  }

  if(mutatorGenotype)
    popParams[0].mutation = mu * popParams[0].birth * popParams[0].numMutablePos;
  else
    popParams[0].mutation = mu * popParams[0].numMutablePos;

  W_f_st(popParams[0]);
  R_f_st(popParams[0]);

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  genot_out.push_back(Genotypes[0]);
  popSizes_out.push_back(popParams[0].popSize);
  index_out.push_back(outNS_i);
  uniqueGenotypes.insert(Genotypes[0].to_ullong());
  time_out.push_back(currentTime);

  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;

  sampleTotPopSize.push_back(popParams[0].popSize);
  sampleLargestPopSize.push_back(popParams[0].popSize);
  sampleMaxNDr.push_back(count_NDrivers(Genotypes[0], numDrivers));
  sampleNDrLargestPop.push_back(sampleMaxNDr[0]);


  while(!simulsDone) {
    iter++;
    tSample = std::min(timeNextPopSample, finalTime);

    if(iter == 1) { 
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample, 
					 ti_dbl_min, ti_e3);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      popParams[0].timeLastUpdate = currentTime;
    } else {
      if(to_update == 1) {
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample, 
					   ti_dbl_min, ti_e3);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);
	popParams[u_1].timeLastUpdate = currentTime;

      } else if(to_update == 2) {
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

      } else { 
	for(size_t i = 0; i < popParams.size(); i++) {
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample, ti_dbl_min, ti_e3);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].timeLastUpdate = currentTime;
	  
	}
      }
    }
    if(forceSample) {
      tSample = currentTime;
      timeNextPopSample = currentTime;
    } 
    
    
    getMinNextMutationTime4(nextMutant, minNextMutationTime, 
			    mapTimes);
  
    
    if( minNextMutationTime <= tSample ) {
      currentTime = minNextMutationTime;

      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
      }      
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2;
      }
      runningWallTime = difftime(time(NULL), start_time);
      if( runningWallTime > maxWallTime ) {
	forceSample = true;
	simulsDone = true;
      }
      getMutatedPos_bitset(mutatedPos, numMutablePosParent, //r,
			   ran_generator,
			   mutablePos,
			   Genotypes[nextMutant], 
			   numGenes);

      newGenotype = Genotypes[nextMutant];
      newGenotype.set(mutatedPos);
      
      new_sp_bitset(sp, newGenotype, Genotypes);

      if(sp == numSpecies) {
	++numSpecies;
	init_tmpP(tmpParam);
	tmpParam.popSize = 1;

	fitness(tmpParam, popParams[nextMutant], mutatedPos, 
		restrictTable,
		typeCBN, newGenotype, birthRate, s, death,
		numDrivers, typeFitness, genTime,
		adjust_fitness_B, sh, adjust_fitness_MF);
	

	if(tmpParam.birth > 0.0) {
	  tmpParam.numMutablePos = numMutablePosParent - 1;
	  if(mutatorGenotype)
	    tmpParam.mutation = mu * tmpParam.birth * tmpParam.numMutablePos;
	  else
	    tmpParam.mutation = mu * tmpParam.numMutablePos;
	  if (tmpParam.mutation > 1 )
	    Rcpp::Rcout << "WARNING: mutation > 1\n";
	  if (numMutablePosParent == 1) 
	    Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
	  W_f_st(tmpParam);
	  R_f_st(tmpParam);
	  tmpParam.timeLastUpdate = -99999.99999; 
	  popParams.push_back(tmpParam);
	  Genotypes.push_back(newGenotype);
	  to_update = 2;
	} else {
	  --sp;
	  --numSpecies;
	  to_update = 1;
	}
      } else {	
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2_st(popParams[sp], currentTime);
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	}
      }
      u_1 = nextMutant;
      u_2 = static_cast<int>(sp);
    } else { 
      to_update = 3; 
      currentTime = tSample;
      sample_all_pop_P(sp_to_remove, 
		       popParams, Genotypes, tSample);
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
				      lastMaxDr,
				      done_at,
				      Genotypes, popParams, 
				      currentTime,
				      numDrivers,
				      keepEvery,
				      detectionSize,
				      finalTime,
				      endTimeEvery,
				      finalDrivers,
				      verbosity); 
      computeMcFarlandError(e1, n_0, n_1, tps_0, tps_1, 
			    typeFitness, totPopSize, K, initSize);
    
      if(simulsDone)
	break; 

      if( (typeFitness == "mcfarlandlog") ) {
	updateRatesMcFarlandLog(popParams, adjust_fitness_MF,
				K, totPopSize);
      }
      forceSample = false;
    }
  }
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species


  // FIXME: all this is ugly and could be a single function
  // up to call to IntegerMatrix
  std::vector<unsigned long long> genot_out_ullong(genot_out.size());
  genot_out_to_ullong(genot_out_ullong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ullong);
  std::vector<unsigned long long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);
  IntegerMatrix returnGenotypes(numGenes, uniqueGenotypes_vector.size());
  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);

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
  NumericMatrix outNS(outNS_r, outNS_c);  
  if(create_outNS) {
    reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ullong, 
		     popSizes_out, 
		     index_out, time_out);
    
  } else {
    outNS(0, 0) = -99;
  }

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

  NumericMatrix perSampleStats(outNS_i + 1, 5);
  fill_SStats(perSampleStats, sampleTotPopSize, sampleLargestPopSize,
	      sampleLargestPopProp, sampleMaxNDr, sampleNDrLargestPop);
  
  return 
    List::create(Named("pops.by.time") = outNS,
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
		 // Named("outi") = outNS_i + 1, // silly. Use the real number of samples. FIXME
		 Named("HittedWallTime") = (runningWallTime > maxWallTime),
		 // Named("iRunningWallTime") = runningWallTime,
		 // Named("oRunningWallTime") = difftime(time(NULL), start_time),
		 // Named("ti_dbl_min") = ti_dbl_min,
		 // Named("ti_e3") = ti_e3,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 Named("OccurringDrivers") = occurringDrivers,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = List::create(Named("errorMF") = 
					       returnMFE(e1, K, 
							 typeFitness),
					       Named("errorMF_size") = e1,
					       Named("errorMF_n_0") = n_0,
					       Named("errorMF_n_1") = n_1)
		 );

  END_RCPP
    
    }







