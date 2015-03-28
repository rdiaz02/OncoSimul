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



// The new semantics: add monotonicity or not. Will think how.

// Restrctions:
//  1, 2, 3, 4, - 6 - monot

// or

//  1:2:3:4,6,monot

// where the 1:2:3:4 are a module, so if any is hit, all is hit. Allows
// combining modules with/without monotonicity.

// FIXME: pass restrictTable by reference?
// FIXME: do not use "mutatorGenotype" but "mutProp" or similar.

// I will use Bozic1, but do not start from init mutant
// as Tomassetti, 2013 paper indicates lots of passengers
// before tumor.

// Choose a set of fast parameters.
// Start from init size larger than 1, something like 100.


// Add again running wall time

// Discrete time?
// Arbitrary genomes?


// Can I reproduce Bozic's numerical results?




// Use simple trees with few restrictions: shallow ones
//      Could it be that oncog trees capture things that are not there?
//        i.e., they impose too complex structure?

// Could it be that there is no restriction, but coadyuvation?
//  something like Helm's model?

// Could mixtures be the results of little sampling, of incomplete data?


// Use model in Mather?

// Check rates are OK in terms of mutation?


// About exceptions in binom, ti_dbl_min, etc:
// easier to trigger if compiled with -O2 instead of -Ofast;
// still, easier in Xeon than in the Opteron.



// enums and switch: http://www.codeguru.com/cpp/cpp/cpp_mfc/article.php/c4067/Switch-on-Strings-in-C.htm
// just plain enums



// check mutation rates in bozic against new sets of births and deaths

// death in bozic2/(60 * mut), 20 drivers top
//            0.1   0.02
//   1e-5     25    140
//   5e-6     50    278
//   1e-6     250

//  if bozic1
//            0.1   0.02
//   1e-5     200   1100
//   5e-6     400   2200



// Check I can reproduce Bozic's patterns
// Show, in terms of drivers, drivers at end, etc, same as mine
//     with linear model
// Use linear model

// - Something of comparing with driver/passenger stuff in Bozic?
// - Launch a bunch of bozic's that will leave output: more than one driver
// - Same as above, but with rt0.60 (60 drivers): log, linear, linear mutator
// - Linear no mutator with rt0
// - Linear mutator with rt0
// - Log with rt0
// - Check patterns


// - check we reproduce patterns: launch and save and plot (use sampleEvery of 5)
// Do it with the four trees and prepare "for real" scripts.



// - Prepare launching script; maybe no mc.preschedule.
// - Check runs fit in patterns of expansion, waves, etc, with Bozic and Beeren
// Output from real runs: just output info, no redundant stuff

// - Run gcc with profiling for those and optimize
// - Launch

// Why not use number of genes? Because what matters are the drivers. If I
// get enough passengers, I can always cut down after the simulations. And
// if passengers (actually, total) number genes large enough, unlikely to
// get into "no possitions left to mutate", (Although simuls faster if
// fewer positions, of course, for a certain definition of "fast".)


// Why not start from a driver, as Bozic do? Because then we are not
// allowing for that driver to have been prceeded by a passenger.


// Mutation: should we leave it unchanged?
// NO: that only makes sense if we allow back mutation. If we do not,
// that creates another problem: we mutation from a species to itself
// and this is a mess. What do we update? When? How?


// - Rethink no back mutation: what are sprouffske, bozic, and beerenwinkel doing?

// Rates of birth and death in no mutator:

//            Birth                 Death
// Linear   whatever I write     
// Bozic     [0.5, 1)              (0, 0.5]
// Exp        >= 1                1 (I set it to this in R, not C++)


//   Rates of mutation as funct num genes
//                40               60
// 2e-4         .008             .012
// 1e-4         .004             .006
// 2e-5         .0008            .0012
// 1e-5         .0004            .0006
// 1e-6         .00004           .00006


// Clean up code in general:
//  - recheck time used by R
//  - minimize debugging and verbose code, or wrap in DEBUG statements
//  - remove all unused code from former functions
//  - use "static"??
//  - profile again, using all current fitness functions
//  - compile after profiling with gcc, after deciding function to use

// Versions:
// P is better than O in every case (and both always better or much, much
// [> 10x]) better than M, K, etc. So use P.  Sampling 10 seems generally
// good, although occasionally a very slow result.  20 is generally
// slow. 5 is about as 10, but larger object sizes.  2 can be about as 5,
// or slightly worse, with much larger object sizes generally.  1 is often
// slow and very large objects.


// FIXME: install and use AMD libraries and recompile R, etc?
// http://blogs.amd.com/developer/2012/04/23/gcc-4-7-is-available-with-support-for-amd-opteron%E2%84%A2-6200-series-and-amd-fx-series-processors/
// Open64: http://devgurus.amd.com/message/1282903#1282903
// With gcc and new flags, I see little or no difference.

// FIXME: keep track of ancestor: A string (we keep concatenating ancestors)
//        added to structure? Pass to R as a different object, not 
//        returnParams

// FIXME: use enums for typeFitness and typeCBN. At least check values are
// within allowed ones. At start of Algo5.


// FIXME: openMP?? The loop at 5.2 and at 5.9 
// but I call R code for random number generator.


#include "OncoSimul.h"
#include <limits>
#include <iostream>
// #include <gsl/gsl_rng.h> // here? in the .h
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

// #include <exception>
#include <stdexcept>



//#define DEBUGZ
// #define DEBUGV
//#define DEBUGW

#ifdef DEBUGW
#define ASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << "\n\n";		\
    }									\
  }
#else
#define ASSERT(x);
//#define ARMA_NO_DEBUG
#endif


#ifdef DEBUGW
#define STOPASSERT(x) {							\
    if (! (x)) {							\
      Rcpp::Rcout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	Rcpp::Rcout << " on line " << __LINE__  << std::endl;		\
	throw std::out_of_range("STOPASSERT");				\
    }									\
  }
#else
#define STOPASSERT(x);
#endif


#define DP2(x) {Rcpp::Rcout << "\n DEBUG2: Value of " << #x << " = " << x << std::endl;}
//#define DP2(x);

#ifdef DEBUGV
#define DP(x) {Rcpp::Rcout << "\n DEBUG: Value of " << #x << " = " << x << std::endl;}
#else
#define DP(x);
#endif

// To track if mutation is really much smaller than birth/death
#define MIN_RATIO_MUTS
#ifdef MIN_RATIO_MUTS
// There is really no need for these to be globals?
// Unless I wanted to use them inside some function. So leave as globals.
double g_min_birth_mut_ratio = DBL_MAX;
double g_min_death_mut_ratio = DBL_MAX;
double g_tmp1 = DBL_MAX;
#endif

// Memory limits, from
// https://marylou.byu.edu/wiki/How+do+I+limit+the+amount+of+memory+that+my+program+uses%3F
// In Linux, and if we are sys admins, I'd rather do this via disabling
// memory overcommiting
// (e.g. http://stackoverflow.com/questions/12582793/limiting-memory-usage-in-r-under-linux,
// and this explanation http://www.win.tue.nl/~aeb/linux/lk/lk-9.html and
// this other thread
// http://serverfault.com/questions/141988/avoid-linux-out-of-memory-application-teardown
// )
// But that might not be possible. Note that ulimit unlikely to work, as
// we can run under multicore, etc (but see the next two, for possible approaches
// which I have not tested 
//  http://coldattic.info/shvedsky/pro/blogs/a-foo-walks-into-a-bar/posts/40
//  https://github.com/pshved/timeout).

// So I allow to limit the RAM. Set to 0 for no effect whatsoever.
// Apparently little cost:
// compare 
// test-speed-max-ram-with-limit.R test-speed-max-ram.R

// Will this work under Windows? Probably not.
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



// Simple custom exception for exceptions that lead to re-runs.
class rerunExcept: public std::runtime_error {
public:
  rerunExcept(const std::string &s) :
    std::runtime_error(s) {}
};


void here(std::string x) {
  Rcpp::Rcout << "\n DEBUG: HERE at " << x << std::endl;
}

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
  double absfitness; //convenient for Beerenwinkel
  int numMutablePos; //for mutator if need update of mutation
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
    
  // my expression, which I think is better       
  //  lpM = (R * Ct + St * (2.0 * death - W ))/(R * Ct + St * (W - 2.0 * growth));
  // theirs, in paper and code
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

// FIXME: if this is expenssive, use doubles, not long doubles.
// after all, we might end up with DBL_MIN anyway if large
// pop size (> 10^9) and runif very close to 1.

// Some comparisons with the exponential (beeren) fitness
// using long doubled showed problems in ~ 25% of cases
// (42 warnings total out of 179 completed C++ executions).
// Using just double, I get 63 out of 205, so ~ 30%.
static double ti_nextTime_tmax_2_st(const spParamsP& spP,
				    const double& currentTime,
				    const double& tSample,
				    int& ti_dbl_min,
				    int& ti_e3) {
  // Following the logic of the code by Mather in
  // findNextMutationTime

  // We return the nextMutationTime or a value larger than the
  // max length of the period (tSample)

  // I also change names rr, r, to match those in Mather r1, r.
  
  // However, I pass mutation, and split computation to avoid numerical problems
  // I was getting ti == 0 and ti < 0 in the other versions with large N.
  using namespace Rcpp ;

  double r1;
  double ti;
  double pM;

  // FIXME: should never happen
  if(spP.popSize <= 0.0) {
    throw std::range_error("ti: popSize <= 0. spP.popSize = "
			   + std::to_string(spP.popSize));
  }
  // long double invpop = 1/spP.popSize;
  // long double r;

  double invpop = 1/spP.popSize;
  double r;


  const double epsilon = 10.0;

  // W < 0 is a signal that mutation is zero, and thus ti is Inf
  if(spP.mutation == 0) { //   spP.W <= -90.0) {
    ti = tSample + 2.0 * epsilon;
    // yes, this is silly but to differentiate from
    // r < pM without further info
    // and to lead to finite value in loop for min.
    //ti = std::numeric_limits<double>::infinity();
  } else {

    RNGScope scope;
    r1 = ::Rf_runif(0.0, 1.0);
    // this was in the original Mather code, but I doubt
    // it really makes it more stable, and seems more expensive
    // r = exp((1.0 / n) * log(r1));
    // r = pow(r1, 1.0/spP.popSize); //what I do
    // r = exp( invpop * log(static_cast<long double>(r1))  );
    //r = pow(static_cast<long double>(r1), invpop); // what I do
    r = pow(r1, invpop); // what I do
    pM = pM_f_st(tSample - currentTime, spP);

    if( r < pM) {// time to mutation longer that this time period
      ti = tSample + epsilon;
    } else {
      // Expand numerator and denominatior, the term for W and simplify.
      // Then, express as (1- r) (which is, inclussively, between 0 and 1)
      // and then multiply by -1 to take the log of each
      
      // long double tmp2 =  2.0L * spP.mutation;
      // long double tmp = (spP.birth - spP.death) - spP.mutation;
      // long double oneminusr = 1.0L - r;
      // long double numerator =  oneminusr * (tmp + spP.R) + tmp2;
      // long double denominator = oneminusr * (tmp - spP.R ) + tmp2;

      double tmp2 =  2.0L * spP.mutation;
      double tmp = (spP.birth - spP.death) - spP.mutation;
      double oneminusr = 1.0L - r;
      double numerator =  oneminusr * (tmp + spP.R) + tmp2;
      double denominator = oneminusr * (tmp - spP.R ) + tmp2;

      // numerator =  (1.0 - r) * (spP.R + spP.birth - spP.death - spP.mutation) 
      // 	+ 2.0 * spP.mutation;
      // denominator = (1.0 - r) * (spP.birth - spP.death - spP.mutation - spP.R ) 
      // 	+ 2.0 * spP.mutation;

      // FIXME? is it really necessary to use log(-a) - log(-b) or could
      // I just use log(a/b), where a and b are -numerator and -denominator?
      // use the log of ratio, in case negative signs in numerator or denom.
      // long double invspr = 1.0L/spP.R;
      // ti = static_cast<double>(invspr * (log(numerator) - log(denominator)));
      double invspr = 1.0L/spP.R;
      ti = invspr * log(numerator/denominator);


      //ti = (1.0/spP.R) * (log(numerator) - log(denominator));
	
      //eq. 11
      // ti = (1.0/R) * (log( -1 * (r * (R - W + 2.0 * growth) - W - R + 2.0 * death )) -
      //  		      log( -1 * (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death )));

      // ti = (1.0/R) * log( (r * (R - W + 2.0 * growth) - W - R + 2.0 * death) /
      //               (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death));
      // Rcpp::Rcout << "\n this is ti = " << ti << "\n";
      if(ti < 0.0) {
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);

	Rcpp::Rcout << "\n ERROR: ti: eq.11 < 0 \n";
	// Rcpp::Rcout << "\n R = " << R;
	// Rcpp::Rcout << "\n W = " << W;
	// Rcpp::Rcout << "\n r1 = " << r1;
	// Rcpp::Rcout << "\n r = " << r;
	// Rcpp::Rcout << "\n n = " << n;
	// Rcpp::Rcout << "\n mu = " << mu;
	// Rcpp::Rcout << "\n growth = " << growth;
	// Rcpp::Rcout << "\n death = " << death << "\n";
	Rcpp::Rcout << "\n numerator = " << numerator;
	Rcpp::Rcout << "\n denominator = " << denominator;
	Rcpp::Rcout << "\n is r > 1? " << (r > 1.0) << "\n";
	Rcpp::Rcout << "\n is r < 0? " << (r < 0.0) << "\n";
	Rcpp::Rcout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	throw std::range_error("ti: eq.11 < 0");
      } 
      if( !std::isfinite(ti) ) {
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);
	double numerator2 = r * (spP.R - spP.W + 2.0 * spP.birth) - 
	  spP.W - spP.R + 2.0 * spP.death;
	double denominator2 = r * (-spP.R - spP.W + 2.0 * spP.birth) - 
	  spP.W + spP.R + 2.0 * spP.death;
	double ti2 = invspr * log(numerator2/denominator2);

	if(std::abs(ti - ti2) > 1e-5) {
	  DP2(ti);
	  DP2(ti2);
	  DP2(numerator);
	  DP2(numerator2);
	  DP2(denominator);
	  DP2(denominator2);
	  DP2(invspr);
	  DP2(r);
	  DP2(r1);
	  DP2(tmp);
	  DP2(tmp2);
	  DP2(oneminusr);
	  //print_spP(spP);
	}
	//	Rcpp::Rcout << "\n ERROR: ti not finite \n";
	// Rcpp::Rcout << "\n R = " << R;
	// Rcpp::Rcout << "\n W = " << W;
	Rcpp::Rcout << "\n r1 = " << r1;
	Rcpp::Rcout << "\n r = " << r;
	// Rcpp::Rcout << "\n n = " << n;
	// Rcpp::Rcout << "\n growth = " << growth;
	// Rcpp::Rcout << "\n death = " << death << "\n";
	Rcpp::Rcout << "\n numerator = " << numerator;
	Rcpp::Rcout << "\n denominator = " << denominator;
	Rcpp::Rcout << "\n ti2 = " << ti2;
	Rcpp::Rcout << "\n numerator2 = " << numerator2;
	Rcpp::Rcout << "\n denominator2 = " << denominator2;

	Rcpp::Rcout << "\n is r > 1? " << (r > 1.0) << "\n";
	Rcpp::Rcout << "\n is r < 0? " << (r < 0.0) << "\n";
	Rcpp::Rcout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	Rcpp::Rcout << "\n tmp = " << tmp << "\n";
	Rcpp::Rcout << "\n tmp2 = " << tmp2 << "\n";
	Rcpp::Rcout << "\n eq12 = " << eq12 << "\n";
	print_spP(spP);
	throw std::range_error("ti: ti not finite");
      }
      if(ti == 0.0) {
#ifdef DEBUGV	
	// FIXME: pass verbosity as argument, and return the warning
	// if set to more than 0?
	Rcpp::Rcout << "\n\n\n WARNING: ti == 0. Setting it to DBL_MIN \n"; 
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);

	Rcpp::Rcout << "\n tmp2 = " << tmp2;
	Rcpp::Rcout << "\n tmp = " << tmp;
	Rcpp::Rcout << "\n invspr = " << invspr;
	Rcpp::Rcout << "\n invpop = " << invpop;
	// Rcpp::Rcout << "\n R = " << R;
	// Rcpp::Rcout << "\n W = " << W;
	// Rcpp::Rcout << "\n r1 = " << r1;
	// Rcpp::Rcout << "\n r = " << r;
	// Rcpp::Rcout << "\n n = " << n;
	// Rcpp::Rcout << "\n growth = " << growth;
	// Rcpp::Rcout << "\n death = " << death << "\n";
	Rcpp::Rcout << "\n numerator = " << numerator;
	Rcpp::Rcout << "\n denominator = " << denominator;
	Rcpp::Rcout << "\n numerator == denominator? " << 
	  (numerator == denominator);
	Rcpp::Rcout << "\n is r > 1? " << (r > 1.0);
	Rcpp::Rcout << "\n is r < 0? " << (r < 0.0);
	Rcpp::Rcout << "\n r = " << r;
	Rcpp::Rcout << "\n is r == 1? " << (r == 1.0L);
	Rcpp::Rcout << "\n oneminusr = " << oneminusr;
	Rcpp::Rcout << "\n is oneminusr == 0? " << (oneminusr == 0.0L);
	Rcpp::Rcout << "\n r1 = " << r1;
	Rcpp::Rcout << "\n is r1 == 1? " << (r1 == 1.0);
	Rcpp::Rcout << "\n is eq12 < r? " << (eq12 < r);
#endif
	++ti_dbl_min;
	ti = DBL_MIN;
	// Beware of this!!  throw std::range_error("ti set to DBL_MIN");
	// Do not exit. Record it. We check for it now in R code. Maybe
	// abort simulation and go to a new one?  FIXME
	// Rcpp::Rcout << "ti set to DBL_MIN\n";
	// Yes, abort because o.w. we can repeat it many, manu times
	// throw std::range_error("ti set to DBL_MIN");
	throw rerunExcept("ti set to DBL_MIN");
      }
      if(ti < 0.001) ++ti_e3;
      ti += currentTime;
    } 
  }
  return ti;
}

static double Algo2_st(const spParamsP& spP,
		       const double& ti) {

  // beware the use of t: now as it used to be, as we pass the value
  // and take the diff in here: t is the difference

  using namespace Rcpp ;
  double t = ti - spP.timeLastUpdate;
  
  // if (t == 0 ) {
  //   Rcpp::Rcout << "\n Entered Algo2 with t = 0\n" <<
  //     "    Is this a forced sampling case?\n";
  //     return num;
  // }

  if (spP.popSize == 0.0) {
    #ifdef DEBUGW
    Rcpp::Rcout << "\n Entered Algo2 with pop size = 0\n";
    #endif
    return 0.0;
  }
  
  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);
  double m; // the holder for the binomial

  double rnb; // holder for neg. bino. So we can check.
  double retval; //So we can check 

    if( (1.0 - pe/pm) > 1.0) {
      Rcpp::Rcout << "\n ERROR: Algo 2: (1.0 - pe/pm) > 1.0\n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2:  1 - pe/pm > 1");
    }

    if( (1.0 - pe/pm) < 0.0 ) {
      // Rcpp::Rcout << "\n ERROR: Algo 2, (1.0 - pe/pm) < 0.0 \n";
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: 1 - pe/pm < 0");
    }

    if( pb > 1.0 ) {
      // Rcpp::Rcout << "\n WARNING: Algo 2, pb > 1.0 \n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb > 1 ");
    }

    if( pb < 0.0 ) {
      // Rcpp::Rcout << "\n WARNING: Algo 2, pb < 0.0 \n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb < 0");
    }
    //}


  if( pe == pm ) {
    // Should never happen. Exact identity??
    Rcpp::Rcout << "\n WARNING: Algo 2: pe == pm \n" ;
	      // << "; pm = " << pm  << "; pe = " 
	      // << pe << " pe == 0? " << (pe == 0) << "\n";
      // t << "; R = " << R 
      // << "; W = " << W << "; death = " << death 
      // << "; growth = " << growth << "; pm = " << pm 
      // << "; pe = " << pe << std::endl;
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
    #ifdef DEBUGW // just checking
      if(m != 0.0) 
	Rcpp::Rcout << "\n WARNING: Algo 2: 0.0 < m < 0.5" <<std::endl;
    #endif    
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
    DP2(rnb); DP2(m); DP2(pe); DP2(pm);
    print_spP(spP);
    throw std::range_error("Algo 2: retval not finite");
  }
  if( std::isnan(retval) )  {
    DP2(rnb); DP2(m); DP2(pe); DP2(pm);
    print_spP(spP);
    throw std::range_error("Algo 2: retval is NaN");
  }
  return retval;
}

static double Algo3_st(const spParamsP& spP, const double& t){
  
  using namespace Rcpp ;

  // double pm = pM_f(t, spP.R, spP.W, spP.death, spP.birth);
  // double pe = pE_f(pm, spP.W, spP.death, spP.birth);
  // double pb = pB_f(pe, spP.death, spP.birth);

  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);

  double m; // the holder for the binomial
  double retval;
  double rnb;

  if( (1.0 - pe/pm) > 1.0) {
    // Rcpp::Rcout << "\n ERROR: Algo 3: (1.0 - pe/pm) > 1.0\n"; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3:  1 - pe/pm > 1");
  }
  
  if( (1.0 - pe/pm) < 0.0 ) {
    // Rcpp::Rcout << "\n ERROR: Algo 3, (1.0 - pe/pm) < 0.0\n ";
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: 1 - pe/pm < 0");
  }
  
  if( pb > 1.0 ) {
    // Rcpp::Rcout << "\n WARNING: Algo 3, pb > 1.0\n "; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb > 1 ");
  }
  
  if( pb < 0.0 ) {
    // Rcpp::Rcout << "\n WARNING: Algo 3, pb < 0.0\n "; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb < 0");
  }
  
  if( pe == pm ) {
    // Should never happen. Exact identity??
    Rcpp::Rcout << "\n WARNING: Algo 3: pm == pe\n"; 
    // << "; t = " << 
    // 	 t << "; R = " << R 
    // 	      << "; W = " << W << "; death = " << death 
    // 	      << "; growth = " << growth << "; pm = " << pm 
    // 	      << "; pb = " << pb << std::endl;
    
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
    DP2(rnb); DP2(m); DP2(pe); DP2(pm);
    print_spP(spP);
    // Rcpp::Rcout << "\n ERROR: Algo 3, retval not finite\n ";
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb 
	      // << "; m = " << m << " ; rnb = " << rnb << std::endl;
    throw std::range_error("Algo 3: retval not finite");
  }
  if( !std::isfinite(retval) )  {
    DP2(rnb); DP2(m); DP2(pe); DP2(pm);
    print_spP(spP);
    throw std::range_error("Algo 3: retval is NaN");
  }
  return retval;
}




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


// this is the log of the ratio of death rates
// so the the difference of the successie death rates, if using
// the log version.
static double returnMFE(double& e1,
			const double& K,
			const std::string& typeFitness) {
  if((typeFitness == "mcfarland0") || (typeFitness == "mcfarlandlog"))
    return log(e1);
  else if(typeFitness == "mcfarland")
    return ((1.0/K) * e1);
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
				  const std::string& typeFitness,
				  const double& totPopSize,
				  const double& K){
  //				  const double& initSize) {
  // static double tps_0 = initSize;
  // static double tps_1 = 0.0;

  if( (typeFitness == "mcfarland0") ||
      (typeFitness == "mcfarland") || 
      (typeFitness == "mcfarlandlog") ) {
    double etmp;
    tps_1 = totPopSize;
    if(typeFitness == "mcfarland")
      etmp = std::abs( tps_1 - (tps_0 + 1) );
    else {
      if( (tps_0 + 1.0) > tps_1 ) 
	etmp = (K + tps_0 + 1.0)/(K + tps_1);
      else
	etmp = (K + tps_1)/(K + tps_0 + 1);
    }
    if(etmp > e1) {
      e1 = etmp;
      n_0 = tps_0;
      n_1 = tps_1;
    }
    tps_0 = tps_1;
  }
}


static void updateRatesMcFarland(std::vector<spParamsP>& popParams,
				 double& adjust_fitness_MF,
				 const double& K,
				 const double& totPopSize){

  adjust_fitness_MF = totPopSize/K;

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}


static void updateRatesMcFarlandLog(std::vector<spParamsP>& popParams,
				    double& adjust_fitness_MF,
				    const double& K,
				    const double& totPopSize){

  adjust_fitness_MF = log1p(totPopSize/K);

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}


// McFarland0 uses: - penalty as log(1 + N/K), and puts
// that in the birth rate.
static void updateRatesMcFarland0(std::vector<spParamsP>& popParams,
				  double& adjust_fitness_MF,
				  const double& K,
				  const double& totPopSize,
				  const int& mutatorGenotype,
				  const double& mu){
  
  adjust_fitness_MF = 1.0 / log1p(totPopSize/K);

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].birth = adjust_fitness_MF * popParams[i].absfitness;
    if(mutatorGenotype) {
      popParams[i].mutation = mu * popParams[i].birth * 
	popParams[i].numMutablePos;
    } else if(popParams[i].birth / popParams[i].mutation < 20) {
      Rcpp::Rcout << "\n WARNING: birth/mutation < 20";
      Rcpp::Rcout << "\n Birth = " << popParams[i].birth 
		<< ";  mutation = " << popParams[i].mutation << "\n";
    }
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}





static void updateRatesBeeren(std::vector<spParamsP>& popParams,
			      double& adjust_fitness_B,
			      const double& initSize,
			      const double& currentTime,
			      const double& alpha,
			      const double& totPopSize,
			      const int& mutatorGenotype,
			      const double& mu){

  double average_fitness = 0.0; // average_fitness in Zhu
  double weighted_sum_fitness = 0.0;
  double N_tilde;
  
  for(size_t i = 0; i < popParams.size(); ++i) {
    weighted_sum_fitness += (popParams[i].absfitness * popParams[i].popSize);
  }
  
  average_fitness = (1.0/totPopSize) * weighted_sum_fitness;
  N_tilde =  initSize * exp(alpha * average_fitness * currentTime);
  adjust_fitness_B = N_tilde/weighted_sum_fitness; 

  if(adjust_fitness_B < 0) {
    throw std::range_error("adjust_fitness_B < 0");
  }
  
  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].birth = adjust_fitness_B * popParams[i].absfitness;
    if(mutatorGenotype) {
      popParams[i].mutation = mu * popParams[i].birth * 
	popParams[i].numMutablePos;
    } else if(popParams[i].birth / popParams[i].mutation < 20) {
      Rcpp::Rcout << "\n WARNING: birth/mutation < 20";
      Rcpp::Rcout << "\n Birth = " << popParams[i].birth 
		<< ";  mutation = " << popParams[i].mutation << "\n";
    }
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}


// FIXME: split into two functions and 
// specialize for types of fitness if this takes
// a good chunck of time.
// Specialized now for fitness as in Bozic, continuous,
// model 1, to deal with non-satisfied drivers.


static void fitness(spParamsP& tmpP,
		    const spParamsP& parentP,
		    const int& mutatedPos, 
		    Rcpp::IntegerMatrix restrictTable,
		    const std::string& typeCBN,
		    const Genotype64& newGenotype,
		    const double& birthRate, 
		    const double& s,
		    // const double& death,
		    const int& numDrivers,
		    const std::string& typeFitness,
		    const double& genTime,
		    const double& adjust_fitness_B,
		    const double& sh,
		    const double& adjust_fitness_MF) {

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
    } else if (typeFitness == "bozic2") {
      double pp = pow( 1.0 - s, sumDriversMet) * 
	pow( 1.0 + sh, sumDriversNoMet);
      tmpP.birth = (1.0/genTime) * (1.0 - 0.5 * pp );
      tmpP.death = (0.5/genTime) * pp;
    } else if(typeFitness == "beerenwinkel") {
      // like Datta et al., 2013
      tmpP.absfitness = pow(1.0 + s, sumDriversMet) * 
	pow( 1.0 - sh, sumDriversNoMet);
      tmpP.birth = adjust_fitness_B * tmpP.absfitness;
    } else if(typeFitness == "mcfarland0") {
      tmpP.absfitness = pow(1.0 + s, sumDriversMet) / 
	pow( 1.0 + sh, sumDriversNoMet);
      tmpP.birth = adjust_fitness_MF * tmpP.absfitness;
    } else if(typeFitness == "mcfarland") {
      tmpP.birth = pow(1.0 + s, sumDriversMet) / 
	pow( 1.0 + sh, sumDriversNoMet);
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

    } else if (typeFitness == "log") {
      tmpP.birth = birthRate + s * log1p(sumDriversMet) - 
	sh * log(1 + sumDriversNoMet);
    } else { // linear
      tmpP.birth = birthRate + s * static_cast<double>(sumDriversMet) - 
	sh * static_cast<double>(sumDriversNoMet);
    } 
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

static inline void mapTimes_updateP(std::multimap<double, int>& mapTimes,
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


static inline void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes) {
  // we want minNextMutationTime and nextMutant
  nextMutant = mapTimes.begin()->second;
  minNextMutationTime = mapTimes.begin()->first;
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
					    const double& endTimeEvery,
					    const int& detectionDrivers,
					    const int& verbosity,
					    const double& fatalPopSize = 1e15) {
  // Fill out, but also compute totPopSize
  // and return sample summaries for popsize, drivers.
  
  // static int lastMaxDr = 0; // preserves value across calls to Algo5 from R.
  // so can not use it.
  bool storeThis = false;
  totPopSize = 0.0;
  
   // DP2(lastMaxDr);
  // DP2(detectionDrivers);
  // DP2(currentTime);
  // DP2((lastStoredSample + endTimeEvery));
  // DP2(detectionSize);

  int tmp_ndr = 0;
  int max_ndr = 0;

  for(size_t i = 0; i < popParams.size(); ++i) {
    totPopSize += popParams[i].popSize;
    tmp_ndr = count_NDrivers(Genotypes[i], NumDrivers);
    if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
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

    // if ( ( currentTime >= (lastStoredSample + endTimeEvery) ) &&
    // 	 ( (totPopSize >= detectionSize) ||
    // 	   (lastMaxDr >= detectionDrivers) ) ) {
    //   simulsDone = true;
    // }
  
  // Beware: this can lead to never stopping if
  // decreases in popSize or drivers

  // Logic: if a period k you meet any condition, recheck again at k +
  // endTimeEvery, and if conditions met exit. Prevents exiting if you
  // reach the cancer state almost by chance. But this is way too
  // paranoid. The idea is of application mainly for McF and Beeren
  // models, so we do not bail out as soon as just a single cell with one
  // new driver. But this makes things very slow.

  // Thus, never pass an endTimeEvery > 0, but use detectionDrivers = 1 +
  // intended final Drivers.

  
  if(endTimeEvery > 0) {
    if(done_at <= 0 ) {
      if( (totPopSize >= detectionSize) ||
	   (lastMaxDr >= detectionDrivers)  )
	done_at = currentTime + endTimeEvery;
    } else if (currentTime >= done_at) {
      if( (totPopSize >= detectionSize) ||
	  (lastMaxDr >= detectionDrivers)  ) {
	simulsDone = true;
	reachDetection = true;
      }
      else
	done_at = -9;
    }
  } else if( (totPopSize >= detectionSize) ||
	     (lastMaxDr >= detectionDrivers) )  {	
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


static inline void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
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
      strDrivers += (comma + std::to_string(i + 1)); //SSTR(i + 1));
      // strDrivers += (comma + SSTR(i + 1));
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
			     const double& tSample){

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
	Algo2_st(popParams[i], tSample);
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



static void precissionLoss(){
  // We are storing population sizes as doubles.
  // Should not loose any precission up to 2^53 - 1
  // (e.g., http://stackoverflow.com/a/1848762)
  // but double check if optims break it.
  // Note that the original code by Mather stores it as int.

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



static void innerBNB(const int& numGenes,
		     const double& initSize,
		     const double& K,
		     const double& alpha,
		     const std::string& typeCBN,
		     const double& genTime,
		     const std::string& typeFitness,
		     const int& mutatorGenotype,
		     const double& mu,
		     const double& sh,
		     const double& s,
		     const double& death,
		     const double& birthRate,
		     const double& keepEvery,
		     const double& sampleEvery,		     
		     const int& numDrivers,
		     const int& initMutant,
		     const time_t& start_time,
		     const double& maxWallTime,
		     const double& finalTime,
		     const double& detectionSize,
		     const double& endTimeEvery,
		     const int& detectionDrivers,
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


  
  // Beerenwinkel
  double adjust_fitness_B = -std::numeric_limits<double>::infinity();
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
  if(initMutant >= 0) {
    popParams[0].numMutablePos = numGenes - 1;
    Genotypes[0].set(initMutant);
    if(typeFitness == "beerenwinkel") {
      popParams[0].death = 1.0; //note same is in McFarland.
      // But makes sense here; adjustment in beerenwinkel is via fitness
      
      // initialize to prevent birth/mutation warning with Beerenwinkel
      // when no mutator. O.w., the defaults
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0 + s;
      updateRatesBeeren(popParams, adjust_fitness_B, initSize,
			currentTime, alpha, initSize, 
			mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland0") {
      // death equal to birth of a non-mutant.
      popParams[0].death = log1p(totPopSize/K); // log(2.0), except rare cases
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0 + s;
      updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
			    totPopSize,
			    mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland") {
      popParams[0].death = totPopSize/K;
      popParams[0].birth = 1.0 + s;
    } else if(typeFitness == "mcfarlandlog") {
      popParams[0].death = log1p(totPopSize/K);
      popParams[0].birth = 1.0 + s;
    } else if(typeFitness == "bozic1") {
      tmpParam.birth =  1.0;
      tmpParam.death = -99.9;
    } else if (typeFitness == "bozic2") {
      tmpParam.birth =  -99;
      tmpParam.death = -99;
    } else if (typeFitness == "exp") {
      tmpParam.birth =  -99;
      tmpParam.death = death;
    } else { // linear or log
      tmpParam.birth =  -99;
      tmpParam.death = death;
    } 
    if( (typeFitness != "beerenwinkel") && (typeFitness != "mcfarland0") 
	&& (typeFitness != "mcfarland") && (typeFitness != "mcfarlandlog")) // wouldn't matter
      fitness(popParams[0], tmpParam, initMutant, restrictTable,
	      typeCBN, Genotypes[0], birthRate, s, numDrivers, 
	      typeFitness, genTime, adjust_fitness_B, sh,
	      adjust_fitness_MF);
    // we pass as the parent the tmpParam; it better initialize
    // everything right, or that will blow. Reset to init
    init_tmpP(tmpParam);
  } else {
    popParams[0].numMutablePos = numGenes;
    if(typeFitness == "beerenwinkel") {
      popParams[0].death = 1.0;
      // initialize to prevent birth/mutation warning with Beerenwinkel
      // when no mutator. O.w., the defaults
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0;
      updateRatesBeeren(popParams, adjust_fitness_B, initSize,
			currentTime, alpha, initSize, 
			mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland0") {
      popParams[0].death = log1p(totPopSize/K);
      if(!mutatorGenotype)
	popParams[0].mutation = mu * popParams[0].numMutablePos;
      popParams[0].absfitness = 1.0;
      updateRatesMcFarland0(popParams, adjust_fitness_MF, K, 
			    totPopSize,
			    mutatorGenotype, mu);
    } else if(typeFitness == "mcfarland") {
      popParams[0].birth = 1.0;
      popParams[0].death = totPopSize/K;
      // no need to call updateRates
    } else if(typeFitness == "mcfarlandlog") {
      popParams[0].birth = 1.0;
      popParams[0].death = log1p(totPopSize/K);
      // no need to call updateRates
    } else if(typeFitness == "bozic1") {
      popParams[0].birth = 1.0;
      popParams[0].death = 1.0;
    } else if (typeFitness == "bozic2") {
      popParams[0].birth = 0.5/genTime;
      popParams[0].death = 0.5/genTime;
    } else if (typeFitness == "exp") {
      popParams[0].birth = 1.0;
      popParams[0].death = death;
    } else { // linear or log
      popParams[0].birth = birthRate;
      popParams[0].death = death;
    }
  }


  // // FIXME debug
  //     Rcpp::Rcout << " popSize[0]  at 2 ";
  //     print_spP(popParams[0]);
  //     // end debug
  

  
  // these lines (up to, and including, R_F_st)
  // not needed with mcfarland0 or beerenwinkel
  if(mutatorGenotype)
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
      // now total number of different species
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
		typeCBN, newGenotype, birthRate, s,
		numDrivers, typeFitness, genTime,
		adjust_fitness_B, sh, adjust_fitness_MF);
	

	if(tmpParam.birth > 0.0) {
	  tmpParam.numMutablePos = numMutablePosParent - 1;
	  if(mutatorGenotype)
	    tmpParam.mutation = mu * tmpParam.birth * tmpParam.numMutablePos;
	  //	    tmpParam.mutation = mu * tmpParam.birth * (numMutablePosParent - 1);
	  else
	    tmpParam.mutation = mu * tmpParam.numMutablePos;
	    //tmpParam.mutation = mu * (numMutablePosParent - 1);
	  if (tmpParam.mutation > 1 )
	    Rcpp::Rcout << "WARNING: mutation > 1\n";
	  if (numMutablePosParent == 1) 
	    Rcpp::Rcout << "Note: mutation = 0; no positions left for mutation\n";
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
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
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
	    Algo2_st(popParams[sp], currentTime);
	  if(verbosity >= 2) {
	    Rcpp::Rcout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  }
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	}
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	//popParams[sp].Flag = true;
      }
      //   ***************  5.7 ***************
      // u_2 irrelevant if to_update = 1;
      u_1 = nextMutant;
      u_2 = static_cast<int>(sp);
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
				      reachDetection,
				      lastMaxDr,
				      done_at,
				      Genotypes, popParams, 
				      currentTime,
				      numDrivers,
				      keepEvery,
				      detectionSize,
				      finalTime,
				      endTimeEvery,
				      detectionDrivers,
				      verbosity); //keepEvery is for thinning
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

      if( (typeFitness == "beerenwinkel") ) {
	updateRatesBeeren(popParams, adjust_fitness_B,
			  initSize, currentTime, alpha, totPopSize,
			  mutatorGenotype, mu);
      } else if( (typeFitness == "mcfarland0") ) {
	updateRatesMcFarland0(popParams, adjust_fitness_MF,
			     K, totPopSize,
			     mutatorGenotype, mu);
      } else if( (typeFitness == "mcfarland") ) {
	updateRatesMcFarland(popParams, adjust_fitness_MF,
			     K, totPopSize);
      } else if( (typeFitness == "mcfarlandlog") ) {
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

SEXP BNB_Algo5(SEXP restrictTable_,
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
	       SEXP detectionDrivers_,
	       SEXP onlyCancer_,
	       SEXP errorHitWallTime_,
	       SEXP maxNumTries_,
	       SEXP errorHitMaxTries_
	       ) {

  BEGIN_RCPP
    using namespace Rcpp;
  precissionLoss();
  const IntegerMatrix restrictTable(restrictTable_);
  const int numDrivers = as<int>(numDrivers_);
  const int numGenes = as<int>(numGenes_);
  const std::string typeCBN = as<std::string>(typeCBN_);
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
  const int verbosity = as<int>(verbose_);
  // const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // species
  const int seed = as<int>(seed_gsl_);
  const long maxram = as<int>(maxram_);
  const int mutatorGenotype = as<int>(mutatorGenotype_);
  const int initMutant = as<int>(initMutant_);
  const double maxWallTime = as<double>(maxWallTime_);
  const double keepEvery = as<double>(keepEvery_);

  
  const double alpha = as<double>(alpha_);
  const double sh = as<double>(sh_); // coeff for fitness
  // if a driver without dependencies. Like in Datta et al., 2013.
  const double K = as<double>(K_); //for McFarland
  const double endTimeEvery = as<double>(endTimeEvery_); 
  const int detectionDrivers = as<int>(detectionDrivers_); 
  const double genTime = 4.0; // should be a parameter. For Bozic only.
  // const bool errorFinalTime = as<bool>(errorFinalTime_);
  const bool errorHitWallTime = as<bool>(errorHitWallTime_);
  const bool onlyCancer = as<bool>(onlyCancer_);
  const int maxNumTries = as<int>(maxNumTries_);
  const bool errorHitMaxTries = as<bool>(errorHitMaxTries_);
  
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
    throw std::range_error("This version only accepts up to 64 genes.");

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
	       alpha,
	       typeCBN,
	       genTime,
	       typeFitness,
	       mutatorGenotype,
	       mu,
	       sh,
	       s,
	       death,
	       birthRate,
	       keepEvery,
	       sampleEvery,		     
	       numDrivers,
	       initMutant,
	       start_time,
	       maxWallTime,
	       finalTime,
	       detectionSize,
	       endTimeEvery,
	       detectionDrivers,
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
	List::create(Named("other") =
		     List::create(Named("UnrecoverExcept") = true,
				  Named("ExceptionMessage") = e.what()));
    } catch (...) {
      Rcpp::Rcout << "\n Unknown unrecoverable exception. Aborting. \n";
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
		       Named("other") =
		       List::create(Named("UnrecoverExcept") = false));
      }
    } else if(numRuns > maxNumTries) {
      //  hittedMaxTries
      hittedMaxTries = true;
      Rcpp::Rcout << "\n Hitted maxtries. Exiting.";
      runAgain = false;
      if(errorHitMaxTries) {
	Rcpp::Rcout << "\n Hitting max tries is regarded as an error. \n";
	return
	  List::create(Named("HittedMaxTries") = true,
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
  IntegerMatrix returnGenotypes(numGenes, uniqueGenotypes_vector.size());
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
  NumericMatrix outNS(outNS_r, outNS_c);  
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

  NumericMatrix perSampleStats(outNS_i + 1, 5);
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
		 Named("other") = List::create(Named("attemptsUsed") = numRuns,
					       Named("errorMF") = 
					       returnMFE(e1, K, 
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

  END_RCPP
    
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



// static inline void fitness_linear_std_bitset64(double& birth,
// 					       double& death,
// 					       const Genotype64& Genotype,
// 					       const double& birthRate, 
// 					       const double& s, const int& numDrivers) {
//   // Crucial: drivers are always first in genotype
//   // For this approach, death should always be 1, and birth is irrelevant
//   int totalMut = 0;
//   for(int i = 0; i < numDrivers; ++i) {
//     totalMut += Genotype[i];
//   }
//   birth = birthRate + s * static_cast<double>(totalMut);
//   // d is not touched
// }


// static inline double fitness_beeren_std_bitset64(double& birth,
// 						 double& death,
// 						 const Genotype64& Genotype,
// 						 //const double& birthRate, 
// 						 const double& s, const int& numDrivers) {
//   // Crucial: drivers are always first in genotype
//   // An exponential-like, inspired in Beerenwinkel et al., 2007
//   int totalMut = 0;
//   for(int i = 0; i < numDrivers; ++i) {
//     totalMut += Genotype[i];
//   }
//   birth = pow( 1 + s, totalMut );
// }


// static inline void fitness_log_std_bitset64(double& birth,
// 					    double& death, 
// 					    const Genotype64& Genotype,
// 					    const double& birthRate,
// 					    const double& s, const int& numDrivers) {
//   // birth will often be 0.25
//   int totalMut = 0;
//   for(int i = 0; i < numDrivers; ++i) {
//     totalMut += Genotype[i];
//   }
//   birth =  birthRate + s * log( 1 + totalMut );
// }


// // plot(function(mut, s = 0.1) { (1+s) ^ mut}, 0, 10)
// // Another log-like, if I wanted to try:
// //  0.25 + s * log(mut + 1)
// //  plot(function(mut, s = 0.1) {0.25 + s * (log (mut + 1))}, 0, 15, main = "s = 0.1")



// static inline void fitness_bozic1_bitset64(double& death,
// 					   const Genotype64& Genotype,
// 					   const double& s, const int& numDrivers) {
//   // The first continous-time model
//   // p.5 of their supplementary material

//   int totalMut = 0;
//   for(int i = 0; i < numDrivers; ++i) {
//     totalMut += Genotype[i];
//   }
//   death = pow( 1.0 - s, totalMut); 
// }




// static double fitness_CBN_bitset64(const int& mutatedPos, 
// 				   Rcpp::IntegerMatrix restrictTable,
// 				   const double& fitnessParent, 
// 				   const std::string& typeCBN,
// 				   const Genotype64& Genotype,
// 				   const double& birthRate, 
// 				   const double& s, 
// 				   const int& numDrivers,
// 				   const std::string typeFitness) {
    
//   using namespace Rcpp ;

//   // I could simplify, by having the restrictTable only
//   // for drivers that DO depend. Would make it faster and exclude
//   // a check below. But less clear.

//   // will later become an argument
//   // double fitnessNo = 0.0;

//   int numDependencies;
//   int sumPresent = 0;
//   int outFitnessYes = 0;

//   if(mutatedPos >= numDrivers) { //the new mutation is a passenger
//     return fitnessParent;
//   } else {
//     const Rcpp::IntegerMatrix::Column thisRestrict = restrictTable(_, mutatedPos);
//     numDependencies = thisRestrict[1];
// #ifdef DEBUGW
//     if(thisRestrict[0] != mutatedPos ) {
//       Rcpp::Rcout << std::endl << "thisRestrict[0] = " << thisRestrict[0] 
// 		<< "; mutatedPos  = " << mutatedPos  << std::endl;
//       throw std::range_error("FitnessCBN: thisRestrict[0] != mutatedPos ");	
//     }
//     if(Genotype[mutatedPos] != 1){
//       Rcpp::Rcout << " mutatedPos = " << mutatedPos << std::endl;
//       throw std::range_error("FitnessCBN: genotype(mutatedPos) != 1");
//     }
// #endif
    
//     if(!numDependencies) { // mutated is driver with no deps
//       outFitnessYes = 1;
//     } else {
//       for(int i = 2; i < (2 + numDependencies); i++) {
// 	sumPresent += Genotype[ thisRestrict[i] ];
//       }
//       if(typeCBN == "Multiple") {
//         if(sumPresent) outFitnessYes = 1;
//       } else{ // if(typeCBN == "CBN")
//         if(sumPresent == numDependencies)
//           outFitnessYes = 1; 
//       }
//     }

//     // FIXME: this could be much improved!
//     // There is a part, above, for dependencies, and a part, below,
//     // for fitness.


    

//     if(outFitnessYes) {
//       if(typeFitness == "bozic1")
// 	fitness_bozic_bitset64(birth, death, Genotype, s, numDrivers);
//       else if (typeFitness == "beeren")
// 	fitness_beeren_std_bitset64(birth, death, Genotype, s, numDrivers);
//       else if (typeFitness == "log")
// 	fitness_log_std_bitset64(Genotype, birthRate, s, numDrivers);
//       else  
// 	fitness_linear_std_bitset64(birth, death, Genotype, birthRate, s, numDrivers);
//     } else {
//       death = 1.0;
//       birth = 0.0; // this is what really matters so that
//       // the pop does not get added.
//       // Line with comment "fitness is 0"
//     }
//   }
// }

// static inline double fitness_bozic_bitset64(const Genotype64& Genotype,
// 				   const double& s, const int& numDrivers) {
//   // This is the wrong one, the discrete time one
  

//   int totalMut = 0;
//   for(int i = 0; i < numDrivers; ++i) {
//     totalMut += Genotype[i];
//   }
//   return 0.5 * pow( 1.0 - s, totalMut); // This is death, not birth
//   // return 1.0 - 0.5 * pow( 1.0 - s, totalMut); 

// }
// plot(function(mut, s = 0.05) {1 - 0.5 * (1 - s)^mut} , 0, 15) ## birth
// plot(function(mut, s = 0.05) { 0.5 * ( (1 - s)^mut)}, 0, 15) ## death

// next two could be void
// inline double W_f_st(const spParamsP& spP){
//   return spP.death + spP.birth + spP.mutation;
// }

// inline double R_f_st(const spParamsP& spP) {
//   return sqrt( pow( spP.birth - spP.death, 2) + 
// 	       ( 2 * (spP.birth + spP.death) + 
// 		 spP.mutation) * spP.mutation );
// }

// does not work, as "<" not defined on bitset.
// see possible solution here http://en.allexperts.com/q/C-1040/inserting-bitset-set.htm
// but I'd rather make conversion explicit and under control
// inline void find_unique_genotypes_v0(std::set<Genotype64>& uniqueGenotypes,
// 				  const std::vector<Genotype64>& genot_out) {
//   for(unsigned i = 0; i < genot_out.size(); ++i) 
//     uniqueGenotypes.insert( genot_out[i] );
// }


// FIXME00: do not use const_iterator?
// void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
// 			     const std::vector<double>& nextMutTime) {
//   // we want minNextMutationTime and nextMutant
//   // does not work if species with popSize == 0

//   std::vector<double>::const_iterator pbeg = nextMutTime.begin();
//   std::vector<double>::const_iterator pend = nextMutTime.end();

//   std::vector<double>::const_iterator pt_pos_min =
//     std::min_element(pbeg, pend);

//   nextMutant = std::distance(pbeg, pt_pos_min);
//   minNextMutationTime = nextMutTime[nextMutant];
// }

// void resize_outNS(arma::mat& outNS, const int& outNS_i,
// 		  const int& numSpecies) {
//   int type_resize = 0;

//   if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
//   if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

//   if(type_resize == 0) return;
//   else if(type_resize == 1) 
//     outNS.resize(2 * numSpecies + 2, outNS.n_cols);
//   else if (type_resize == 2)
//     outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
//   else if (type_resize == 3)
//     outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);
// }

// static inline void create_returnGenotypes(Rcpp::IntegerMatrix& returnGenotypes,
// 					  const int& numGenes,
// 					  const std::vector<unsigned long>& uniqueGenotypesV){
// #ifdef DEBUGV
//   Rcpp::Rcout << "\n entering create_returnGenotypes\n";
// #endif
//   DP(uniqueGenotypesV.size());
//   DP(numGenes);

//   for(size_t i = 0; i < uniqueGenotypesV.size(); ++i) {
//     Genotype64 tmpbs(uniqueGenotypesV[i]);
//     for(int j = 0; j < numGenes; ++j) {
//       returnGenotypes(i, j) = tmpbs[j];
//     }
//   }

// #ifdef DEBUGV
//   Rcpp::Rcout << "\n exiting create_returnGenotypes\n";
// #endif

// }


// // cannot pass the structure for parent, for the case
// // of initial individual
// static void fitness(spParamsP& tmpP,
// 		    const double& parentBirth,
// 		    const double& parentDeath,
// 		    const int& mutatedPos, 
// 		    Rcpp::IntegerMatrix restrictTable,
// 		    const std::string& typeCBN,
// 		    const Genotype64& newGenotype,
// 		    const double& birthRate, 
// 		    const double& s,
// 		    const double& death,
// 		    const int& numDrivers,
// 		    const std::string typeFitness) {
  
//   // FIXME00: This is ugly!!

//   if(typeFitness == "bozic1") {
//     // if bozic, always death, so we pass death of parent
//     tmpP.death = fitness_CBN_bitset64(mutatedPos,
// 				      restrictTable,
// 				      parentDeath,
// 				      typeCBN,
// 				      newGenotype,
// 				      birthRate,
// 				      s,
// 				      numDrivers,
// 				      typeFitness);
//     // FIXME: there are cases when it is zero
//     tmpP.birth = 1 - tmpP.death;
//   } else if(typeFitness == "bozic2") {

//   } else {
//     tmpP.birth = fitness_CBN_bitset64(mutatedPos,
// 				      restrictTable,
// 				      parentBirth,
// 				      typeCBN,
// 				      newGenotype,
// 				      birthRate,
// 				      s,
// 				      numDrivers,
// 				      typeFitness);
//     tmpP.death = death;
//   }
// }

// fitness used up to version 1.0.15.
// static void fitness(spParamsP& tmpP,
// 		    const spParamsP& parentP,
// 		    const int& mutatedPos, 
// 		    Rcpp::IntegerMatrix restrictTable,
// 		    const std::string& typeCBN,
// 		    const Genotype64& newGenotype,
// 		    const double& birthRate, 
// 		    const double& s,
// 		    const double& death,
// 		    const int& numDrivers,
// 		    const std::string typeFitness,
// 		    const double& genTime,
// 		    const double& adjust_fitness_B,
// 		    const double& fitness_unmet) {

//   using namespace Rcpp;
//   // Two pieces: split into two functions??
//   //    - checking restrictions
//   //    - returning actual fitness according


//   int numDependencies;
//   int sumPresent = 0;
//   int outFitnessYes = 0;
//   //outFitnessYes = 1 means restrictions are met

//   // Set to default values.
//   //   - allows immediate exit if passenger
//   //   - sets correct defaults for all other fitness functs
//   tmpP.birth = parentP.birth;
//   tmpP.death = parentP.death;
//   tmpP.absfitness = parentP.absfitness;

//   //      **** Are driver constraints met ***
//   // I could simplify, by having the restrictTable only
//   // for drivers that DO depend. Would make it faster and exclude
//   // a check below. But less clear.

//   if(mutatedPos >= numDrivers) { //the new mutation is a passenger
//     return;
//   } else {
//     const Rcpp::IntegerMatrix::Column thisRestrict = restrictTable(_, mutatedPos);
//     numDependencies = thisRestrict[1];
// #ifdef DEBUGW
//     if(thisRestrict[0] != mutatedPos ) {
//       Rcpp::Rcout << std::endl << "thisRestrict[0] = " << thisRestrict[0] 
// 		<< "; mutatedPos  = " << mutatedPos  << std::endl;
//       throw std::range_error("FitnessCBN: thisRestrict[0] != mutatedPos ");	
//     }
//     if(newGenotype[mutatedPos] != 1){
//       Rcpp::Rcout << " mutatedPos = " << mutatedPos << std::endl;
//       throw std::range_error("FitnessCBN: genotype(mutatedPos) != 1");
//     }
// #endif
    
//     if(!numDependencies) { // mutated is driver with no deps
//       outFitnessYes = 1;
//     } else {
//       for(int i = 2; i < (2 + numDependencies); i++) {
// 	sumPresent += newGenotype[ thisRestrict[i] ];
//       }
//       if(typeCBN == "Multiple") {
//         if(sumPresent) outFitnessYes = 1;
//       } else{ // if(typeCBN == "CBN")
//         if(sumPresent == numDependencies)
//           outFitnessYes = 1; 
//       }
//     }

    
//     if(outFitnessYes) {
//       int totalMut = 0;
//       for(int i = 0; i < numDrivers; ++i) {
// 	totalMut += newGenotype[i];
//       }
//       //tmpP.totalMut = totalMut; 
//       if(typeFitness == "beerenwinkel") {
// 	tmpP.absfitness = pow(1.0 + s, totalMut);
// 	tmpP.birth = adjust_fitness_B * tmpP.absfitness;
//       } else if(typeFitness == "bozic1") {
// 	tmpP.death = pow( 1.0 - s, totalMut); 
//       } else if (typeFitness == "bozic2") {
// 	double pp = pow( 1.0 - s, totalMut);
// 	tmpP.birth = (1.0/genTime) * (1.0 - 0.5 * pp );
// 	tmpP.death = (0.5/genTime) * pp;
//       } else if (typeFitness == "exp") {
// 	tmpP.birth = pow( 1 + s, totalMut );
//       } else if (typeFitness == "log") {
// 	tmpP.birth =  birthRate + s * log( 1 + totalMut );
//       }
//       else { // linear
// 	tmpP.birth = birthRate + s * static_cast<double>(totalMut);
//       } 
//     } else { // restrictions not met
//       if(fitness_unmet < 0) {
// 	tmpP.death = 1.0;
// 	tmpP.birth = 0.0; // this is what really matters so that
// 	// the pop does not get added.
// 	// Line with comment "fitness is 0"
//       } else { // penalization
// 	if(fitness_unmet == 1.0) { // common special case
// 	  return;
// 	  // above we did:
// 	  // tmpP.death = parentDeath;
// 	  // tmpP.birth = parentBirth;
// 	} else {
// 	  if( (typeFitness == "beerenwinkel") ) {
// 	    tmpP.absfitness = parentP.absfitness * fitness_unmet;
// 	    tmpP.birth = parentP.birth * fitness_unmet;
// 	  } else if ( (typeFitness == "exp") ||
// 		      (typeFitness == "log") ||
// 		      (typeFitness == "linear") ) {
// 	    tmpP.birth = parentP.birth * fitness_unmet;
// 	  } else { // both bozics are the same
// 	    tmpP.death = parentP.death * (1.0/fitness_unmetn);
// 	  }
// 	}
//       }
//     }
//   } 
// }


// Why the above was wrong

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
