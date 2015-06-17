// Clean up code in general:
//  - recheck time used by R
//  - minimize debugging and verbose code, or wrap in DEBUG statements
//  - remove all unused code from former functions
//  - use "static"??
//  - profile again, using all current fitness functions
//  - compile after profiling with gcc, after deciding function to use



//  - Explore params
//  - Check alternatives to Bozic:
//    + density-dep


// P is better than O in every case (and both always better or much, much
// [> 10x]) better than M, K, etc. So use P.  Sampling 10 seems generally
// good, although occasionally a very slow result.  20 is generally
// slow. 5 is about as 10, but larger object sizes.  2 can be about as 5,
// or slightly worse, with much larger object sizes generally.  1 is often
// slow and very large objects.


// Use templates
// Put timers in specific parts of the code.
// Clean up multimap: 
//        - profile with gooble prof in AMDs
//        - profile with different sampling periods
//        - run valgrind





// FIXME: with bitset, fix printing so only the relevant genes




// FIXME: install and use AMD libraries and recompile R, etc?
// http://blogs.amd.com/developer/2012/04/23/gcc-4-7-is-available-with-support-for-amd-opteron%E2%84%A2-6200-series-and-amd-fx-series-processors/
// Open64: http://devgurus.amd.com/message/1282903#1282903
// With gcc and new flags, I see little or no difference.


// getMinNextMutationTime: ~ 20%

// new_sp_gmp: ~ 40% but long int version is about twice as fast (and
// new_sp only takes about 5%) of time.

// So for now, use version H which is about 3 times as fast as version I
// (even faster in the AMD), and about twice as fast as version G. Version
// H, though, only allows for genotypes up to 60 genes. Oh well; will provide 
// automatic jump to G if more than those (or to I, if later it is faster)

// FIXME: max_remove: that is ugly and way too large.


// Time speed of removal of zero species (not very worth it: profiling
// suggests less than 7% of time in Intel, 8% in AMD ) [Leaving debugging
// library when linking apparently no performance penalty]



// Remove all zero, even if non-zero previously


// FIXME: how to return memory to R even if we C-c?


// FIXME: keep track of ancestor: A string (we keep concatenating ancestors)
//        added to structure? Pass to R as a different object, not 
//        returnParams

// FIXME: use enums for typeFitness and typeCBN. At least check values are
// within allowed ones. At start of Algo5.


// FIXME: openMP?? The loop at 5.2 and at 5.9 
// but I call R code for random number generator.


// FIXME: remove species?
//       - For drivers, I definitely want to know their fate to plot it
//       - For the rest maybe not, if just sampling at end.
//       - removing from an std::vector is expensive (more than carrying them??)
//   So:
//       - When a species (any, with or without drivers) becomes popSize 0
//         move it to a point beyond the last non-zero species
//       - Use a counter (which is not numSpecies)
//       - In 5.2 and 5.9, cycle up to that counter.
//       - Actually, two counters: numExist and numZero (or lastIndexZero)
// 



#include "matherR.h"
#include <limits>
#include <iostream>
#include <gsl/gsl_rng.h> // here? in the .h
// #include <gmpxx.h>
// #include <queue>
#include <bitset>
#include <set>
#include <iterator>
#include <map>

// #define DEBUGV
// #define DEBUGW

#ifdef DEBUGW
#define ASSERT(x) {							\
    if (! (x)) {							\
      std::cout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	std::cout << " on line " << __LINE__  << "\n\n";		\
    }									\
  }
#else
#define ASSERT(x);
//#define ARMA_NO_DEBUG
#endif


#ifdef DEBUGW
#define STOPASSERT(x) {							\
    if (! (x)) {							\
      std::cout << "\n\nERROR!! Assertion  " << #x << " failed\n";	\
	std::cout << " on line " << __LINE__  << "\n\n";		\
	throw std::out_of_range("STOPASSERT");				\
    }									\
  }
#else
#define STOPASSERT(x);
#endif


void here(std::string x) {
std::cout << "\n HERE " << x << "\n";
}



// #define myT short int
typedef int myT;

// #include <typeinfo> // FIXME: remove later
// #include "output_any.h" // FIXME: remove later
// #include "output_any.h" // FIXME: remove later

typedef std::bitset<64> Genotype64;



struct spParamsK {
  bool Flag;
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  double nextMutationTime;
  double timeLastUpdate;
};

struct spParamsM {
  bool Flag;
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  //double nextMutationTime;
  double timeLastUpdate;
  std::multimap<double, int>::iterator pv;
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
};


inline double W_f(const double& death, const double& growth, 
		  const double& mu) {
  return death + growth + mu;
}

inline double R_f(const double& death, const double& growth, 
		  const double& mu) {
  return sqrt( pow( growth - death, 2) + ( 2 * (growth + death) + mu) * mu  );
}


inline double W_f_st(const spParamsP& spP){
  return spP.death + spP.birth + spP.mutation;
}

inline double R_f_st(const spParamsP& spP) {
  return sqrt( pow( spP.birth - spP.death, 2) + 
	       ( 2 * (spP.birth + spP.death) + 
		 spP.mutation) * spP.mutation );
}

double pM_f(const double& t, const double& R, const double& W, 
	    const double& death, const double& growth) {

  long double Ct = cosh(R * t/2.0);
  long double St = sinh(R * t/2.0);
  long double lpM = -99.99;
  
  if( (!std::isfinite(Ct) ) || (!std::isfinite(St)) ) {
    throw std::range_error("pM.f: Ct or St too big");
  }
    
  // my expression, which I think is better       
  //  lpM = (R * Ct + St * (2.0 * death - W ))/(R * Ct + St * (W - 2.0 * growth));
  // theirs, in paper and code
  lpM = (R * Ct + 2.0 * death * St - W * St)/(R * Ct - 2.0 * growth * St + W * St);    
  
  double pM = static_cast<double>(lpM);

  if( !std::isfinite(pM) ) {
    throw std::range_error("pM.f: pM not finite");
  }
  if(pM <= 0.0) std::cout << "\nWARNING: pM <= 0.0. Get ready for problems.\n";
  return pM;

}


double pM_f_st(const double& t, 
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
    throw std::range_error("pM.f: pM not finite");
  }
  if(pM <= 0.0) std::cout << "\nWARNING: pM <= 0.0. Get ready for problems.\n";
  return pM;

}


inline double pE_f(const double& pM, const double& W, 
		   const double& death, const double& growth) {
  double pE = (death * (1.0 - pM ) )/(W - death - growth * pM );
  if( !std::isfinite(pE) ) {
    throw std::range_error("pE.f: pE not finite");
  }
  return pE;
}


inline double pE_f_st(double& pM, const spParamsP& spP){
  double pE = (spP.death * (1.0 - pM ) )/(spP.W - spP.death - spP.birth * pM );
  if( !std::isfinite(pE) ) {
    throw std::range_error("pE.f: pE not finite");
  }
  return pE;
}


inline double pB_f(const double& pE, const double& death, 
		   const double& growth) {
  return (growth * pE)/death; 
}

inline double pB_f_st(const double& pE,
		      const spParamsP& spP) {
  return (spP.birth * pE)/spP.death; 
}

double ti_nextTime_tmax_2(const double& R, const double& W, 
			  const double& death, const double& growth,
			  const double& mu,
			  const double& n, const double& currentTime,
			  const double& tSample) {
  // Following the logic of the code by Mather in
  // findNextMutationTime

  // We return the nextMutationTime or a value larger than the
  // max length of the period (tSample)

  // I also change names rr, r, to match those in Mather r1, r.
  
  // However, I pass mutation, and split computation to avoid numerical problems
  // I was getting ti == 0 and ti < 0 in the other versions with large N.
  using namespace Rcpp ;

  double r1;
  double r;
  double ti;
  double pM;
  
  double numerator;
  double denominator;

  const double epsilon = 10.0;

  // W < 0 is a signal that mutation is zero, and thus ti is Inf
  if(W <= -90.0) {
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
    r = pow(r1, 1.0/n);

    pM = pM_f(tSample - currentTime, R, W, death, growth);

    if( r < pM) {// time to mutation longer that this time period
      ti = tSample + epsilon;
    } else {
      // Expand numerator and denominatior, the term for W and simplify.
      // Then, express as (1- r) (which is, inclussively, between 0 and 1)
      // and then multiply by -1 to take the log of each
      numerator =  (1.0 - r) * (R + growth - death - mu) + 2.0 * mu;
      denominator = (1.0 - r) * (growth - death - mu - R ) + 2.0 * mu;

      // FIXME? is it really necessary to use log(-a) - log(-b) or could
      // I just use log(a/b), where a and b are -numerator and -denominator?

      ti = (1.0/R) * (log(numerator) - log(denominator));
	
      //eq. 11
      // ti = (1.0/R) * (log( -1 * (r * (R - W + 2.0 * growth) - W - R + 2.0 * death )) -
      //  		      log( -1 * (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death )));

      // ti = (1.0/R) * log( (r * (R - W + 2.0 * growth) - W - R + 2.0 * death) /
      //               (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death));
      // std::cout << "\n this is ti = " << ti << "\n";
      if(ti < 0.0) {
	
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);

	std::cout << "\n ERROR: ti: eq.11 < 0 \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n mu = " << mu;
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";

	throw std::range_error("ti: eq.11 < 0");
      } 
      if( !std::isfinite(ti) ) {
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);

	std::cout << "\n ERROR: ti not finite \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	throw std::range_error("ti: ti not finite");
      }
      if(ti == 0.0) {
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);
	std::cout << "\n WARNING: ti == 0. Expect problems \n" 
		  << " especially if using forced sampling \n";

	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
      }
      ti += currentTime;
    } 
  }
  return ti;
}

double ti_nextTime_tmax_2_st(const spParamsP& spP,
			     const double& currentTime,
			     const double& tSample) {
  // Following the logic of the code by Mather in
  // findNextMutationTime

  // We return the nextMutationTime or a value larger than the
  // max length of the period (tSample)

  // I also change names rr, r, to match those in Mather r1, r.
  
  // However, I pass mutation, and split computation to avoid numerical problems
  // I was getting ti == 0 and ti < 0 in the other versions with large N.
  using namespace Rcpp ;

  double r1;
  double r;
  double ti;
  double pM;
  
  double numerator;
  double denominator;

  const double epsilon = 10.0;

  // W < 0 is a signal that mutation is zero, and thus ti is Inf
  if(spP.W <= -90.0) {
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
    r = pow(r1, 1.0/spP.popSize);

    pM = pM_f_st(tSample - currentTime, spP);

    if( r < pM) {// time to mutation longer that this time period
      ti = tSample + epsilon;
    } else {
      // Expand numerator and denominatior, the term for W and simplify.
      // Then, express as (1- r) (which is, inclussively, between 0 and 1)
      // and then multiply by -1 to take the log of each
      numerator =  (1.0 - r) * (spP.R + spP.birth - spP.death - spP.mutation) 
	+ 2.0 * spP.mutation;
      denominator = (1.0 - r) * (spP.birth - spP.death - spP.mutation - spP.R ) 
	+ 2.0 * spP.mutation;

      // FIXME? is it really necessary to use log(-a) - log(-b) or could
      // I just use log(a/b), where a and b are -numerator and -denominator?

      ti = (1.0/spP.R) * (log(numerator) - log(denominator));
	
      //eq. 11
      // ti = (1.0/R) * (log( -1 * (r * (R - W + 2.0 * growth) - W - R + 2.0 * death )) -
      //  		      log( -1 * (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death )));

      // ti = (1.0/R) * log( (r * (R - W + 2.0 * growth) - W - R + 2.0 * death) /
      //               (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death));
      // std::cout << "\n this is ti = " << ti << "\n";
      if(ti < 0.0) {
	
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);

	std::cout << "\n ERROR: ti: eq.11 < 0 \n";
	// std::cout << "\n R = " << R;
	// std::cout << "\n W = " << W;
	// std::cout << "\n r1 = " << r1;
	// std::cout << "\n r = " << r;
	// std::cout << "\n n = " << n;
	// std::cout << "\n mu = " << mu;
	// std::cout << "\n growth = " << growth;
	// std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	throw std::range_error("ti: eq.11 < 0");
      } 
      if( !std::isfinite(ti) ) {
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);
	std::cout << "\n ERROR: ti not finite \n";
	// std::cout << "\n R = " << R;
	// std::cout << "\n W = " << W;
	// std::cout << "\n r1 = " << r1;
	// std::cout << "\n r = " << r;
	// std::cout << "\n n = " << n;
	// std::cout << "\n growth = " << growth;
	// std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	throw std::range_error("ti: ti not finite");
      }
      if(ti == 0.0) {
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) / 
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);
	std::cout << "\n WARNING: ti == 0. Expect problems \n" 
		  << " especially if using forced sampling \n";

	// std::cout << "\n R = " << R;
	// std::cout << "\n W = " << W;
	// std::cout << "\n r1 = " << r1;
	// std::cout << "\n r = " << r;
	// std::cout << "\n n = " << n;
	// std::cout << "\n growth = " << growth;
	// std::cout << "\n death = " << death << "\n";
	std::cout << "\n numerator = " << numerator;
	std::cout << "\n denominator = " << denominator;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
      }
      ti += currentTime;
    } 
  }
  return ti;
}

double Algo2(const double& num, const double& t, 
	     const double& R, 
	     const double& W, 
	     const double& death, 
	     const double& growth) {

  using namespace Rcpp ;
  
  // if (t == 0 ) {
  //   std::cout << "\n Entered Algo2 with t = 0\n" <<
  //     "    Is this a forced sampling case?\n";
  //     return num;
  // }

  if (num == 0.0) {
    #ifdef DEBUGW
    std::cout << "\n Entered Algo2 with pop size = 0\n";
    #endif
    return 0.0;
  }

  double pm = pM_f(t, R, W, death, growth);
  double pe = pE_f(pm, W, death, growth);
  double pb = pB_f(pe, death, growth);
  double m; // the holder for the binomial

  double rnb = -99; // holder for neg. bino. So we can check.
  double retval; //So we can check 

    if( (1.0 - pe/pm) > 1.0) {
      std::cout << "\n ERROR: Algo 2: (1.0 - pe/pm) > 1.0" 
		<< " t = " << t << "; R = " << R  
		<<  "; W = " << W << ";\n death = " << death 
		<<  "; growth = " << growth << ";\n pm = " << pm 
		<< "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2:  1 - pe/pm > 1");
    }

    if( (1.0 - pe/pm) < 0.0 ) {
      std::cout << "\n ERROR: Algo 2, (1.0 - pe/pm) < 0.0 "
		<< " t = " << t << "; R = " << R  
		<<  "; W = " << W << ";\n death = " << death 
		<<  "; growth = " << growth << ";\n pm = " << pm 
		<< "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: 1 - pe/pm < 0");
    }

    if( pb > 1.0 ) {
      std::cout << "\n WARNING: Algo 2, pb > 1.0 " 
		<< " t = " << t << "; R = " << R  
		<<  "; W = " << W << ";\n death = " << death 
		<<  "; growth = " << growth << ";\n pm = " << pm 
		<< "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb > 1 ");
    }

    if( pb < 0.0 ) {
      std::cout << "\n WARNING: Algo 2, pb < 0.0 " 
		<< " t = " << t << "; R = " << R  
		<<  "; W = " << W << ";\n death = " << death 
		<<  "; growth = " << growth << ";\n pm = " << pm 
		<< "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb < 0");
    }
    //}


  if( pe == pm ) {
    // Should never happen. Exact identity??
    std::cout << "\n WARNING: Algo 2: pe == pm" << "; t = " << 
      t << "; R = " << R 
      << "; W = " << W << "; death = " << death 
      << "; growth = " << growth << "; pm = " << pm 
      << "; pe = " << pe << std::endl;
    
    return 0.0;
  }

  
  RNGScope scope;
  m = ::Rf_rbinom(num, 1.0 - (pe/pm));
  
  if(m <= 0.5) { // they are integers, so 0 or 1.
    #ifdef DEBUGW // just checking
      if(m != 0.0) 
	std::cout << "\n WARNING: Algo 2: 0.0 < m < 0.5" <<std::endl;
    #endif    
    retval = 0.0;
  } else {
    rnb = ::Rf_rnbinom(m, 1.0 - pb);
    retval = m + rnb;
  }

  if( !std::isfinite(retval) )  {
    std::cout << "\n ERROR: Algo 2, retval not finite "
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb 
	      << "; m = " << m << " ; rnb = " << rnb << std::endl;
    throw std::range_error("Algo 2: retval not finite");

  }

  return retval;
}

double Algo3(const double& num, const double& t, 
	     const double& R, 
	     const double& W, 
	     const double& death, 
	     const double& growth) {
  
  using namespace Rcpp ;


  double pm = pM_f(t, R, W, death, growth);
  double pe = pE_f(pm, W, death, growth);
  double pb = pB_f(pe, death, growth);
  double m; // the holder for the binomial
  double retval;
  double rnb;

  if( (1.0 - pe/pm) > 1.0) {
    std::cout << "\n ERROR: Algo 3: (1.0 - pe/pm) > 1.0" 
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3:  1 - pe/pm > 1");
  }
  
  if( (1.0 - pe/pm) < 0.0 ) {
    std::cout << "\n ERROR: Algo 3, (1.0 - pe/pm) < 0.0 "
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: 1 - pe/pm < 0");
  }
  
  if( pb > 1.0 ) {
    std::cout << "\n WARNING: Algo 3, pb > 1.0 " 
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb > 1 ");
  }
  
  if( pb < 0.0 ) {
    std::cout << "\n WARNING: Algo 3, pb < 0.0 " 
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb < 0");
  }
  
  if( pe == pm ) {
    // Should never happen. Exact identity??
    std::cout << "\n WARNING: Algo 3: pm == pe" << "; t = " << 
      t << "; R = " << R 
	      << "; W = " << W << "; death = " << death 
	      << "; growth = " << growth << "; pm = " << pm 
	      << "; pb = " << pb << std::endl;
    
    return 0.0;
  }

  RNGScope scope;
  m = ::Rf_rbinom(num - 1.0, 1.0 - (pe/pm));
  rnb = ::Rf_rnbinom(m + 2.0, 1.0 - pb);
  retval = m + 1 + rnb;

  if( !std::isfinite(retval) )  {
    std::cout << "\n ERROR: Algo 3, retval not finite "
	      << " t = " << t << "; R = " << R  
	      <<  "; W = " << W << ";\n death = " << death 
	      <<  "; growth = " << growth << ";\n pm = " << pm 
	      << "; pe = " << pe << "; pb = " << pb 
	      << "; m = " << m << " ; rnb = " << rnb << std::endl;
    throw std::range_error("Algo 3: retval not finite");
    
  }

  return retval;
}


double Algo2_st(const spParamsP& spP,
		const double& ti) {

  // beware the use of t: now as it used to be, as we pass the value
  // and take the diff in here: t is the difference

  using namespace Rcpp ;
  double t = ti - spP.timeLastUpdate;
  
  // if (t == 0 ) {
  //   std::cout << "\n Entered Algo2 with t = 0\n" <<
  //     "    Is this a forced sampling case?\n";
  //     return num;
  // }

  if (spP.popSize == 0.0) {
    #ifdef DEBUGW
    std::cout << "\n Entered Algo2 with pop size = 0\n";
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
      std::cout << "\n ERROR: Algo 2: (1.0 - pe/pm) > 1.0\n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2:  1 - pe/pm > 1");
    }

    if( (1.0 - pe/pm) < 0.0 ) {
      std::cout << "\n ERROR: Algo 2, (1.0 - pe/pm) < 0.0 \n";
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: 1 - pe/pm < 0");
    }

    if( pb > 1.0 ) {
      std::cout << "\n WARNING: Algo 2, pb > 1.0 \n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb > 1 ");
    }

    if( pb < 0.0 ) {
      std::cout << "\n WARNING: Algo 2, pb < 0.0 \n"; 
		// << " t = " << t << "; R = " << R  
		// <<  "; W = " << W << ";\n death = " << death 
		// <<  "; growth = " << growth << ";\n pm = " << pm 
		// << "; pe = " << pe << "; pb = " << pb << std::endl;
      throw std::range_error("Algo 2: pb < 0");
    }
    //}


  if( pe == pm ) {
    // Should never happen. Exact identity??
    std::cout << "\n WARNING: Algo 2: pe == pm \n" ;
      // 	      << "; t = " << 
      // t << "; R = " << R 
      // << "; W = " << W << "; death = " << death 
      // << "; growth = " << growth << "; pm = " << pm 
      // << "; pe = " << pe << std::endl;
    return 0.0;
  }

  
  RNGScope scope;
  m = ::Rf_rbinom(spP.popSize, 1.0 - (pe/pm));
  
  if(m <= 0.5) { // they are integers, so 0 or 1.
    #ifdef DEBUGW // just checking
      if(m != 0.0) 
	std::cout << "\n WARNING: Algo 2: 0.0 < m < 0.5" <<std::endl;
    #endif    
    retval = 0.0;
  } else {
    rnb = ::Rf_rnbinom(m, 1.0 - pb);
    retval = m + rnb;
  }

  if( !std::isfinite(retval) )  {
    std::cout << "\n ERROR: Algo 2, retval not finite\n ";
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb 
	      // << "; m = " << m << " ; rnb = " << rnb << std::endl;
    throw std::range_error("Algo 2: retval not finite");
  }
  return retval;
}

double Algo3_st(const spParamsP& spP, const double& t){
  
  using namespace Rcpp ;

  double pm = pM_f(t, spP.R, spP.W, spP.death, spP.birth);
  double pe = pE_f(pm, spP.W, spP.death, spP.birth);
  double pb = pB_f(pe, spP.death, spP.birth);
  double m; // the holder for the binomial
  double retval;
  double rnb;

  if( (1.0 - pe/pm) > 1.0) {
    std::cout << "\n ERROR: Algo 3: (1.0 - pe/pm) > 1.0\n"; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3:  1 - pe/pm > 1");
  }
  
  if( (1.0 - pe/pm) < 0.0 ) {
    std::cout << "\n ERROR: Algo 3, (1.0 - pe/pm) < 0.0\n ";
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: 1 - pe/pm < 0");
  }
  
  if( pb > 1.0 ) {
    std::cout << "\n WARNING: Algo 3, pb > 1.0\n "; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb > 1 ");
  }
  
  if( pb < 0.0 ) {
    std::cout << "\n WARNING: Algo 3, pb < 0.0\n "; 
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb << std::endl;
    throw std::range_error("Algo 3: pb < 0");
  }
  
  if( pe == pm ) {
    // Should never happen. Exact identity??
    std::cout << "\n WARNING: Algo 3: pm == pe\n"; 
    // << "; t = " << 
    // 	 t << "; R = " << R 
    // 	      << "; W = " << W << "; death = " << death 
    // 	      << "; growth = " << growth << "; pm = " << pm 
    // 	      << "; pb = " << pb << std::endl;
    
    return 0.0;
  }

  RNGScope scope;
  m = ::Rf_rbinom(spP.popSize - 1.0, 1.0 - (pe/pm));
  rnb = ::Rf_rnbinom(m + 2.0, 1.0 - pb);
  retval = m + 1 + rnb;

  if( !std::isfinite(retval) )  {
    std::cout << "\n ERROR: Algo 3, retval not finite\n ";
	      // << " t = " << t << "; R = " << R  
	      // <<  "; W = " << W << ";\n death = " << death 
	      // <<  "; growth = " << growth << ";\n pm = " << pm 
	      // << "; pe = " << pe << "; pb = " << pb 
	      // << "; m = " << m << " ; rnb = " << rnb << std::endl;
    throw std::range_error("Algo 3: retval not finite");
  }
  return retval;
}


inline double fitness_linear_std_bitset64(const Genotype64& Genotype,
					  const double& birthRate, 
					  const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype
  // For this approach, death should always be 1, and birth is irrelevant
  int totalMut = 0;
  for(int i = 0; i < numDrivers; ++i) {
    totalMut += Genotype[i];
  }
  return birthRate + s * static_cast<double>(totalMut);
}


inline double fitness_beeren_std_bitset64(const Genotype64& Genotype,
					  //const double& birthRate, 
					  const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype
  // An exponential-like, inspired in Beerenwinkel et al., 2007
  int totalMut = 0;
  for(int i = 0; i < numDrivers; ++i) {
    totalMut += Genotype[i];
  }
  return pow( 1 + s, totalMut );
}



inline double fitness_bozic_bitset64(const Genotype64& Genotype,
				   const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype

  int totalMut = 0;
  for(int i = 0; i < numDrivers; ++i) {
    totalMut += Genotype[i];
  }
  return 0.5 * pow( 1.0 - s, totalMut); 
  // return 1.0 - 0.5 * pow( 1.0 - s, totalMut); 

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


// using std vector for genotypes
// but two different typeFitness
double fitness_CBN_bitset64(const int& mutatedPos, 
			    Rcpp::IntegerMatrix restrictTable,
			    const double& fitnessParent, 
			    const std::string& typeCBN,
			    const Genotype64& Genotype,
			    const double& birthRate, 
			    const double& s, 
			    const int& numDrivers,
			    const std::string typeFitness) {
    
  using namespace Rcpp ;

  // IMPORTANT: remember that for Bozic we return death, o.w. it is birth.

  // FIXME: check numDrivers == ncol(restrictTable)??
  // or obtain the numDrivers from there??
  // Yes, numDrivers is known before hand. Take if from there
  // in some upper level function that calls this one.


  // FIXME: add checks of no negative numbers in deps: thisRestrict?

  // In fact, I could simplify, by having the restrictTable only
  // for drivers that DO depend. Would make it faster and exclude
  // a check below. But less clear.

  
  // will later become an argument
  // double fitnessNo = 0.0;

  int numDependencies;
  int sumPresent = 0;
  int outFitnessYes = 0;

  // remember positions start at 0!!
  if(mutatedPos >= numDrivers) { //the new mutation is a passenger
    return fitnessParent;
    // If I did not pass fitnessParent then
    //  return(fitnessYes(genotype, birth.rate, s, num.drivers))
  } else {
    const Rcpp::IntegerMatrix::Column thisRestrict = restrictTable(_, mutatedPos);
    numDependencies = thisRestrict[1];
    #ifdef DEBUGW
      if(thisRestrict[0] != mutatedPos ) {
	std::cout << std::endl << "thisRestrict[0] = " << thisRestrict[0] 
		<< "; mutatedPos  = " << mutatedPos  << std::endl;
      throw std::range_error("FitnessCBN: thisRestrict[0] != mutatedPos ");	
      }
      if(Genotype[mutatedPos] != 1){
	std::cout << " mutatedPos = " << mutatedPos << std::endl;
	throw std::range_error("FitnessCBN: genotype(mutatedPos) != 1");
      }
    #endif
    

    if(!numDependencies) {
      // if(DEBUG2) {
      // 	std::cout << " FitnessCBN: exiting at !numDependencies " << std::endl;
      // }
      outFitnessYes = 1;
    } else {
      //outFitnessYes = 0;
      for(int i = 2; i < (2 + numDependencies); i++) {
	// if(DEBUG3) {
	//   if(thisRestrict[i] < 0 ) throw std::out_of_range("restict < 0");
	//   }
	sumPresent += Genotype[ thisRestrict[i] ];
      }
      if(typeCBN == "Multiple") {
        if(sumPresent) outFitnessYes = 1;
      } else{ // if(typeCBN == "CBN")
        if(sumPresent == numDependencies)
          outFitnessYes = 1; // return(fitnessYes(genotype, birth.rate, s, num.drivers))
      }
    }
    // FIXME: this could be much improved!
    // There is a part, above, for dependencies, and a part, below,
    // for fitness.
    if(outFitnessYes) {
      if(typeFitness == "bozic")
	return fitness_bozic_bitset64(Genotype, s, numDrivers);
      else if (typeFitness == "beeren")
	return fitness_beeren_std_bitset64(Genotype, s, numDrivers);
      else  
	return fitness_linear_std_bitset64(Genotype, birthRate, s, numDrivers);
    } else {
      if(typeFitness == "bozic")
	return 1.0; //death
      else
	return 0.0; //birth
    }
  }
}

// limited benchmarks suggest the following is slower
inline void new_sp_bitset2(unsigned int& sp, const Genotype64& newGenotype,
			   const std::vector<Genotype64>& Genotypes) {
  sp = std::distance(Genotypes.begin(),
		     std::find(Genotypes.begin(), 
			       Genotypes.end(), newGenotype));
}


inline void new_sp_bitset(unsigned int& sp, const Genotype64& newGenotype,
			  const std::vector<Genotype64>& Genotypes) {
  sp = 0;

  for(sp = 0; sp < Genotypes.size(); ++sp) {
    if( newGenotype == Genotypes[sp] )
      break;
  }
}



void getMutatedPos_bitset(int& mutatedPos, int& numMutablePosParent,
			  gsl_rng *r,
			  std::vector<int>& mutablePos,
			  const Genotype64& nextMutantGenotype,
			  // const int& nextMutant,
			  // const std::vector<Genotype64>& Genotypes,
			  const int& numGenes) {
  // We want mutatedPos and numMutablePosParent turned into a function for
  // profiling 

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
    mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePosParent)];
  } else {
    mutatedPos = mutablePos[0];
  } 


#ifdef DEBUGV
      std::cout << "\n numMutablePosParent = " << numMutablePosParent;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
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


inline void mapTimes_update(std::multimap<double, int>& mapTimes,
		     std::vector<spParamsM>& popParams,
		     const int index,
		     const double time) {
  if(popParams[index].timeLastUpdate > -1)
    mapTimes.erase(popParams[index].pv);
  popParams[index].pv = mapTimes.insert(std::make_pair(time, index));
}

inline void mapTimes_updateP(std::multimap<double, int>& mapTimes,
			     std::vector<spParamsP>& popParams,
			     const int index,
			     const double time) {
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




inline bool popParamsK_time_less (const spParamsK& a, const spParamsK& b) {
  return a.nextMutationTime < b.nextMutationTime;
}



void getMinNextMutationTime3(int& nextMutant, double& minNextMutationTime,
			     const std::vector<spParamsK>& popParams) {
  // we want minNextMutationTime and nextMutant
  // does not work if species with popSize == 0

  std::vector<spParamsK>::const_iterator pbeg = popParams.begin();
  std::vector<spParamsK>::const_iterator pend = popParams.end();

  std::vector<spParamsK>::const_iterator pt_pos_min =
    std::min_element(pbeg, pend, popParamsK_time_less);

  nextMutant = std::distance(pbeg, pt_pos_min);
  minNextMutationTime = popParams[nextMutant].nextMutationTime;
}

// FIXME00: do not use const_iterator?
void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::vector<double>& nextMutTime) {
  // we want minNextMutationTime and nextMutant
  // does not work if species with popSize == 0

  std::vector<double>::const_iterator pbeg = nextMutTime.begin();
  std::vector<double>::const_iterator pend = nextMutTime.end();

  std::vector<double>::const_iterator pt_pos_min =
    std::min_element(pbeg, pend);

  nextMutant = std::distance(pbeg, pt_pos_min);
  minNextMutationTime = nextMutTime[nextMutant];
}


void remove_zero_sp_v4(const std::vector<int>& sp_to_remove,
		       std::vector<Genotype64>& Genotypes,
		       std::vector<spParamsK>& popParams) {

  std::vector<spParamsK>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    popParams.erase(popParams_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
  }
}




void remove_zero_sp_v6(const std::vector<int>& sp_to_remove,
		       std::vector<Genotype64>& Genotypes,
		       std::vector<spParamsM>& popParams,
		       std::multimap<double, int>& mapTimes) {

  std::vector<spParamsM>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    mapTimes.erase(popParams[sp_to_remove[j]].pv);
    popParams.erase(popParams_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
  }
}

void remove_zero_sp_v7(const std::vector<int>& sp_to_remove,
		       std::vector<Genotype64>& Genotypes,
		       std::vector<spParamsP>& popParams,
		       std::multimap<double, int>& mapTimes) {

  std::vector<spParamsP>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    mapTimes.erase(popParams[sp_to_remove[j]].pv);
    popParams.erase(popParams_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
  }
}


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





// this should be templated
void totPopSize_and_fill_out_crude(int& outNS_i,
				   double& totPopSize, 
				   std::vector<Genotype64>& genot_out,
				   //std::vector<unsigned long>& sp_id_out,
				   std::vector<double>& popSizes_out,
				   std::vector<int>& index_out,
				   std::vector<double>& time_out,
				   const std::vector<Genotype64>& Genotypes,
				   const std::vector<spParamsK>& popParams, 
				   const double& currentTime) {
  // Actually, fill out, but also compute totPopSize

  #ifdef DEBUGV
      std::cout << "\n Filling up out crude \n";
#endif
      outNS_i++;
      totPopSize = 0.0;

      time_out.push_back(currentTime);

      for(unsigned int i = 0; i < popParams.size(); ++i) {
	genot_out.push_back(Genotypes[i]);
	popSizes_out.push_back(popParams[i].popSize);
	index_out.push_back(outNS_i);

	totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	std::cout << "\n       Species " << i 
		  << ", Genotype = " << Genotypes[i]
		  << ", sp_id = " << Genotypes[i].to_ulong()
		  << ". Pop size = " << popParams[i].popSize ;
#endif
      }
      
#ifdef DEBUGV
      std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
      
      if( !std::isfinite(totPopSize) ) {
	throw std::range_error("totPopSize not finite");
      }
}

void totPopSize_and_fill_out_crude_P(int& outNS_i,
				   double& totPopSize, 
				   std::vector<Genotype64>& genot_out,
				   //std::vector<unsigned long>& sp_id_out,
				   std::vector<double>& popSizes_out,
				   std::vector<int>& index_out,
				   std::vector<double>& time_out,
				   const std::vector<Genotype64>& Genotypes,
				   const std::vector<spParamsP>& popParams, 
				   const double& currentTime) {
  // Actually, fill out, but also compute totPopSize

  #ifdef DEBUGV
      std::cout << "\n Filling up out crude \n";
#endif
      outNS_i++;
      totPopSize = 0.0;

      time_out.push_back(currentTime);

      for(unsigned int i = 0; i < popParams.size(); ++i) {
	genot_out.push_back(Genotypes[i]);
	popSizes_out.push_back(popParams[i].popSize);
	index_out.push_back(outNS_i);

	totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	std::cout << "\n       Species " << i 
		  << ", Genotype = " << Genotypes[i]
		  << ", sp_id = " << Genotypes[i].to_ulong()
		  << ". Pop size = " << popParams[i].popSize ;
#endif
      }
      
#ifdef DEBUGV
      std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
      
      if( !std::isfinite(totPopSize) ) {
	throw std::range_error("totPopSize not finite");
      }
}


void totPopSize_and_fill_out_crude_M(int& outNS_i,
				   double& totPopSize, 
				   std::vector<Genotype64>& genot_out,
				   //std::vector<unsigned long>& sp_id_out,
				   std::vector<double>& popSizes_out,
				   std::vector<int>& index_out,
				   std::vector<double>& time_out,
				   const std::vector<Genotype64>& Genotypes,
				   const std::vector<spParamsM>& popParams, 
				   const double& currentTime) {
  // Actually, fill out, but also compute totPopSize

  #ifdef DEBUGV
      std::cout << "\n Filling up out crude \n";
#endif
      outNS_i++;
      totPopSize = 0.0;

      time_out.push_back(currentTime);

      for(unsigned int i = 0; i < popParams.size(); ++i) {
	genot_out.push_back(Genotypes[i]);
	popSizes_out.push_back(popParams[i].popSize);
	index_out.push_back(outNS_i);

	totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	std::cout << "\n       Species " << i 
		  << ", Genotype = " << Genotypes[i]
		  << ", sp_id = " << Genotypes[i].to_ulong()
		  << ". Pop size = " << popParams[i].popSize ;
#endif
      }
      
#ifdef DEBUGV
      std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
      
      if( !std::isfinite(totPopSize) ) {
	throw std::range_error("totPopSize not finite");
      }
}
 
//   //FIXME00: try to make const as many as possible 
// inline void reshape_to_outNS(Rcpp::NumericMatrix& outNS,
// 			     std::vector<Genotype64>& uniqueGenotV,
// 			     std::vector<Genotype64>& genot_out,
// 			     std::vector<double>& popSizes_out,
// 			     std::vector<int>& index_out,
// 			     std::vector<double>& time_out){
  
//   std::vector<Genotype64>::iterator ff;
//   std::vector<Genotype64>::iterator fbeg = uniqueGenotV.begin();
//   std::vector<Genotype64>::iterator fend = uniqueGenotV.end();

//   int row;
  
//   for(int i = 0; i < genot_out.size(); ++i) {
//     ff = lower_bound(fbeg, fend, genot_out[i]);
//     row = std::distance(fbeg, ff);
//     // row = std::distance(fbeg, 
//     // 			lower_bound(fbeg, fend, genot_out[i]));
			
//     // row is the species, column the time period
//     outNS(row + 1, index_out[i]) =  popSizes_out[i];
//   }

//   for(int j = 0; j < time_out.size(); ++j)
//     outNS(1, j) = time_out[j];

// }


inline void reshape_to_outNS(Rcpp::NumericMatrix& outNS,
			     const std::vector<unsigned long>& uniqueGenotV,
			     const std::vector<unsigned long>& genot_out_ul,
			     const std::vector<double>& popSizes_out,
			     const std::vector<int>& index_out,
			     const std::vector<double>& time_out){
  
  std::vector<unsigned long>::const_iterator fbeg = uniqueGenotV.begin();
  std::vector<unsigned long>::const_iterator fend = uniqueGenotV.end();

  int row;

  for(unsigned int i = 0; i < genot_out_ul.size(); ++i) {
    row = std::distance(fbeg, lower_bound(fbeg, fend, genot_out_ul[i]) );
    // row is the species, column the time period
    outNS(row + 1, index_out[i]) =  popSizes_out[i];
  }

  for(unsigned int j = 0; j < time_out.size(); ++j)
    outNS(0, j) = time_out[j];

}



// does not work, as "<" not defined on bitset.
// see possible solution here http://en.allexperts.com/q/C-1040/inserting-bitset-set.htm
// but I'd rather make conversion explicit and under control
// inline void find_unique_genotypes_v0(std::set<Genotype64>& uniqueGenotypes,
// 				  const std::vector<Genotype64>& genot_out) {
//   for(unsigned i = 0; i < genot_out.size(); ++i) 
//     uniqueGenotypes.insert( genot_out[i] );
// }

inline void find_unique_genotypes(std::set<unsigned long>& uniqueGenotypes,
				  const std::vector<unsigned long>& genot_out_l) {
  for(unsigned i = 0; i < genot_out_l.size(); ++i) 
    uniqueGenotypes.insert( genot_out_l[i] );
}

inline void genot_out_to_ulong(std::vector<unsigned long>& go_l,
			       const std::vector<Genotype64>& go) {
  for(unsigned i = 0; i < go.size(); ++i)
    go_l[i] = go[i].to_ulong();
}

// inline void find_unique_genotypes2(std::set<Genotype64>& uniqueGenotypes,
// 				   int& totalNumSpecies,
// 				   const std::vector<Genotype64>& genot_out) {
//   static unsigned last_g_i = 0;

//   for(last_g_i; last_g_i < genot_out.size(); ++last_g_i) 
//     uniqueGenotypes.insert( genot_out[last_g_i] );

//   totalNumSpecies = uniqueGenotypes.size();
// }

// inline void uniqueGenotypes_to_vector(std::vector<Genotype64>& ugV,
// 				   const std::set<Genotype64>& uniqueGenotypes) {
//   ugV.assign(uniqueGenotypes.begin(), uniqueGenotypes.end() );
// }

inline void uniqueGenotypes_to_vector(std::vector<unsigned long>& ugV,
				      const std::set<unsigned long>& uniqueGenotypes) {
  ugV.assign(uniqueGenotypes.begin(), uniqueGenotypes.end() );
}


inline void create_returnGenotypes(Rcpp::IntegerMatrix& returnGenotypes,
				   const int& numGenes,
				   const std::vector<unsigned long>& uniqueGenotypesV){
  
#ifdef DEBUGV
  std::cout << "\n Here X4 \n";
#endif  

  for(unsigned int i = 0; i < uniqueGenotypesV.size(); ++i) {
    Genotype64 tmpbs(uniqueGenotypesV[i]);
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = tmpbs[j];
    }
  }

#ifdef DEBUGV
  std::cout << "\n Here X5 \n";
#endif  


}


// void create_returnGenotypes(Rcpp::IntegerMatrix& returnGenotypes,
// 			    const int& numGenes,
// 			    const std::vector<Genotype64>& uniqueGenotypesV){
  
//   for(int i = 0; i < uniqueGenotypesV.size(); ++i) {
//     for(int j = 0; j < numGenes; ++j) {
//       returnGenotypes(i, j) = uniqueGenotypesV[i][j];
//     }
//   }
// }



// template this too
void sample_all_pop_K(double& currentTime, 
		      std::vector<int>& sp_to_remove,
		      std::vector<spParamsK>& popParams,
		      const std::vector<Genotype64>& Genotypes,
		      const double& tSample){

  currentTime = tSample;
  sp_to_remove[0] = 0;

  for(unsigned int i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].Flag == false);
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    std::cout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
	      << "\n      Genotype = " << Genotypes[i]
	      << "\n      sp_id = " << Genotypes[i].to_ulong() // sp_id[i]  
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
	      << " \n     species birth " << popParams[i].birth
	      << " \n     species nextMutationTime " 
	      << popParams[i].nextMutationTime;
#endif

    // Account for forceSampling. When 
    // forceSampling, popSize for at least one species
    // was updated in previous loop, so we skip that one
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize = 
	Algo2(popParams[i].popSize,
	      tSample - popParams[i].timeLastUpdate,
	      popParams[i].R,
	      popParams[i].W,
	      popParams[i].death,
	      popParams[i].birth);
    }
    if( popParams[i].popSize <=  0.0 ) {
      // this i has never been non-zero in any sampling time
      sp_to_remove[0]++;
      sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
      std::cout << "\n\n     Removing species i = " << i 
		<< " with sp_id = " << Genotypes[i].to_ulong(); //sp_id[i];
#endif
    } else {
      popParams[i].Flag = true;
    }
#ifdef DEBUGV
    std::cout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
#ifdef DEBUGW	  
    popParams[i].timeLastUpdate = -999999.99999;
#endif
    
  }
}

void sample_all_pop_M(double& currentTime, 
		      std::vector<int>& sp_to_remove,
		      std::vector<spParamsM>& popParams,
		      const std::vector<Genotype64>& Genotypes,
		      const double& tSample){

  currentTime = tSample;
  sp_to_remove[0] = 0;

  for(unsigned int i = 0; i < popParams.size(); i++) {
    STOPASSERT(popParams[i].Flag == false);
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    std::cout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
	      << "\n      Genotype = " << Genotypes[i]
	      << "\n      sp_id = " << Genotypes[i].to_ulong() // sp_id[i]  
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
	      << " \n     species birth " << popParams[i].birth
	      << " \n     species nextMutationTime " 
	      << popParams[i].nextMutationTime;
#endif

    // Account for forceSampling. When 
    // forceSampling, popSize for at least one species
    // was updated in previous loop, so we skip that one
    if(tSample > popParams[i].timeLastUpdate) {
      popParams[i].popSize = 
	Algo2(popParams[i].popSize,
	      tSample - popParams[i].timeLastUpdate,
	      popParams[i].R,
	      popParams[i].W,
	      popParams[i].death,
	      popParams[i].birth);
    }
    if( popParams[i].popSize <=  0.0 ) {
      // this i has never been non-zero in any sampling time
      sp_to_remove[0]++;
      sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
      std::cout << "\n\n     Removing species i = " << i 
		<< " with sp_id = " << Genotypes[i].to_ulong(); //sp_id[i];
#endif
    } else {
      popParams[i].Flag = true;
    }
#ifdef DEBUGV
    std::cout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
#ifdef DEBUGW	  
    popParams[i].timeLastUpdate = -999999.99999;
#endif
    
  }
}


void sample_all_pop_P(double& currentTime, 
		      std::vector<int>& sp_to_remove,
		      std::vector<spParamsP>& popParams,
		      const std::vector<Genotype64>& Genotypes,
		      const double& tSample){

  currentTime = tSample;
  sp_to_remove[0] = 0;

  for(unsigned int i = 0; i < popParams.size(); i++) {
    //STOPASSERT(popParams[i].Flag == false);
    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
    STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
    std::cout << "\n\n     ********* 5.9 ******\n " 
	      << "     Species  = " << i 
	      << "\n      Genotype = " << Genotypes[i]
	      << "\n      sp_id = " << Genotypes[i].to_ulong() // sp_id[i]  
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
	      << " \n     species birth " << popParams[i].birth
	      << " \n     species nextMutationTime " 
	      << popParams[i].nextMutationTime;
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
      sp_to_remove[0]++;
      sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
      std::cout << "\n\n     Removing species i = " << i 
		<< " with sp_id = " << Genotypes[i].to_ulong(); //sp_id[i];
#endif
    } // else {
    //   popParams[i].Flag = true;
    // }
#ifdef DEBUGV
    std::cout << "\n\n   post-update popSize = " 
	      << popParams[i].popSize << "\n";
#endif
#ifdef DEBUGW	  
    popParams[i].timeLastUpdate = -999999.99999;
#endif
    
  }
}


SEXP Algorithm5P(SEXP restrictTable_,
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
		 SEXP typeFitness_) {

  // Based on 5O, but functions now take struct, not doubles, etc via []


  BEGIN_RCPP
  using namespace Rcpp;
    
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
  const int initIt = as<int>(initSize_iter_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // species
  const int seed = as<int>(seed_gsl_);

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
  //  int totalNumSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  //int indexMutatedPos = 0;
  int outNS_i = 0; // the column in the outNS
  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  //int timeL = 1000;
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  //FIXME: this is ugly and will eventually break.
  // or turn into regular vector.
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  // FIXME: uncomment this for package
  // verify we are OK with usigned long
  // if( !(static_cast<double>(std::numeric_limits<unsigned long>::max()) 
  // 	>= pow(2, 64)) )
  //   throw std::range_error("The size of unsigned long is too short.");

  // if(numGenes > 64)  
  //   throw std::range_error("This version only accepts up to 64 genes.");



  // those to update
  bool short_update = true;
  int u_1 = -99;
  int u_2 = -99;

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);

  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  //  std::set<Genotype64> uniqueGenotypes;
  std::set<unsigned long> uniqueGenotypes;
  spParamsP tmpParam; 
  std::vector<spParamsP> popParams(1);
  const int sp_per_period = 5000;

  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

  //Output
  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!

  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  // FIXME00: is m1pos needed?
  // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;


  // 5.1 Initialize 
  //popParams[0].Flag = true;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else if (typeFitness == "beeren") {
    popParams[0].birth = 1.0;
    popParams[0].death = death;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  // FIXME: use st versions
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  Genotypes[0].reset();
  genot_out.push_back(Genotypes[0]);
  popSizes_out.push_back(popParams[0].popSize);
  index_out.push_back(outNS_i);
  uniqueGenotypes.insert(Genotypes[0].to_ulong());
  time_out.push_back(currentTime);
  

  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;
  //totalNumSpecies = 1;

  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      std::cout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    std::cout << " DEBUGV\n";
    std::cout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    // FIXME00: but this is crazy! in most loops, I only need
    // to iterate over just two!!! Not the whole population
    // Create a vector of those to update.
    // to_update(3). First is number: 0 -> whole population
    // 2: two, and then give the pos.
    // A problem with forced sampling??

    // The special case of the first iteration: FIXME
    // Convert this into a function, for profiling?

    // FIXME: unroll completely the case of first iteration, and start
    // later.  Is it worth it?

    if(iter == 1) { // handle special case of first iter
      tmpdouble1 = ti_nextTime_tmax_2_st(popParams[0],
					 currentTime,
					 tSample);
      mapTimes_updateP(mapTimes, popParams, 0, tmpdouble1);
      //popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(short_update) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2_st(popParams[u_1],
					   currentTime,
					   tSample);
	mapTimes_updateP(mapTimes, popParams, u_1, tmpdouble1);

	tmpdouble2 = ti_nextTime_tmax_2_st(popParams[u_2],
					   currentTime,
					   tSample);
	mapTimes_updateP(mapTimes, popParams, u_2, tmpdouble2);
	
	// popParams[u_1].Flag = false;
	// popParams[u_2].Flag = false;
	popParams[u_1].timeLastUpdate = currentTime;
	popParams[u_2].timeLastUpdate = currentTime;
      } else { // we sampled, so update all
	for(unsigned int i = 0; i < popParams.size(); i++) {
	  // FIXME: modify function, and pass only an
	  // element of the structure!
	  tmpdouble1 = ti_nextTime_tmax_2_st(popParams[i],
					     currentTime,
					     tSample);
	  mapTimes_updateP(mapTimes, popParams, i, tmpdouble1);
	  // popParams[i].Flag = false;
	  popParams[i].timeLastUpdate = currentTime;
	  
#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
		    << "\n       popParams[i].nextMutationTime = " 
		    << popParams[i].nextMutationTime
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
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      short_update = true;
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
#ifdef DEBUGW
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
#endif

      // FIXME00: pass only structre and other thing. Do not index params.
      popParams[nextMutant].popSize = Algo3_st(popParams[nextMutant],
					       mutantTimeSinceLastUpdate);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
#ifdef DEBUGV
	//if(verbosity > -2) {
	  // We always warn about this, since interaction with ti==0
	  std::cout << "\n Forced sampling triggered for next loop: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant <<
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
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
#endif
      }
      
      // ************   5.5   ***************
      getMutatedPos_bitset(mutatedPos, numMutablePosParent, r, 
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
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	tmpParam.popSize = 1;

	// FIXME00: This is ugly!!
	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN_bitset64(mutatedPos,
					      restrictTable,
					      popParams[nextMutant].death,
					      typeCBN,
					      newGenotype,
					      birthRate,
					      s,
					      numDrivers,
					      typeFitness);
	  tmpParam.birth = 1 - tmpParam.death;
	} else {
	  tmpParam.birth = fitness_CBN_bitset64(mutatedPos,
						  restrictTable,
						  popParams[nextMutant].birth,
						  typeCBN,
						  newGenotype,
						  birthRate,
						  s,
						  numDrivers,
						  typeFitness);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n Genotype = " << Genotypes[sp];
	  std::cout << "\n sp_id = " << Genotypes[sp].to_ulong();
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent Genotypes = " << Genotypes[nextMutant];
	}
#endif
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	tmpParam.mutation = mu * (numMutablePosParent - 1);

	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	}
	// to easily find out if new species
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors

	// tmpParam.Flag = true;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (Genotype = " << Genotypes[sp] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")"
		    << "\n from species "  <<   nextMutant
		    << " (Genotypes = " << Genotypes[nextMutant] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")";
	}

	// FIXME00: the if can be removed??
	      // FIXME00: pass only structre and other thing. Do not index params.

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2_st(popParams[sp], currentTime);
	  if(verbosity >= 2) {
	    std::cout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  }
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	  if(verbosity >= 2) {
	    std::cout << "\n Mutation to an extinct species\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	//popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      // popParams[nextMutant].Flag = true;

      // to make output identical to previous versions
      if(nextMutant < sp) {
      	u_1 = nextMutant;
      	u_2 = sp;
      } else {
      	u_1 = sp;
      	u_2 = nextMutant;
      }

      // u_1 = nextMutant;
      // u_2 = sp;
      // pop of receiving mutant flagged above 
    } else { //       *********** We are sampling **********
      short_update = false;
      if(verbosity >= 2) {
	std::cout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  std::cout << " at time " << tSample << "\n";
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }

      sample_all_pop_P(currentTime, sp_to_remove, 
		       popParams, Genotypes, tSample);

      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove[0])
	remove_zero_sp_v7(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();
      
      totPopSize_and_fill_out_crude_P(outNS_i, totPopSize, genot_out, 
				      //sp_id_out,
				      popSizes_out, index_out,
				      time_out, Genotypes, popParams, 
				      currentTime);

      //find_unique_genotypes2(uniqueGenotypes, totalNumSpecies, genot_out);
      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || (tSample >= finalTime)) 
	simulsDone = true;
	
      forceSample = false;
#ifdef DEBUGV
      std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species


  // FIXME: all this is ugly and could be a single function
  std::vector<unsigned long> genot_out_ulong(genot_out.size());
  genot_out_to_ulong(genot_out_ulong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ulong);
  std::vector<unsigned long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);


  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)

  NumericMatrix outNS(uniqueGenotypes.size() + 1, outNS_i + 1);
  reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ulong, 
		   popSizes_out, 
		   index_out, time_out);
  
#ifdef DEBUGV
  std::cout << "\n Here X004 \n";
#endif  


  IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);

#ifdef DEBUGV
  std::cout << "\n Here X04 \n";
#endif  

  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);
  
#ifdef DEBUGV
  std::cout << "\n Here X6 \n";
#endif  


  return List::create(Named("NumSpecies") = uniqueGenotypes.size(), 
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = outNS,
		      Named("Genotypes") = returnGenotypes,
		      Named("FinalTime") = currentTime,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  END_RCPP
}


SEXP Algorithm5O(SEXP restrictTable_,
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
		 SEXP typeFitness_) {
  // Based on 5M, but loop only over flagged, keeping a different structure


  BEGIN_RCPP
  using namespace Rcpp;
    
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
  const int initIt = as<int>(initSize_iter_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // species
  const int seed = as<int>(seed_gsl_);

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
  //  int totalNumSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  // int indexMutatedPos = 0;
  int outNS_i = 0; // the column in the outNS
  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  //int timeL = 1000;
  
  double tmpdouble1 = 0.0;
  double tmpdouble2 = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  // FIXME: uncomment this for package
  // verify we are OK with usigned long
  // if( !(static_cast<double>(std::numeric_limits<unsigned long>::max()) 
  // 	>= pow(2, 64)) )
  //   throw std::range_error("The size of unsigned long is too short.");

  // if(numGenes > 64)  
  //   throw std::range_error("This version only accepts up to 64 genes.");



  // those to update
  bool short_update = true;
  int u_1 = -99;
  int u_2 = -99;

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);

  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  //  std::set<Genotype64> uniqueGenotypes;
  std::set<unsigned long> uniqueGenotypes;
  spParamsM tmpParam; 
  std::vector<spParamsM> popParams(1);
  const int sp_per_period = 5000;

  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

  //Output
  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!

  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  // FIXME00: is m1pos needed?
  // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;


  // 5.1 Initialize 
  popParams[0].Flag = true;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else if (typeFitness == "beeren") {
    popParams[0].birth = 1.0;
    popParams[0].death = death;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  Genotypes[0].reset();
  genot_out.push_back(Genotypes[0]);
  popSizes_out.push_back(popParams[0].popSize);
  index_out.push_back(outNS_i);
  uniqueGenotypes.insert(Genotypes[0].to_ulong());
  time_out.push_back(currentTime);
  

  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;
  //totalNumSpecies = 1;

  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      std::cout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    std::cout << " DEBUGV\n";
    std::cout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    // FIXME00: but this is crazy! in most loops, I only need
    // to iterate over just two!!! Not the whole population
    // Create a vector of those to update.
    // to_update(3). First is number: 0 -> whole population
    // 2: two, and then give the pos.
    // A problem with forced sampling??

    // The special case of the first iteration: FIXME
    // Convert this into a function, for profiling?

    // FIXME: unroll completely the case of first iteration, and start
    // later.  Is it worth it?

    if(iter == 1) { // handle special case of first iter
      tmpdouble1 = ti_nextTime_tmax_2(popParams[0].R,
				      popParams[0].W,
				      popParams[0].death,
				      popParams[0].birth,
				      popParams[0].mutation,
				      popParams[0].popSize,
				      currentTime,
				      tSample);
      mapTimes_update(mapTimes, popParams, 0, tmpdouble1);
      popParams[0].Flag = false;
      popParams[0].timeLastUpdate = currentTime;
    } else { // any other iter
      if(short_update) {
	// we did not sample in previous period.
	tmpdouble1 = ti_nextTime_tmax_2(popParams[u_1].R,
					popParams[u_1].W,
					popParams[u_1].death,
					popParams[u_1].birth,
					popParams[u_1].mutation,
					popParams[u_1].popSize,
					currentTime,
					tSample);
	mapTimes_update(mapTimes, popParams, u_1, tmpdouble1);

	tmpdouble2 = ti_nextTime_tmax_2(popParams[u_2].R,
					popParams[u_2].W,
					popParams[u_2].death,
					popParams[u_2].birth,
					popParams[u_2].mutation,
					popParams[u_2].popSize,
					currentTime,
					tSample);
	mapTimes_update(mapTimes, popParams, u_2, tmpdouble2);
	
	popParams[u_1].Flag = false;
	popParams[u_2].Flag = false;
	popParams[u_1].timeLastUpdate = currentTime;
	popParams[u_2].timeLastUpdate = currentTime;
      } else { // we sampled, so update all
	for(unsigned int i = 0; i < popParams.size(); i++) {
	  // FIXME: modify function, and pass only an
	  // element of the structure!
	  tmpdouble1 = ti_nextTime_tmax_2(popParams[i].R,
					  popParams[i].W,
					  popParams[i].death,
					  popParams[i].birth,
					  popParams[i].mutation,
					  popParams[i].popSize,
					  currentTime,
					  tSample);
	  mapTimes_update(mapTimes, popParams, i, tmpdouble1);
	  popParams[i].Flag = false;
	  popParams[i].timeLastUpdate = currentTime;
	  
#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
		    << "\n       popParams[i].nextMutationTime = " 
		    << popParams[i].nextMutationTime
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
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      short_update = true;
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
#ifdef DEBUGW
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
#endif

      // FIXME00: pass only structre and other thing. Do not index params.
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    popParams[nextMutant].death,
					    popParams[nextMutant].birth);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
	if(verbosity > -2) {
	  // We always warn about this, since interaction with ti==0
	  std::cout << "\n Forced sampling triggered for next loop: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant <<
	    " at iteration " << iter << "\n";
	}
      }      
      
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet 
      // now total number of different species
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      getMutatedPos_bitset(mutatedPos, numMutablePosParent, r, 
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
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	tmpParam.popSize = 1;

	// FIXME00: This is ugly!!
	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN_bitset64(mutatedPos,
					      restrictTable,
					      popParams[nextMutant].death,
					      typeCBN,
					      newGenotype,
					      birthRate,
					      s,
					      numDrivers,
					      typeFitness);
	  tmpParam.birth = 1 - tmpParam.death;
	} else {
	  tmpParam.birth = fitness_CBN_bitset64(mutatedPos,
						  restrictTable,
						  popParams[nextMutant].birth,
						  typeCBN,
						  newGenotype,
						  birthRate,
						  s,
						  numDrivers,
						  typeFitness);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n Genotype = " << Genotypes[sp];
	  std::cout << "\n sp_id = " << Genotypes[sp].to_ulong();
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent Genotypes = " << Genotypes[nextMutant];
	}
#endif
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	tmpParam.mutation = mu * (numMutablePosParent - 1);

	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	}
	// to easily find out if new species
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors

	tmpParam.Flag = true;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (Genotype = " << Genotypes[sp] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")"
		    << "\n from species "  <<   nextMutant
		    << " (Genotypes = " << Genotypes[nextMutant] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")";
	}

	// FIXME00: the if can be removed??
	      // FIXME00: pass only structre and other thing. Do not index params.

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  popParams[sp].death,
		  popParams[sp].birth);
	  if(verbosity >= 2) {
	    std::cout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  }
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	  if(verbosity >= 2) {
	    std::cout << "\n Mutation to an extinct species\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      popParams[nextMutant].Flag = true;
      // to make output identical to previous versions
      if(nextMutant < sp) {
      	u_1 = nextMutant;
      	u_2 = sp;
      } else {
      	u_1 = sp;
      	u_2 = nextMutant;
      }
      // u_1 = nextMutant;
      // u_2 = sp;
      // pop of receiving mutant flagged above 
    } else { //       *********** We are sampling **********
      short_update = false;
      if(verbosity >= 2) {
	std::cout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  std::cout << " at time " << tSample << "\n";
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }

      sample_all_pop_M(currentTime, sp_to_remove, 
		       popParams, Genotypes, tSample);

      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove[0])
	remove_zero_sp_v6(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();
      
      totPopSize_and_fill_out_crude_M(outNS_i, totPopSize, genot_out, 
				      //sp_id_out,
				      popSizes_out, index_out,
				      time_out, Genotypes, popParams, 
				      currentTime);

      //find_unique_genotypes2(uniqueGenotypes, totalNumSpecies, genot_out);
      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || (tSample >= finalTime)) 
	simulsDone = true;
	
      forceSample = false;
#ifdef DEBUGV
      std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species


  // FIXME: all this is ugly and could be a single function
  std::vector<unsigned long> genot_out_ulong(genot_out.size());
  genot_out_to_ulong(genot_out_ulong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ulong);
  std::vector<unsigned long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);


  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)

  NumericMatrix outNS(uniqueGenotypes.size() + 1, outNS_i + 1);
  reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ulong, 
		   popSizes_out, 
		   index_out, time_out);
  
#ifdef DEBUGV
  std::cout << "\n Here X004 \n";
#endif  


  IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);

#ifdef DEBUGV
  std::cout << "\n Here X04 \n";
#endif  

  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);
  
#ifdef DEBUGV
  std::cout << "\n Here X6 \n";
#endif  


  return List::create(Named("NumSpecies") = uniqueGenotypes.size(), 
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = outNS,
		      Named("Genotypes") = returnGenotypes,
		      Named("FinalTime") = currentTime,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  END_RCPP
}



SEXP Algorithm5M(SEXP restrictTable_,
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
		 SEXP typeFitness_) {
  // Based on 5K, but time is kept ordered using multimap


  BEGIN_RCPP
  using namespace Rcpp;
    
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
  const int initIt = as<int>(initSize_iter_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // species
  const int seed = as<int>(seed_gsl_);

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
  //  int totalNumSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  // int indexMutatedPos = 0;
  int outNS_i = 0; // the column in the outNS
  unsigned int sp = 0;
  // int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  // int timeL = 1000;
  
  double tmpdouble = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  // FIXME: uncomment this for package
  // verify we are OK with usigned long
  // if( !(static_cast<double>(std::numeric_limits<unsigned long>::max()) 
  // 	>= pow(2, 64)) )
  //   throw std::range_error("The size of unsigned long is too short.");

  // if(numGenes > 64)  
  //   throw std::range_error("This version only accepts up to 64 genes.");


  
  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);

  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  //  std::set<Genotype64> uniqueGenotypes;
  std::set<unsigned long> uniqueGenotypes;
  spParamsM tmpParam; 
  std::vector<spParamsM> popParams(1);
  const int sp_per_period = 5000;

  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

  //Output
  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!

  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);

  // FIXME00: is m1pos needed?
  // multimap to hold nextMutationTime
  std::multimap<double, int> mapTimes;
  //std::multimap<double, int>::iterator m1pos;


  // 5.1 Initialize 
  popParams[0].Flag = true;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else if (typeFitness == "beeren") {
    popParams[0].birth = 1.0;
    popParams[0].death = death;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);

  popParams[0].pv = mapTimes.insert(std::make_pair(-999, 0));

  Genotypes[0].reset();
  genot_out.push_back(Genotypes[0]);
  popSizes_out.push_back(popParams[0].popSize);
  index_out.push_back(outNS_i);
  uniqueGenotypes.insert(Genotypes[0].to_ulong());
  time_out.push_back(currentTime);
  

  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;
  //totalNumSpecies = 1;

  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      std::cout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    std::cout << " DEBUGV\n";
    std::cout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    // FIXME00: but this is crazy! in most loops, I only need
    // to iterate over just two!!! Not the whole population
    // Create a vector of those to update.
    // to_update(3). First is number: 0 -> whole population
    // 2: two, and then give the pos.
    // A problem with forced sampling??

    // Convert this into a function, for profiling?
    for(unsigned int i = 0; i < popParams.size(); i++) {
      if( popParams[i].Flag )  {
	tmpdouble = ti_nextTime_tmax_2(popParams[i].R,
				       popParams[i].W,
				       popParams[i].death,
				       popParams[i].birth,
				       popParams[i].mutation,
				       popParams[i].popSize,
				       currentTime,
				       tSample);
	mapTimes_update(mapTimes, popParams, i, tmpdouble);
	popParams[i].Flag = false;
	popParams[i].timeLastUpdate = currentTime;

#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
		    << "\n       popParams[i].nextMutationTime = " 
		    << popParams[i].nextMutationTime
		    << " \n     species R " << popParams[i].R
		    << " \n     species W " << popParams[i].W
		    << " \n     species death " << popParams[i].death
		    << " \n     species birth " << popParams[i].birth;
#endif
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
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
#ifdef DEBUGW
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
#endif
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    popParams[nextMutant].death,
					    popParams[nextMutant].birth);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
	if(verbosity > -2) {
	  // We always warn about this, since interaction with ti==0
	  std::cout << "\n Forced sampling triggered for next loop: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant <<
	    " at iteration " << iter << "\n";
	}
      }      
      
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet 
      // now total number of different species
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      getMutatedPos_bitset(mutatedPos, numMutablePosParent, r, 
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
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	tmpParam.popSize = 1;

	// FIXME00: This is ugly!!
	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN_bitset64(mutatedPos,
					      restrictTable,
					      popParams[nextMutant].death,
					      typeCBN,
					      newGenotype,
					      birthRate,
					      s,
					      numDrivers,
					      typeFitness);
	  tmpParam.birth = 1 - tmpParam.death;
	} else {
	  tmpParam.birth = fitness_CBN_bitset64(mutatedPos,
						  restrictTable,
						  popParams[nextMutant].birth,
						  typeCBN,
						  newGenotype,
						  birthRate,
						  s,
						  numDrivers,
						  typeFitness);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n Genotype = " << Genotypes[sp];
	  std::cout << "\n sp_id = " << Genotypes[sp].to_ulong();
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent Genotypes = " << Genotypes[nextMutant];
	}
#endif
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	tmpParam.mutation = mu * (numMutablePosParent - 1);

	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	}
	// to easily find out if new species
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors

	tmpParam.Flag = true;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (Genotype = " << Genotypes[sp] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")"
		    << "\n from species "  <<   nextMutant
		    << " (Genotypes = " << Genotypes[nextMutant] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")";
	}

	// FIXME00: the if can be removed??
	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  popParams[sp].death,
		  popParams[sp].birth);
	  if(verbosity >= 2) {
	    std::cout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  }
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	  if(verbosity >= 2) {
	    std::cout << "\n Mutation to an extinct species\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      popParams[nextMutant].Flag = true;
      // pop of receiving mutant flagged above 

    } else { //       *********** We are sampling **********
      if(verbosity >= 2) {
	std::cout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  std::cout << " at time " << tSample << "\n";
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }

      sample_all_pop_M(currentTime, sp_to_remove, 
		       popParams, Genotypes, tSample);

      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove[0])
	remove_zero_sp_v6(sp_to_remove, Genotypes, popParams, mapTimes);

      numSpecies = popParams.size();
      
      totPopSize_and_fill_out_crude_M(outNS_i, totPopSize, genot_out, 
				      //sp_id_out,
				      popSizes_out, index_out,
				      time_out, Genotypes, popParams, 
				      currentTime);

      //find_unique_genotypes2(uniqueGenotypes, totalNumSpecies, genot_out);
      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || (tSample >= finalTime)) 
	simulsDone = true;
	
      forceSample = false;
#ifdef DEBUGV
      std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species


  // FIXME: all this is ugly and could be a single function
  std::vector<unsigned long> genot_out_ulong(genot_out.size());
  genot_out_to_ulong(genot_out_ulong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ulong);
  std::vector<unsigned long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);


  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)

  NumericMatrix outNS(uniqueGenotypes.size() + 1, outNS_i + 1);
  reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ulong, 
		   popSizes_out, 
		   index_out, time_out);
  
#ifdef DEBUGV
  std::cout << "\n Here X004 \n";
#endif  


  IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);

#ifdef DEBUGV
  std::cout << "\n Here X04 \n";
#endif  

  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);
  
#ifdef DEBUGV
  std::cout << "\n Here X6 \n";
#endif  


  return List::create(Named("NumSpecies") = uniqueGenotypes.size(), 
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = outNS,
		      Named("Genotypes") = returnGenotypes,
		      Named("FinalTime") = currentTime,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  END_RCPP
}


SEXP Algorithm5K(SEXP restrictTable_,
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
		 SEXP typeFitness_) {
  // Based on 5J, using bitsets and removing all 
  // those species that become zero at any sampling.


  BEGIN_RCPP
  using namespace Rcpp;
    
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
  const int initIt = as<int>(initSize_iter_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  double ratioForce = as<double>(ratioForce_); // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  int speciesFS = as<int>(speciesFS_); // to force sampling when too many 
  // species
  const int seed = as<int>(seed_gsl_);

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
  //  int totalNumSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  //int indexMutatedPos = 0;
  int outNS_i = 0; // the column in the outNS
  unsigned int sp = 0;
  //int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  //int timeL = 1000;
  
  // double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  // FIXME: uncomment this for package
  // verify we are OK with usigned long
  // if( !(static_cast<double>(std::numeric_limits<unsigned long>::max()) 
  // 	>= pow(2, 64)) )
  //   throw std::range_error("The size of unsigned long is too short.");

  // if(numGenes > 64)  
  //   throw std::range_error("This version only accepts up to 64 genes.");


  
  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);

  Genotype64 newGenotype;
  std::vector<Genotype64> Genotypes(1);
  //  std::set<Genotype64> uniqueGenotypes;
  std::set<unsigned long> uniqueGenotypes;
  spParamsK tmpParam; 
  std::vector<spParamsK> popParams(1);
  const int sp_per_period = 5000;

  popParams.reserve(sp_per_period);
  Genotypes.reserve(sp_per_period);

  std::vector<int>mutablePos(numGenes); // could be inside getMuatedPos_bitset

  //Output
  std::vector<Genotype64> genot_out;
  std::vector<double> popSizes_out;
  std::vector<int> index_out; 
  std::vector<double> time_out; //only one entry per period!

  genot_out.reserve(initSp);
  popSizes_out.reserve(initSp);
  index_out.reserve(initSp);
  time_out.reserve(initIt);


  // 5.1 Initialize 
  popParams[0].Flag = true;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else if (typeFitness == "beeren") {
    popParams[0].birth = 1.0;
    popParams[0].death = death;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, 
		       popParams[0].mutation);

  Genotypes[0].reset();
  genot_out.push_back(Genotypes[0]);
  popSizes_out.push_back(popParams[0].popSize);
  index_out.push_back(outNS_i);
  uniqueGenotypes.insert(Genotypes[0].to_ulong());
  time_out.push_back(currentTime);

  timeNextPopSample = currentTime + sampleEvery;
  numSpecies = 1;
  //totalNumSpecies = 1;

  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL )) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
      }
    }

    //  ************   5.2   ***************
    if(verbosity >= 2)
      std::cout <<"\n\n\n*** Looping through 5.2. Iter = " << iter << " \n";

    tSample = std::min(timeNextPopSample, finalTime);

#ifdef DEBUGV
    std::cout << " DEBUGV\n";
    std::cout << "\n ForceSample? " << forceSample 
	      << "  tSample " << tSample 
	      << "  currentTime " << currentTime;
#endif

    
    for(unsigned int i = 0; i < popParams.size(); i++) {
      if( popParams[i].Flag )  {
	popParams[i].nextMutationTime = ti_nextTime_tmax_2(popParams[i].R,
							   popParams[i].W,
							   popParams[i].death,
							   popParams[i].birth,
							   popParams[i].mutation,
							   popParams[i].popSize,
							   currentTime,
							   tSample);
	popParams[i].Flag = false;
	popParams[i].timeLastUpdate = currentTime;

#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " 
		    << "     Species  = " << i 
		    << "\n       genotype =  " << Genotypes[i] 
		    << "\n       popSize = " << popParams[i].popSize 
		    << "\n       currentTime = " << currentTime 
		    << "\n       popParams[i].nextMutationTime = " 
		    << popParams[i].nextMutationTime
		    << " \n     species R " << popParams[i].R
		    << " \n     species W " << popParams[i].W
		    << " \n     species death " << popParams[i].death
		    << " \n     species birth " << popParams[i].birth;
#endif
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
    

   
    getMinNextMutationTime3(nextMutant, minNextMutationTime, popParams);
 
    if(verbosity >= 2) {
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;

      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
#ifdef DEBUGW
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
#endif
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    popParams[nextMutant].death,
					    popParams[nextMutant].birth);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
	if(verbosity > -2) {
	  // We always warn about this, since interaction with ti==0
	  std::cout << "\n Forced sampling triggered for next loop: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant <<
	    " at iteration " << iter << "\n";
	}
      }      
      
      // Check also for numSpecies, and force sampling if needed
      // This is very different from the other algos, as we do not yet 
      // now total number of different species
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      getMutatedPos_bitset(mutatedPos, numMutablePosParent, r, 
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
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	tmpParam.popSize = 1;

	// FIXME00: This is ugly!!
	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN_bitset64(mutatedPos,
					      restrictTable,
					      popParams[nextMutant].death,
					      typeCBN,
					      newGenotype,
					      birthRate,
					      s,
					      numDrivers,
					      typeFitness);
	  tmpParam.birth = 1 - tmpParam.death;
	} else {
	  tmpParam.birth = fitness_CBN_bitset64(mutatedPos,
						  restrictTable,
						  popParams[nextMutant].birth,
						  typeCBN,
						  newGenotype,
						  birthRate,
						  s,
						  numDrivers,
						  typeFitness);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n Genotype = " << Genotypes[sp];
	  std::cout << "\n sp_id = " << Genotypes[sp].to_ulong();
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent Genotypes = " << Genotypes[nextMutant];
	}
#endif
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	tmpParam.mutation = mu * (numMutablePosParent - 1);

	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(tmpParam.death, tmpParam.birth, tmpParam.mutation);
	}
#ifdef DEBUGW
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors
#endif
	tmpParam.Flag = true;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (Genotype = " << Genotypes[sp] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")"
		    << "\n from species "  <<   nextMutant
		    << " (Genotypes = " << Genotypes[nextMutant] 
		    << "; sp_id = " << Genotypes[sp].to_ulong() << ")";
	}

	// FIXME00: the if can be removed??
	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  popParams[sp].death,
		  popParams[sp].birth);
	  if(verbosity >= 2) {
	    std::cout << "\n New popSize = " << popParams[sp].popSize << "\n";
	  }
	} else {
	  throw std::range_error("\n popSize == 0 but existing? \n");
	  if(verbosity >= 2) {
	    std::cout << "\n Mutation to an extinct species\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
#endif
	popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      popParams[nextMutant].Flag = true;
      // pop of receiving mutant flagged above 

    } else { //       *********** We are sampling **********
      if(verbosity >= 2) {
	std::cout <<"\n We are SAMPLING";   
	if(tSample < finalTime) {
	  std::cout << " at time " << tSample << "\n";
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }

      sample_all_pop_K(currentTime, sp_to_remove, 
		       popParams, Genotypes, tSample);

      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove[0])
	remove_zero_sp_v4(sp_to_remove, Genotypes, popParams);

      numSpecies = popParams.size();
      
      totPopSize_and_fill_out_crude(outNS_i, totPopSize, genot_out, //sp_id_out,
				    popSizes_out, index_out,
				    time_out, Genotypes, popParams, 
				    currentTime);

      //find_unique_genotypes2(uniqueGenotypes, totalNumSpecies, genot_out);
      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || (tSample >= finalTime)) 
	simulsDone = true;
	
      forceSample = false;
#ifdef DEBUGV
      std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }
  // FIXME: do I want to move this right after out_crude
  // and do it incrementally? I'd have also a counter of total unique species


  // FIXME: all this is ugly and could be a single function
  std::vector<unsigned long> genot_out_ulong(genot_out.size());
  genot_out_to_ulong(genot_out_ulong, genot_out);
  find_unique_genotypes(uniqueGenotypes, genot_out_ulong);
  std::vector<unsigned long> uniqueGenotypes_vector(uniqueGenotypes.size());
  uniqueGenotypes_to_vector(uniqueGenotypes_vector, uniqueGenotypes);


  // The out.ns in R code; holder of popSizes over time
  // The first row is time, then the genotypes (in column major)

  NumericMatrix outNS(uniqueGenotypes.size() + 1, outNS_i + 1);
  reshape_to_outNS(outNS, uniqueGenotypes_vector, genot_out_ulong, 
		   popSizes_out, 
		   index_out, time_out);
  
#ifdef DEBUGV
  std::cout << "\n Here X004 \n";
#endif  


  IntegerMatrix returnGenotypes(uniqueGenotypes_vector.size(), numGenes);

#ifdef DEBUGV
  std::cout << "\n Here X04 \n";
#endif  

  create_returnGenotypes(returnGenotypes, numGenes, uniqueGenotypes_vector);
  
#ifdef DEBUGV
  std::cout << "\n Here X6 \n";
#endif  


  return List::create(Named("NumSpecies") = uniqueGenotypes.size(), 
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = outNS,
		      Named("Genotypes") = returnGenotypes,
		      Named("FinalTime") = currentTime,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  END_RCPP
}





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

























