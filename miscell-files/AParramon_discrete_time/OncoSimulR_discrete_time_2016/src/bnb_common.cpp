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



#include "bnb_common.h"
#include "new_restrict.h" // for the TypeModel enum
#include <Rcpp.h>


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

double pM_f_st(const double& t, 
	       const spParamsP& spP){
  // For interpretation, recall, from suppl. mat. of their paper, p.2 that
  // p M (t)^n0 = G(0, t) is the probability that a mutation has not yet occurred.
  
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

double ti_nextTime_tmax_2_st(const spParamsP& spP,
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
    
#ifdef _WIN32     
    throw std::range_error("ti: popSize <= 0. spP.popSize = "
			   + SSTR(spP.popSize));
#endif
    
#ifndef _WIN32 
        throw std::range_error("ti: popSize <= 0. spP.popSize = "
			   + std::to_string(spP.popSize));
#endif
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
	// abort simulation and go to a new one?  
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

double Algo2_st(const spParamsP& spP,
		const double& ti,
		const int& mutationPropGrowth) { // need mutPropGrowth to
						 // know if we should
						 // throw

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


  if( (spP.mutation == 0.0) &&
      !(spP.birth <= 0 && mutationPropGrowth) ) {
    Rcpp::Rcout << "\n Entered Algo2 with mutation rate = 0\n";
    if( spP.numMutablePos != 0 ) 
      throw std::range_error("mutation = 0 with numMutable != 0?");
  }


  // double pm, pe, pb;
  double m; // the holder for the binomial

  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);

  // if(spP.numMutablePos == 0) {
  //   // Just do the math. In this case mutation rate is 0. Thus, pM (eq. 8
  //   // in paper) is 1 necessarily. And pE is birth/death rates (if you do
  //   // things sensibly, not multiplying numbers as in computing pE
  //   // below). And pB is also 1.

  //   // This will blow up below if death > birth, as then pe/pm > 1 and 1 -
  //   // pe/pm < 0. But I think this would just mean this is extinct?
  //   pm = 1;
  //   pe = spP.death/spP.birth;
  //   pb = 1;
  // } else {
  //   pm = pM_f_st(t, spP);
  //   pe = pE_f_st(pm, spP);
  //   pb = pB_f_st(pe, spP);
  // }

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
      Rcpp::Rcout << "\n ERROR: Algo 2, (1.0 - pe/pm) < 0.0 \n"
		  << " t = " << t << "; R = " << spP.R  
		  <<  "; W = " << spP.W << ";\n death = " << spP.death 
		  <<  "; growth = " << spP.birth << ";\n pm = " << pm 
		  << "; pe = " << pe << "; pb = " << pb << std::endl;
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
    rnb = ::Rf_rnbinom(m, 1.0 - pb); // this is the correct ONE
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

double Algo3_st(const spParamsP& spP, const double& t){
  
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




void precissionLoss(){
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

void init_tmpP(spParamsP& tmpParam) {
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



// this is the log of the ratio of death rates
// so the the difference of the successie death rates, if using
// the log version.
// double returnMFE(double& e1,
// 			const double& K,
// 			const std::string& typeFitness) {
//   if((typeFitness == "mcfarland0") || (typeFitness == "mcfarlandlog"))
//     return log(e1);
//   else if(typeFitness == "mcfarland")
//     return ((1.0/K) * e1);
//   else
//     return -99;
// }


// double returnMFE(double& e1,
// 			const double& K,
// 			const TypeModel typeModel) {
//   if((typeModel == TypeModel::mcfarland0) || (typeModel == TypeModel::mcfarlandlog))
//     return log(e1);
//   else if(typeModel == TypeModel::mcfarland)
//     return ((1.0/K) * e1);
//   else
//     return -99;
// }

double returnMFE(double& e1,
		 // const double& K,
		 const std::string& typeFitness) {
  if(typeFitness == "mcfarlandlog")
    return log(e1);
  else
    return -99;
}

double returnMFE(double& e1,
		 // const double& K,
			const TypeModel typeModel) {
  if(typeModel == TypeModel::mcfarlandlog)
    return log(e1);
  else
    return -99;
}



// FIXME But I'd probably want a percent error, compared to the death rate
// something like (log(1+N1/K) - log(1+N2/K))/(log(1+N1/K))


// void computeMcFarlandError(double& e1,
// 				  double& n_0,
// 				  double& n_1,
// 				  double& tps_0,
// 				  double& tps_1,
// 				  const std::string& typeFitness,
// 				  const double& totPopSize,
// 				  const double& K){
//   //				  const double& initSize) {
//   // static double tps_0 = initSize;
//   // static double tps_1 = 0.0;

//   if( (typeFitness == "mcfarland0") ||
//       (typeFitness == "mcfarland") || 
//       (typeFitness == "mcfarlandlog") ) {
    
//     double etmp;
//     tps_1 = totPopSize;
//     if(typeFitness == "mcfarland")
//       etmp = std::abs( tps_1 - (tps_0 + 1) );
//     else {
//       if( (tps_0 + 1.0) > tps_1 ) 
// 	etmp = (K + tps_0 + 1.0)/(K + tps_1);
//       else
// 	etmp = (K + tps_1)/(K + tps_0 + 1);
//     }
//     if(etmp > e1) {
//       e1 = etmp;
//       n_0 = tps_0;
//       n_1 = tps_1;
//     }
//     tps_0 = tps_1;
//   }
// }

void computeMcFarlandError(double& e1,
			   double& n_0,
			   double& n_1,
			   double& tps_0,
			   double& tps_1,
			   const std::string& typeFitness,
			   const double& totPopSize,
			   const double& K){

  if(typeFitness == "mcfarlandlog")  {
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


// void computeMcFarlandError(double& e1,
// 				  double& n_0,
// 				  double& n_1,
// 				  double& tps_0,
// 				  double& tps_1,
// 				  const TypeModel typeModel,
// 				  const double& totPopSize,
// 				  const double& K){
//   //				  const double& initSize) {
//   // static double tps_0 = initSize;
//   // static double tps_1 = 0.0;

//   if( (typeModel == TypeModel::mcfarland0) ||
//       (typeModel == TypeModel::mcfarland) || 
//       (typeModel == TypeModel::mcfarlandlog) ) {
//     double etmp;
//     tps_1 = totPopSize;
//     if(typeModel == TypeModel::mcfarland)
//       etmp = std::abs( tps_1 - (tps_0 + 1) );
//     else {
//       if( (tps_0 + 1.0) > tps_1 ) 
// 	etmp = (K + tps_0 + 1.0)/(K + tps_1);
//       else
// 	etmp = (K + tps_1)/(K + tps_0 + 1);
//     }
//     if(etmp > e1) {
//       e1 = etmp;
//       n_0 = tps_0;
//       n_1 = tps_1;
//     }
//     tps_0 = tps_1;
//   }
// }

void computeMcFarlandError(double& e1,
			   double& n_0,
			   double& n_1,
			   double& tps_0,
			   double& tps_1,
			   const TypeModel typeModel,
			   const double& totPopSize,
			   const double& K){

  if( typeModel == TypeModel::mcfarlandlog ) {
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


// void updateRatesMcFarland(std::vector<spParamsP>& popParams,
// 				 double& adjust_fitness_MF,
// 				 const double& K,
// 				 const double& totPopSize){

//   adjust_fitness_MF = totPopSize/K;

//   for(size_t i = 0; i < popParams.size(); ++i) {
//     popParams[i].death = adjust_fitness_MF;
//     W_f_st(popParams[i]);
//     R_f_st(popParams[i]);
//   }
// }


void updateRatesMcFarlandLog(std::vector<spParamsP>& popParams,
			     double& adjust_fitness_MF,
			     const double& K,
			     const double& totPopSize){

  // from original log(1 + totPopSize/K)
  adjust_fitness_MF = log1p(totPopSize/K);
  
  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}


// // McFarland0 uses: - penalty as log(1 + N/K), and puts
// // that in the birth rate.
// void updateRatesMcFarland0(std::vector<spParamsP>& popParams,
// 				  double& adjust_fitness_MF,
// 				  const double& K,
// 				  const double& totPopSize,
// 				  const int& mutationPropGrowth,
// 				  const double& mu){
  
//   adjust_fitness_MF = 1.0 / log1p(totPopSize/K);

//   for(size_t i = 0; i < popParams.size(); ++i) {
//     popParams[i].birth = adjust_fitness_MF * popParams[i].absfitness;
//     if(mutationPropGrowth) {
//       popParams[i].mutation = mu * popParams[i].birth * 
// 	popParams[i].numMutablePos;
//     } else if(popParams[i].birth / popParams[i].mutation < 20) {
//       Rcpp::Rcout << "\n WARNING: birth/mutation < 20";
//       Rcpp::Rcout << "\n Birth = " << popParams[i].birth 
// 		<< ";  mutation = " << popParams[i].mutation << "\n";
//     }
//     W_f_st(popParams[i]);
//     R_f_st(popParams[i]);
//   }
// }

// void updateRatesBeeren(std::vector<spParamsP>& popParams,
// 			      double& adjust_fitness_B,
// 			      const double& initSize,
// 			      const double& currentTime,
// 			      const double& alpha,
// 			      const double& totPopSize,
// 			      const int& mutationPropGrowth,
// 			      const double& mu){

//   double average_fitness = 0.0; // average_fitness in Zhu
//   double weighted_sum_fitness = 0.0;
//   double N_tilde;
  
//   for(size_t i = 0; i < popParams.size(); ++i) {
//     weighted_sum_fitness += (popParams[i].absfitness * popParams[i].popSize);
//   }
  
//   average_fitness = (1.0/totPopSize) * weighted_sum_fitness;
//   N_tilde =  initSize * exp(alpha * average_fitness * currentTime);
//   adjust_fitness_B = N_tilde/weighted_sum_fitness; 

//   if(adjust_fitness_B < 0) {
//     throw std::range_error("adjust_fitness_B < 0");
//   }
  
//   for(size_t i = 0; i < popParams.size(); ++i) {
//     popParams[i].birth = adjust_fitness_B * popParams[i].absfitness;
//     if(mutationPropGrowth) {
//       popParams[i].mutation = mu * popParams[i].birth * 
// 	popParams[i].numMutablePos;
//     } else if(popParams[i].birth / popParams[i].mutation < 20) {
//       Rcpp::Rcout << "\n WARNING: birth/mutation < 20";
//       Rcpp::Rcout << "\n Birth = " << popParams[i].birth 
// 		<< ";  mutation = " << popParams[i].mutation << "\n";
//     }
//     W_f_st(popParams[i]);
//     R_f_st(popParams[i]);
//   }
// }





void mapTimes_updateP(std::multimap<double, int>& mapTimes,
		      std::vector<spParamsP>& popParams,
		      const int index,
		      const double time) {
  // Update the map times <-> indices
  // Recall this is the map of nextMutationTime and index of species
  // First, remove previous entry, then insert.
  // But if we just created the species, nothing to remove from the map.
  if(popParams[index].timeLastUpdate > -1)
    mapTimes.erase(popParams[index].pv);
  popParams[index].pv = mapTimes.insert(std::make_pair(time, index));
}


void getMinNextMutationTime4(int& nextMutant, double& minNextMutationTime,
			     const std::multimap<double, int>& mapTimes) {
  // we want minNextMutationTime and nextMutant
  nextMutant = mapTimes.begin()->second;
  minNextMutationTime = mapTimes.begin()->first;
}


void fill_SStats(Rcpp::NumericMatrix& perSampleStats,
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
