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

#include <cfloat> 
#include "bnb_common.h"
#include <Rcpp.h>



// Can be useful for debugging
void print_mapTimes(std::multimap<double, int>& mapTimes) {
  Rcpp::Rcout << "\n Printing mapTimes\n";
  for(auto elem : mapTimes) {
    Rcpp::Rcout << elem.first << "\t " << elem.second << "\n";
  }
}

// Can be useful for debugging
void print_initMutant(const std::vector < std::vector<int> >& initMutant) {
  Rcpp::Rcout <<"\n This is initMutant\n";
  for(size_t i = 0; i != initMutant.size(); ++i) {
    Rcpp::Rcout << "Init Mutant " << i
		<< ". Number of mutated genes: "
		<< initMutant[i].size()
		<<". Mutated genes: ";
    for(auto const &g : initMutant[i]) {
      Rcpp::Rcout << g << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "Finished printing initMutant \n";
}

// Can be useful for debugging
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
	      <<"\n numMutablePos = " << spP.numMutablePos
	      <<"\n";
}

// Can be useful for debugging mapFvarsValues
void print_map_string_double(const std::map<std::string, double>& efv) {
  Rcpp::Rcout << "\n Printing evalFVars_struct\n";
  for(auto elem : efv) {
    Rcpp::Rcout << elem.first << "\t " << elem.second << "\n";
  }
}

// Can be useful for debugging
void print_Genotype(const Genotype& ge) {
  Rcpp::Rcout << "\n Printing Genotype";
  Rcpp::Rcout << "\n\t\t order effects genes:";
  for(auto const &oo : ge.orderEff) Rcpp::Rcout << " " << oo;
  Rcpp::Rcout << "\n\t\t epistasis and restriction effects genes:";
  for(auto const &oo : ge.epistRtEff) Rcpp::Rcout << " " << oo;
  Rcpp::Rcout << "\n\t\t non interaction genes :";
  for(auto const &oo : ge.rest) Rcpp::Rcout << " " << oo;
  Rcpp::Rcout << "\n\t\t fitness landscape genes :";
  for(auto const &oo : ge.flGenes) Rcpp::Rcout << " " << oo;
  Rcpp::Rcout << std::endl;
}


// Compute pM of Mather's algorithm
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


// Find out time to next mutation
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

  // FIXMEMaybe: Should we short-circuit this when tSample = 0? Then, pM is 1,
  // and we could return ti = tSample + epsilon;
  
  using namespace Rcpp ;

  double r1;
  double ti;
  double pM;

  // FIXMEMaybe: should never happen
  if(spP.popSize <= 0.0) {
    throw std::range_error("ti: popSize <= 0. spP.popSize = "
			   + std::to_string(spP.popSize));

    // FIXMEclean: when we pass checks in BioC, Windows
    // #ifdef _WIN32
    //     throw std::range_error("ti: popSize <= 0. spP.popSize = "
    // 			   + SSTR(spP.popSize));
    // #endif
    // #ifndef _WIN32
    //         throw std::range_error("ti: popSize <= 0. spP.popSize = "
    // 			   + std::to_string(spP.popSize));
    // #endif
  }

  double invpop = 1/spP.popSize;
  double r;


  const double epsilon = 10.0;

  // W < 0 is a signal that mutation is zero, and thus ti is Inf
  if(spP.mutation == 0) { // FIXME: isn't this dead code?
    ti = tSample + 2.0 * epsilon;
    // yes, this is silly but to differentiate from
    // r < pM without further info
    // and to lead to finite value in loop for min.
    //ti = std::numeric_limits<double>::infinity();
  } else {

    RNGScope scope;
    r1 = ::Rf_runif(0.0, 1.0);
    // This was in the original Mather code, but I doubt
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

      // Is it really necessary to use log(-a) - log(-b) or could
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
      if(ti < 0.0) {
	double eq12 = pow( (spP.R - spP.W + 2.0 * spP.death) /
			   (spP.R + spP.W - 2.0 * spP.birth) , spP.popSize);

	Rcpp::Rcout << "\n ERROR: ti: eq.11 < 0 \n";
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

	DP_ti_notfinite;
	throw std::range_error("ti: ti not finite");
      }
      if((ti == 0.0) || (ti <= DBL_MIN)) {
	DEBUG_ti;
	++ti_dbl_min;
	ti = DBL_MIN;
	// Beware of this!!  
	// Abort because o.w. we can repeat it many, manu times
	Rcpp::Rcout << "         ti set to DBL_MIN: spP.popSize = " << spP.popSize << "\n";
	// It seems poSize over 1e8, and even 3.5 e7 can trigger this exception (depending
	// on mutation rate, of course)
	throw rerunExcept("ti set to DBL_MIN");
      }
      if(ti < (2*DBL_MIN)) ++ti_e3; // Counting how often this happens.
      // Can be smaller than the ti_dbl_min count
      ti += currentTime;
      // But we can still have issues here if the difference is too small
      if( (ti <= currentTime) ) {
	// Rcpp::Rcout << "\n (ti <= currentTime): expect problems\n";
	throw rerunExcept("ti <= currentTime");
      }
    }
  }
  return ti;
}

// Algorithm 2 of Mather
double Algo2_st(const spParamsP& spP,
		const double& ti,
		const int& mutationPropGrowth) { // need mutPropGrowth to
						 // know if we should
						 // throw

  // beware the use of t: now as it used to be, as we pass the value
  // and take the diff in here: t is the difference

  using namespace Rcpp ;
  double t = ti - spP.timeLastUpdate;

  if (spP.popSize == 0.0) {
    Rcpp::Rcout << "\n Entered Algo2 with pop size = 0\n";
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
    throw std::range_error("Algo 2: pb > 1 ");
  }

  if( pb < 0.0 ) {
    throw std::range_error("Algo 2: pb < 0");
  }
  //}


  if( pe == pm ) {
    // Should never happen. Exact identity??
    Rcpp::Rcout << "\n WARNING: Algo 2: pe == pm \n" ;
    return 0.0;
  }


  RNGScope scope;
  m = ::Rf_rbinom(spP.popSize, 1.0 - (pe/pm));
  // we can get issues with rbinom and odd numbers > 1e15
  // see "example-binom-problems.cpp"
  if(m <= 0.5) { // they are integers, so 0 or 1.
    if(m != 0.0)  Rcpp::Rcout << "\n WARNING: Algo 2: 0.0 < m < 0.5" <<std::endl;
    retval = 0.0;
  } else {
    rnb = ::Rf_rnbinom(m, 1.0 - pb); // this is the correct one
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

  double pm = pM_f_st(t, spP);
  double pe = pE_f_st(pm, spP);
  double pb = pB_f_st(pe, spP);

  double m; // the holder for the binomial
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
  // we can get issues with rbinom and odd numbers > 1e15
  // see "example-binom-problems.cpp"
  rnb = ::Rf_rnbinom(m + 2.0, 1.0 - pb);

  retval = m + 1 + rnb;

  if( !std::isfinite(retval) )  {
    DP2(rnb); DP2(m); DP2(pe); DP2(pm);
    print_spP(spP);
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


// initalize to absurd values an element of spParams
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



// For McF error 
// Get a -99 where there should be no error because of model

double returnMFE_new(double& en1,
		     const TypeModel typeModel) {
  if((typeModel == TypeModel::mcfarlandlog) ||
     (typeModel == TypeModel::mcfarlandlog_d) )
    return en1;
  else
    return -99;
}



void computeMcFarlandError_new(double& em1,
			       double& em1sc, // scaled
			       double& totPopSize_previous,
			       double& DA_previous,
			       const TypeModel typeModel,
			       const double& totPopSize,
			       const double& K){
  // Simple logic:
  // compute difference between successive death
  // rates, and also scale. Period.

  if( (typeModel == TypeModel::mcfarlandlog) ||
      (typeModel == TypeModel::mcfarlandlog_d)) {
    double etmp, etmpsc, DC;
    etmp = 0.0;
    etmpsc = 0.0;
    DC = -999999999; // o.w., we get a warning for possible uninitialized usage
    if(typeModel == TypeModel::mcfarlandlog) {
      DC = log1p(totPopSize/K);
    } else if (typeModel == TypeModel::mcfarlandlog_d) {
      DC = std::max(1.0, log1p(totPopSize/K));
    }
    if( std::abs(totPopSize - totPopSize_previous) < 1 ) {
      etmp = 0.0;
    } else {
      etmp = std::abs(DC - DA_previous);
      etmpsc = etmp/DA_previous;
    }
    if(etmp > em1) em1 = etmp;
    if(etmpsc > em1sc) em1sc = etmpsc;
    DA_previous = DC;
    totPopSize_previous = totPopSize;
  }
}

void updateRatesMcFarlandLog(std::vector<spParamsP>& popParams,
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


// Yes, identical to previous one, but with death rate a minimal value of 1
void updateRatesMcFarlandLog_D(std::vector<spParamsP>& popParams,
			       double& adjust_fitness_MF,
			       const double& K,
			       const double& totPopSize){

  adjust_fitness_MF = std::max(1.0,log1p(totPopSize/K));
  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}



void updateRatesFDFMcFarlandLog(std::vector<spParamsP>& popParams,
				const std::vector<Genotype>& Genotypes,
				const fitnessEffectsAll& fitnessEffects,
				double& adjust_fitness_MF,
				const double& K,
				const double& totPopSize,
				const double& currentTime) {

  const std::vector<spParamsP>& lastPopParams = popParams;

  adjust_fitness_MF = log1p(totPopSize/K);
  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    popParams[i].birth =
      prodFitness(evalGenotypeFitness(Genotypes[i],
				      fitnessEffects, Genotypes, lastPopParams, currentTime));
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}

// Yes, identical to previous one, but with death rate a minimal value of 1
void updateRatesFDFMcFarlandLog_D(std::vector<spParamsP>& popParams,
				  const std::vector<Genotype>& Genotypes,
				  const fitnessEffectsAll& fitnessEffects,
				  double& adjust_fitness_MF,
				  const double& K,
				  const double& totPopSize,
				  const double& currentTime) {

  const std::vector<spParamsP>& lastPopParams = popParams;

  // Min death rate is 1.0
  adjust_fitness_MF = std::max(1.0,log1p(totPopSize/K));
  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death = adjust_fitness_MF;
    popParams[i].birth =
      prodFitness(evalGenotypeFitness(Genotypes[i],
				      fitnessEffects, Genotypes, lastPopParams, currentTime));
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }

}


void updateRatesFDFExp(std::vector<spParamsP>& popParams,
		       const std::vector<Genotype>& Genotypes,
		       const fitnessEffectsAll& fitnessEffects,
		       const double& currentTime) {
  
  const std::vector<spParamsP>& lastPopParams = popParams;

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].birth =
      prodFitness(evalGenotypeFitness(Genotypes[i],
				      fitnessEffects, Genotypes, lastPopParams, currentTime));
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }
}

void updateRatesFDFBozic(std::vector<spParamsP>& popParams,
			 const std::vector<Genotype>& Genotypes,
			 const fitnessEffectsAll& fitnessEffects,
			 const double& currentTime) {
  
  const std::vector<spParamsP>& lastPopParams = popParams;

  for(size_t i = 0; i < popParams.size(); ++i) {
    popParams[i].death =
      prodDeathFitness(evalGenotypeFitness(Genotypes[i],
					   fitnessEffects, Genotypes, lastPopParams, currentTime));
    W_f_st(popParams[i]);
    R_f_st(popParams[i]);
  }

}


// Update birth or death rates as appropriate
// Exp and Bozic are not updated, unless FDF
void updateBirthDeathRates(std::vector<spParamsP>& popParams,
			   const std::vector<Genotype>& Genotypes,
			   const fitnessEffectsAll& fitnessEffects,
			   double& adjust_fitness_MF,
			   const double& K,
			   const double& totPopSize,
			   const double& currentTime,
			   const TypeModel typeModel) {

  if(!fitnessEffects.frequencyDependentFitness) {
    if (typeModel == TypeModel::mcfarlandlog)
      updateRatesMcFarlandLog(popParams, adjust_fitness_MF, K, totPopSize);
    if(typeModel == TypeModel::mcfarlandlog_d)
      updateRatesMcFarlandLog_D(popParams, adjust_fitness_MF, K, totPopSize);
  } else { // FDF
    if( typeModel == TypeModel::mcfarlandlog)  
      updateRatesFDFMcFarlandLog(popParams, Genotypes, fitnessEffects,
				 adjust_fitness_MF, K, totPopSize, currentTime);
    
    else if( typeModel == TypeModel::mcfarlandlog_d) 
      updateRatesFDFMcFarlandLog_D(popParams, Genotypes, fitnessEffects,
				   adjust_fitness_MF, K, totPopSize, currentTime);
	  
    else if(typeModel == TypeModel::exp) 
      updateRatesFDFExp(popParams, Genotypes, fitnessEffects, currentTime);

    else if(typeModel == TypeModel::bozic1)
      updateRatesFDFBozic(popParams, Genotypes, fitnessEffects, currentTime);
    else throw std::invalid_argument("this ain't a valid typeModel");
  } 
}



// Update the map times <-> indices and the popParams.pv
void mapTimes_updateP(std::multimap<double, int>& mapTimes,
		      std::vector<spParamsP>& popParams,
		      const int index,
		      const double time) {
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
    perSampleStats(i, 1) = sampleLargestPopSize[i]; // Never used in R FIXMEMaybe: remove!!
    perSampleStats(i, 2) = sampleLargestPopProp[i]; // Never used in R
    perSampleStats(i, 3) = static_cast<double>(sampleMaxNDr[i]);
    perSampleStats(i, 4) = static_cast<double>(sampleNDrLargestPop[i]);
  }
}

// Do not use this routinely. Too expensive and not needed.
void detect_ti_duplicates(const std::multimap<double, int>& m,
			  const double ti,
			  const int species) {

  double maxti = m.rbegin()->first;
  if((ti < maxti) && (m.count(ti) > 1)) {
    Rcpp::Rcout << "\n *** duplicated ti for species " << species << "\n";

    std::multimap<double, int>::const_iterator it = m.lower_bound(ti);
    std::multimap<double, int>::const_iterator it2 = m.upper_bound(ti);

    while(it != it2) {
      Rcpp::Rcout << "\tgenotype: " << (it->second) << "; time: " <<
	(it->first) << "\n";
      ++it;
    }
    Rcpp::Rcout << "\n\n\n";
  }
}

// Give message output of the state of the simulation
void message1(const int verbosity, const std::string message,
	      const int iteration, const double currentTime,
	      const unsigned int numSpecies,
	      const double totalPopulationSize,
	      const double timeNextPopSample,
	      const double minNextMutationTime) {
  if(verbosity >= 2)
    Rcpp::Rcout << "\n\n Verbose message at " << message
		<< ". Iteration = " << iteration
		<< ". currentTime =" << currentTime
		<< ". numSpecies = " << numSpecies
		<< ". totalPopulationSize " << totalPopulationSize
		<< ". timeNextPopSample " << timeNextPopSample
		<< ". minNextMutationTime " << minNextMutationTime
		<< "\n";
}

// another verbose message about the simulation
void messageNewSpecies(const int verbosity,
		       const int iteration, 
		       const unsigned int numSpecies,
		       const int nextMutant) {
  if(verbosity >= 2)
    Rcpp::Rcout <<"\n     Creating new species   " << (numSpecies - 1)
		<< " from species "  <<   nextMutant
		<<" at iteration " << iteration << "\n";
}

// more messages about the state of the simulation
void vvmessageNewSpecies(const int verbosity,
			 const unsigned int sp,
			 const Genotype& newGenotype,
			 const Genotype& parentGenotype,
			 const spParamsP& tmpParam,
			 const spParamsP& parentParam) {
  if(verbosity >= 2) {
    Rcpp::Rcout << " \n\n\n Looking at NEW species " << sp << " at creation";
    Rcpp::Rcout << "\n New Genotype :";
    print_Genotype(newGenotype);
    Rcpp::Rcout << "\n Parent Genotype :";
    print_Genotype(parentGenotype);
    Rcpp::Rcout << "\n birth of sp = " << tmpParam.birth;
    Rcpp::Rcout << "\n death of sp = " << tmpParam.death;
    Rcpp::Rcout << "\n parent birth = " << parentParam.birth;
    Rcpp::Rcout << "\n parent death = " << parentParam.death;
    Rcpp::Rcout << "\n\n popParams parent: \n";
    print_spP(parentParam);
    Rcpp::Rcout << "\n\npopParams child: \n";
    print_spP(tmpParam);
  }    
}

void messageSampling(const int verbosity,
		     const double tSample,
		     const double finalTime,
		     std::vector<spParamsP>& popParams) {
  if(verbosity >= 2) {
    Rcpp::Rcout <<"\n We are SAMPLING";
    if(tSample < finalTime) {
      Rcpp::Rcout << " at time " << tSample << "\n";
    } else
      Rcpp::Rcout <<". We reached finalTime " << finalTime << "\n";
    Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size() << "\n";
  }
}


void messagePostSampling(const int verbosity,
			 std::vector<spParamsP>& popParams,
			 const double totPopSize) {
  if(verbosity >= 2) 
    Rcpp::Rcout << "\n popParams.size() before sampling " << popParams.size()
		<< "\n totPopSize after sampling " << totPopSize << "\n";
}


// set dummyMutationRate
double setDummyMutationRate(std::vector<double> mu, const int verbosity) {
  double mymindummy = 1.0e-11; //1e-10
  double targetmindummy = 1.0e-10; //1e-9
  double minmu = *std::min_element(mu.begin(), mu.end());
  // Very small, but no less than mymindummy, for numerical issues.
  // We can probably go down to 1e-13. 1e-16 is not good as we get lots
  // of pE.f not finite. 1e-15 is probably too close, and even if no pE.f
  // we can get strange behaviors.
  // Ensures between mymindummy and targetmindummy
  double dummyMutationRate = std::max(std::min(minmu/1.0e4, targetmindummy),
				      mymindummy);

  // This should very rarely happen:
  if(minmu <= dummyMutationRate) { 
    double newdd = minmu/10.0;
    Rcpp::Rcout << "WARNING: the smallest mutation rate is "
		<< "<= " << dummyMutationRate << ". That is a really small value "
		<< "(per-base mutation rate in the human genome is"
		<< " ~ 1e-11 to 1e-9). "
		<< "Setting dummyMutationRate to your min/10 = "
		<< newdd
		<< ". There can be numerical problems later.\n";
    dummyMutationRate = newdd;
  }

  if(verbosity >= 2) {
    Rcpp::Rcout << "\n dummyMutationRate set at " << dummyMutationRate
		<< ".  That is the smallest possible mutation rate and the one"
		<< " for the null event.";
  }
  return dummyMutationRate;
}



// Initialize the population from initMutants and initSize
//     Set birth/death rates of populations, store if required, etc.
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
	     ) {

  numSpecies = initSize.size();
  totPopSize = std::accumulate(initSize.begin(), initSize.end(), 0.0);
  
 
  // Placeholders for genotype and params created at    initialization
  Genotype newGenotype = wtGenotype();
  spParamsP tmpParam;
  init_tmpP(tmpParam);

  // Create temporal copy: we were passed a const. If there are no initMutants,
  // we create the initMutant with just the wildtype
  
  std::vector < std::vector<int> > tmpInitMutant;
  if(initMutant.size() == 0) {
    // Empty genotype
    tmpInitMutant = {{}};
  } else {
    tmpInitMutant = initMutant;
  }


  int numGenesInitMut = -99;
  const int numLoci = fitnessEffects.allGenes.size();

  // Loop twice: create genotypes and fill up population sizes in
  // popParams. Then compute fitness (which might be affected by identity
  // and population size of other genotypes), W, R, mutation, etc

  // We do not remove non-viable genotypes. 
  
  //  Fill up Genotypes and popParams
  for(size_t m = 0; m != tmpInitMutant.size(); ++m ) {
    init_tmpP(tmpParam);
    newGenotype = wtGenotype();
    std::vector<int> this_initMutant = tmpInitMutant[m];
    // Create the new genotype by adding mutations into a wtGenotype
    newGenotype = createNewGenotype(wtGenotype(),
				    this_initMutant, 
				    fitnessEffects,
				    ran_gen,
				    false);
    
    numGenesInitMut = newGenotype.orderEff.size() +
      newGenotype.epistRtEff.size() + newGenotype.rest.size() +
      newGenotype.flGenes.size();

    if( (!(numGenesInitMut == 0)) != (!(newGenotype == wtGenotype()))  )
      throw std::logic_error("InitMutant: Either a WT genotype without 0 mutations or a non-WT with 0 mutations");

    
    tmpParam.numMutablePos = numLoci - numGenesInitMut;
    // Next unreachable since caught in R.
    // But just in case, since it would lead to seg fault.
    if(tmpParam.numMutablePos < 0)
      throw std::invalid_argument("initMutant's genotype has more genes than are possible.");

    
    tmpParam.popSize = initSize[m];
    popParams.push_back(tmpParam);
    Genotypes.push_back(newGenotype);
    // birth, death, W, R, absfitness: updated when fitness, below
    // mutation: below, when calling mutationFromScratch
    // pv: when calling mapTimes_updateP
    // timeLastUpdate updated below, right after pv

  }

  tmpInitMutant.clear();
   
  // Assign fitness, mutation, W, R. In that order. Then pv and add to POM.
  // The last two are specific of BNB. Beware if using a different algorithm
  for(size_t m = 0; m != popParams.size(); ++m ) {
    if(typeModel == TypeModel::bozic1) {
      popParams[m].death =
	prodDeathFitness(evalGenotypeFitness(Genotypes[m],
					     fitnessEffects, Genotypes, popParams,
					     currentTime));
      popParams[m].birth = 1.0;
      if(!fitnessEffects.frequencyDependentFitness &&
	 (popParams[m].numMutablePos == numLoci ) &&
	 (popParams[m].death != 1.0)) throw std::logic_error("WT initMutant in non-FDF must have death rate 1 with this model");
      // Usual runs without initMutant can use this
      if( (popParams[m].death == 1.0) &&
	  (initMutant.size() != 0)) Rcpp::Rcout << "Init Mutant with death == 1.0\n";
      if(popParams[m].death >  99) Rcpp::warning("Init Mutant with death > 99");
    } else {
      
      popParams[m].birth =
	prodFitness(evalGenotypeFitness(Genotypes[m],
					fitnessEffects, Genotypes, popParams,
					currentTime));
      if (typeModel == TypeModel::exp) 
	popParams[m].death = death; // passed from R; set at 1
      if (typeModel == TypeModel::mcfarlandlog)
	popParams[m].death = log1p(totPopSize/K);
      if (typeModel == TypeModel::mcfarlandlog_d)
	popParams[m].death = std::max(1.0, log1p(totPopSize/K));
      
      if(!fitnessEffects.frequencyDependentFitness &&
	 (popParams[m].numMutablePos == numLoci ) &&
	 (popParams[m].birth != 1.0) ) throw std::logic_error("WT initMutant in non-FDF must have birth rate 1 with this model");
      if((popParams[m].birth == 1.0) &&
	 (initMutant.size() != 0) ) Rcpp::Rcout << "Init Mutant with birth == 1.0\n";
      if(popParams[m].birth == 0.0) Rcpp::warning("Init Mutant with birth == 0.0");
    }

    popParams[m].mutation = mutationFromScratch(mu, popParams[m], Genotypes[m],
						fitnessEffects, mutationPropGrowth,
						full2mutator, muEF,
						Genotypes, popParams,
						currentTime,
						dummyMutationRate);
    W_f_st(popParams[m]);
    R_f_st(popParams[m]);
  }

  // POM and storing output
  // Code repeated from end of  nr_totPopSize_and_fill_out_crude_P
  // (except the keepEvery condition)
  // Finding max_ndr from code at start of that function
  
  if( keepEvery > 0 ) {
    // We keep the first genotype(s) ONLY if we are storing more than one.
    lastStoredSample = currentTime;
    outNS_i++;
    int ndr_lp = 0;
    double l_pop_s = 0.0;
    int largest_clone = -99;

    int tmp_ndr = 0;
    int max_ndr = 0;

    time_out.push_back(currentTime);

    for(size_t i = 0; i < popParams.size(); ++i) {
      genot_out.push_back(Genotypes[i]);
      popSizes_out.push_back(popParams[i].popSize);
      index_out.push_back(outNS_i);

      tmp_ndr = getGenotypeDrivers(Genotypes[i], fitnessEffects.drv).size();
      if(tmp_ndr > max_ndr) max_ndr = tmp_ndr;
      
      if(popParams[i].popSize > l_pop_s) {
	l_pop_s = popParams[i].popSize;
	ndr_lp = getGenotypeDrivers(Genotypes[i], fitnessEffects.drv).size();
	largest_clone = i;
      }
    }
    
    sampleTotPopSize.push_back(totPopSize); //totPopSize computed above
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
  } else {
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

  if(verbosity > 2) {
    Rcpp::Rcout << "\n Population right after initialization\n";
    for(size_t i = 0; i < popParams.size(); ++i) {
      print_spP(popParams[i]);
    }
  }

}
