struct spParamsF {
  bool Flag;
  int pos_in_NS;
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  double nextMutationTime;
  double timeLastUpdate;
};


void getMinNextMutationTime(int& nextMutant, double& minNextMutationTime,
			     const std::vector<spParamsF>& popParams) {
  
  // we want minNextMutationTime and nextMutant
  // turned into a function for profiling
  minNextMutationTime = std::numeric_limits<double>::infinity();
  nextMutant = -99; 
  
  for(int i = 0; i < popParams.size(); i++) {
    if(popParams[i].popSize > 0.0) {
      if(popParams[i].nextMutationTime < minNextMutationTime) {
	nextMutant = i;
	minNextMutationTime = popParams[i].nextMutationTime;
      }     
    }
  }
}

void remove_zero_sp_v2(const std::vector<int>& sp_to_remove,
		       std::vector<unsigned long>& sp_id,
		       std::vector<std::vector<myT> >& Genotypes,
		       std::vector<spParamsF>& popParams) {

  std::vector<unsigned long>::iterator sp_begin = sp_id.begin();
  std::vector<std::vector<myT> >::iterator Genotypes_begin = Genotypes.begin();
  std::vector<spParamsF>::iterator popParams_begin = popParams.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    sp_id.erase(sp_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
    popParams.erase(popParams_begin + sp_to_remove[j]);
  }
}


void totPopSize_and_fill_outNS(double& totPopSize, int& outNS_i, 
			       arma::mat& outNS, 
			       const std::vector<spParamsF>& popParams, 
			       //const int& iter,
			       const std::vector<unsigned long>& sp_id,
			       const double& currentTime) {
  // Actually, fill out, but also compute totPopSize

  #ifdef DEBUGV
      std::cout << "\n Filling up outNS \n";
#endif
      
      outNS_i++;
      //outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(0, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      
      for(int i = 0; i < popParams.size(); ++i) {
	outNS(i + 1, outNS_i) = popParams[i].popSize;
	// I need one additional subscripting, unneeded if I keep all ever non-zero
	// as popParams[i].pos_in_NS == i
	// outNS(popParams[i].pos_in_NS + 1, outNS_i) = popParams[i].popSize;
	totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	std::cout << "\n       Species " << i 
		  // << ", pos_in_NS = " << popParams[i].pos_in_NS
		  // << ", equal pos and i? " << (popParams[i].pos_in_NS == i)
		  << ", sp_id = " << sp_id[i] 
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



// using std vector for genotypes
double fitness_CBN_std(const int& mutatedPos, 
		       Rcpp::IntegerMatrix restrictTable,
		       const double& fitnessParent, 
		       const std::string typeCBN,
		       const std::vector<myT>& Genotype,
		       const double& birthRate, 
		       const double& s, 
		       const int& numDrivers) {
    
  using namespace Rcpp ;

  // FIXME: check numDrivers == ncol(restrictTable)??
  // or obtain the numDrivers from there??
  // Yes, numDrivers is known before hand. Take if from there
  // in some upper level function that calls this one.


  // FIXME: add checks of no negative numbers in deps: thisRestrict?

  // In fact, I could simplify, by having the restrictTable only
  // for drivers that DO depend. Would make it faster and exclude
  // a check below. But less clear.

  
  // will later become an argument
  double fitnessNo = 0.0;

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

    if(outFitnessYes) {
      return fitness_linear_std(Genotype, birthRate, s, numDrivers);
    } else {
      return fitnessNo;
    }

  }
}


// using std vector for genotypes
// but two different typeFitness
double fitness_CBN(const int& mutatedPos, 
		       Rcpp::IntegerMatrix restrictTable,
		       const double& fitnessParent, 
		       const std::string typeCBN,
		       const std::vector<myT>& Genotype,
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

    if(outFitnessYes) {
      if(typeFitness == "bozic")
	return fitness_bozic(Genotype, s, numDrivers);
      else
	return fitness_linear_std(Genotype, birthRate, s, numDrivers);
    } else {
      if(typeFitness == "bozic")
	return 1.0; //death
      else
	return 0.0; //birth
    }
  }
}

// this uses Genotype is an std vector of a vector of ints
inline double fitness_linear_std(const std::vector<myT>& Genotype,
				 const double& birthRate, 
				 const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype

  int totalMut = std::accumulate(Genotype.begin(), 
				 Genotype.begin() + numDrivers,
				 0);
  return birthRate + s * static_cast<double>(totalMut);
}


// this uses Genotype is an std vector of a vector of ints
// This returns death!!! not birth
inline double fitness_bozic(const std::vector<myT>& Genotype,
			    const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype

  int totalMut = std::accumulate(Genotype.begin(), 
				 Genotype.begin() + numDrivers,
				 0);
  return 0.5 * pow( 1.0 - s, totalMut); 
  // return 1.0 - 0.5 * pow( 1.0 - s, totalMut); 

}

// Unclear which of the following two is faster.
// Limited benchmarking on Intel suggests longint1 might be just
// a tiny bit faster. (new_sp_longint2 is now in older-code.cpp)
void new_sp_longint1(int& sp, unsigned long& new_id,
		     const std::vector<unsigned long>& sp_id,
		     const int& numSpecies,
		     const int& nextMutant,
		     const int& mutatedPos){
  unsigned long myone = static_cast<unsigned long> (1);
  unsigned long add_to_id = (myone << mutatedPos);
  new_id = add_to_id + sp_id[nextMutant];
    
  sp = 0;

  for(sp = 0; sp < numSpecies; ++sp) {
    if( new_id == sp_id[sp] )
      break;
  }  
  
}

void getMutatedPos(int& mutatedPos, int& numMutablePosParent,
		   gsl_rng *r,
		   std::vector<int>& mutablePos,
		   const int& nextMutant,
		   const std::vector<std::vector<myT> >& Genotypes,
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
    if(!Genotypes[nextMutant][i]) { 
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

void remove_zero_sp_v5(const std::vector<int>& sp_to_remove,
		       std::vector<Genotype64>& Genotypes,
		       std::vector<spParamsK>& popParams,
		       std::vector<double>& nextMutTime) {

  std::vector<spParamsK>::iterator popParams_begin = popParams.begin();
  std::vector<Genotype64>::iterator Genotypes_begin = Genotypes.begin();
  std::vector<double>::iterator nextMutTime_begin = nextMutTime.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    popParams.erase(popParams_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
    nextMutTime.erase(nextMutTime_begin + sp_to_remove[j]);

  }
}

void resize_outNS(arma::mat& outNS, const int& outNS_i,
		  const int& numSpecies) {
  int type_resize = 0;

  if( outNS.n_rows <= (numSpecies + 1) ) type_resize += 1;
  if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

  if(type_resize == 0) return;
  else if(type_resize == 1) 
    outNS.resize(2 * numSpecies + 1, outNS.n_cols);
  else if (type_resize == 2)
    outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
  else if (type_resize == 3)
    outNS.resize(2 * numSpecies + 1, 2 * outNS_i + 2);
}



void getMinNextMutationTime2(int& nextMutant, double& minNextMutationTime,
			     std::vector<spParamsK>& popParams) {
  //nope: popParams can not go as const
  // we want minNextMutationTime and nextMutant
  // does not work if species with popSize == 0

  std::vector<spParamsK>::iterator pt_pos_min =
    std::min_element(popParams.begin(), popParams.end(), popParamsK_time_less);
  nextMutant = std::distance(popParams.begin(), pt_pos_min);
  minNextMutationTime = popParams[nextMutant].nextMutationTime;
}



double ti_nextTime_tmax(const double& R, const double& W, 
			const double& death, const double& growth, 
			const double& n, const double& currentTime,
			double& tSample) {
  // Following the logic of the code by Mather in
  // findNextMutationTime

  // We return the nextMutationTime or a value larger than the
  // max length of the period (tSample)

  // I also change names rr, r, to match those in Mather r1, r.
  using namespace Rcpp ;

  double eq11;
  double r1;
  double r;
  double ti;
  double pM;
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
      //eq. 11

      ti = (1.0/R) * (log( -1 * (r * (R - W + 2.0 * growth) - W - R + 2.0 * death )) -
       		      log( -1 * (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death )));

      // ti = (1.0/R) * log( (r * (R - W + 2.0 * growth) - W - R + 2.0 * death) /
      //               (r * (-R -W + 2.0 * growth) - W + R + 2.0 * death));
      // std::cout << "\n this is ti = " << ti << "\n";
      if(ti < 0.0) {
	
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);

	std::cout << "\n ERROR: ti: eq.11 < 0 \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death << "\n\n";
	throw std::range_error("ti: eq.11 < 0");
      } 
      if( !std::isfinite(ti) ) {
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);
	std::cout << "\n ERROR: ti: not finite \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death;
	throw std::range_error("ti: ti not finite");
      }
      if(ti == 0.0) {
	double eq12 = pow( (R - W + 2.0 * death) / (R + W - 2.0 * growth) , n);
	std::cout << "\n WARNING: ti == 0. Expect problems \n" 
		  << " especially if using forced sampling \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n r1 = " << r1;
	std::cout << "\n is r > 1? " << (r > 1.0) << "\n";
	std::cout << "\n is r < 0? " << (r < 0.0) << "\n";
	std::cout << "\n is eq12 < r? " << (eq12 < r) << "\n";
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death << "\n\n";
      }
      ti += currentTime;
    } 
  }
  return ti;
}



double ti_f(const double& R, const double& W, 
	    const double& death, const double& growth, 
	    const double& n) {
  using namespace Rcpp ;

  double eq11;
  double r;
  double rr;
  double ti;

  // W < 0 is a signal that mutation is zero, and thus ti is Inf
  if(W <= -90.0) {
    ti = std::numeric_limits<double>::infinity();
  } else {

    RNGScope scope;
    r = ::Rf_runif(0.0, 1.0);
    rr = exp((1.0 / n) * log(r));

    //    std::cout << "r = " << r << " rr = " << rr << std::endl;
    
    // eq.12
    if( ((R - W + 2.0 * death)/(R + W - 2.0 * growth)) < rr ) {
      //eq. 11
      ti = (1.0/R) * log( (rr * (R - W + 2.0 * growth) - W - R + 2.0 * death) /
                     (rr * (-R -W + 2.0 * growth) - W + R + 2.0 * death));
      if(ti < 0.0) {
	std::cout << "\n ERROR: ti: eq.11 < 0 \n";
	std::cout << "\n R = " << R;
	std::cout << "\n W = " << W;
	std::cout << "\n r = " << r;
	std::cout << "\n n = " << n;
	std::cout << "\n rr = " << rr;
	std::cout << "\n growth = " << growth;
	std::cout << "\n death = " << death;
	throw std::range_error("ti: eq.11 < 0");
      }
    } else {
      ti = std::numeric_limits<double>::infinity();
    }
  }
  return ti;
}

// Comparing with R code
//ss <- 1; set.seed(ss); d <- 0.3; g <- 0.300005; mu <- 0.005; n <- 300;
//exwrap(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g, n, 99);
//set.seed(ss); ti.original(R.f.c(d, g, mu, 0), W.f.c2(d, g, mu, 0), d, g,
//n)



void new_sp_gmp(int& sp, mpz_t& add_to_id, mpz_t& new_id,
		const std::vector<mpz_class>& sp_id,
		const int& numSpecies,
		const int& nextMutant,
		const int& mutatedPos){
  // We want sp and new_id
  // add_to_id could be created here, but I reuse space. Use static?

  // Turned into function for timing
  mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
  mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());

  sp = 0;

  for(sp = 0; sp < numSpecies; ++sp) {
    if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      break;
  }  
}


void new_sp_genotypes(int& sp, 
		      const int& numSpecies,
		      const int& numGenes,
		      const std::vector<myT>& newGenotype,
		      const std::vector<std::vector<myT> >& Genotypes) {
  // Find out if new species by looking at each position in the genotype.
  int k = 0;
  sp = 0;
  
  while ( (sp < numSpecies) && (k < numGenes) ) {
    if (newGenotype[k] == Genotypes[sp][k]) k++;
    else {
      sp++;
      k = 0;
    }
  }
}


void new_sp_longint2(int& sp, unsigned long& new_id,
		     const std::vector<unsigned long>& sp_id,
		     const int& numSpecies,
		     const int& nextMutant,
		     const int& mutatedPos){
  // using pow
  unsigned long add_to_id = static_cast<unsigned long>(pow(2, mutatedPos));
  new_id = add_to_id + sp_id[nextMutant];
  
  sp = 0;
  
  for(sp = 0; sp < numSpecies; ++sp) {
    if( new_id == sp_id[sp] )
      break;
  }  
  
}

unsigned long id_code(const std::vector<myT>& genotype) {
  // returns the same as what we used with sp_id
  // Could be made faster with accumulate? Not worth it.
  unsigned long myone = static_cast<unsigned long> (1);
  unsigned long a = 0;
  for(int i = 0; i < genotype.size(); ++i) {
    if(genotype[i]) {
      a += (myone << i);
    }
  }
  return a;
}



// void getMinNextMutationTime_pq(int& nextMutant, double& minNextMutationTime,
// 			       std::priority_queue<time_index> pq){
  
//   // we want minNextMutationTime and nextMutant
//   // turned into a function for profiling
//   minNextMutationTime = std::numeric_limits<double>::infinity();
//   nextMutant = -99; 

//   time_index t2; // FIXME: pass it?

//   t2 = pq.top();

  

//   for(int i = 0; i < popParams.size(); i++) {
//     if(popParams[i].popSize > 0.0) {
//       if(popParams[i].nextMutationTime < minNextMutationTime) {
// 	nextMutant = i;
// 	minNextMutationTime = popParams[i].nextMutationTime;
//       }     
//     }
//   }
// }



struct spParams {
  bool Flag;
  double birth;
  double popSize;
  double timeLastUpdate;
  double W;
  double R;
  double tis;
  double nextMutationTime;
double mutation; 
};

// without tis
struct spParamsC {
  bool Flag;
  double birth;
  double popSize;
  double timeLastUpdate;
  double W;
  double R;
  double nextMutationTime;
double mutation; 
};

// with bozic, I do 1- and then again 1 -. 
// But if I don't I'd need to change the logic.

// without tis
struct spParamsD {
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


void remove_zero_sp_v1(const std::vector<int>& sp_to_remove,
		       std::vector<mpz_class>& sp_id,
		       std::vector<std::vector<myT> >& Genotypes,
		       std::vector<spParamsF>& popParams) {

  std::vector<mpz_class>::iterator sp_begin = sp_id.begin();
  std::vector<std::vector<myT> >::iterator Genotypes_begin = Genotypes.begin();
  std::vector<spParamsF>::iterator popParams_begin = popParams.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    sp_id.erase(sp_begin + sp_to_remove[j]);
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
    popParams.erase(popParams_begin + sp_to_remove[j]);
  }
}

void remove_zero_sp_v3(const std::vector<int>& sp_to_remove,
		       std::vector<std::vector<myT> >& Genotypes,
		       std::vector<spParamsF>& popParams) {

  std::vector<std::vector<myT> >::iterator Genotypes_begin = Genotypes.begin();
  std::vector<spParamsF>::iterator popParams_begin = popParams.begin();

  for(int j = sp_to_remove[0]; j > 0 ; --j) {
    Genotypes.erase(Genotypes_begin + sp_to_remove[j]);
    popParams.erase(popParams_begin + sp_to_remove[j]);
  }
}


// struct time_index{
//   float time;
//   int index;
// };

// inline bool operator<(const time_index& a, const time_index& b){
//   return a.time > b.time;
// }



void sample_all_pop(int& nS_nonZero,
		    double& currentTime, 
		    std::vector<int>& sp_to_remove,
		    std::vector<spParamsF>& popParams,
		    const std::vector<unsigned long>& sp_id,
		    const double& tSample){

  currentTime = tSample;
  sp_to_remove[0] = 0;

  for(int i = 0; i < popParams.size(); i++) {
    if(popParams[i].popSize > 0.0) {

      STOPASSERT(popParams[i].Flag == false);
      STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
      STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
#ifdef DEBUGV
      std::cout << "\n\n     ********* 5.9 ******\n " 
		<< "     Species  = " << i 
		<< "\n      Genotype =     " // just so output matches with scroll-all-mode
		<< "\n      sp_id = " << sp_id[i]  
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
      if( (popParams[i].pos_in_NS == -99) && (popParams[i].popSize <=  0.0) ) {
	// this i has never been non-zero in any sampling time
	sp_to_remove[0]++;
	sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
	std::cout << "\n\n     Removing species i = " << i ;
	  //		  << " with sp_id = " << sp_id[i];
#endif
      } else {
	if(popParams[i].pos_in_NS == -99) {
	  popParams[i].pos_in_NS = nS_nonZero;	      
	  ++nS_nonZero;
	}
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
}


void create_return(Rcpp::NumericVector& popSizes,
		   Rcpp::IntegerMatrix& returnGenotypes,
		   arma::mat& outNS, 
		   const int& outNS_i,
		   const int& numSpecies, 
		   const int& numGenes,
		   const std::vector<spParamsF>& popParams,
		   const std::vector<std::vector<myT> >& Genotypes){
  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }
  
  // Resize before return
  //outNS.resize(numSpecies + 2, outNS_i + 1);
  outNS.resize(numSpecies + 1, outNS_i + 1);

  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  
  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  // #ifdef DEBUGW
  // NumericMatrix returnParams(numSpecies, 9);
  // for(int i = 0; i < numSpecies; ++i) {
  //   returnParams(i, 0) = popParams[i].Flag;
  //   returnParams(i, 1) = popParams[i].birth;
  //   returnParams(i, 2) = popParams[i].popSize;
  //   returnParams(i, 3) = popParams[i].timeLastUpdate;
  //   returnParams(i, 4) = popParams[i].W;
  //   returnParams(i, 5) = popParams[i].R;
  //   returnParams(i, 6) = popParams[i].nextMutationTime;
  //   returnParams(i, 7) = popParams[i].mutation;
  //   returnParams(i, 8) = popParams[i].death;

  // }
  // #else
  // std::string returnParams = "NA. Set DEBUGW to return it";
  // #endif
  
  // for(int i = 0; i < popParams.size(); ++i) {
  //   popSizes[i] = popParams[i].popSize;
  // }
}



SEXP Algorithm5L(SEXP restrictTable_,
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
  // Based on 5K, but the min for finding the time is done
  // using another vector just of times. This could be done better,
  // removing time from the struct, but I try version M, with maps
  // which should be a lot faster.

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
  int numSpecies = 0;
  //  int totalNumSpecies = 0;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0; // the column in the outNS
  int sp = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  //double tmpSize = 0.0;
  
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


  //FIXME: nmt
  std::vector<double> nextMutTime(1);
  nextMutTime.reserve(sp_per_period);

  // 5.1 Initialize 
  popParams[0].Flag = true;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
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

    
    for(int i = 0; i < popParams.size(); i++) {
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
	//FIXME:nmt
	nextMutTime[i] = popParams[i].nextMutationTime;

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
    

   
    //FIXME:nmt
    getMinNextMutationTime4(nextMutant, minNextMutationTime, nextMutTime);
    //getMinNextMutationTime3(nextMutant, minNextMutationTime, popParams);
 
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
	nextMutTime.push_back(-99);
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
	remove_zero_sp_v5(sp_to_remove, Genotypes, popParams, nextMutTime);

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


// FIXME00: remove iter from output in NS
SEXP Algorithm5J(SEXP restrictTable_,
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
  // A copy of 5H, with more things turned into functions

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
  int numSpecies = 0;
  int nS_nonZero = 0; // nS_nonZero <= numSpecies;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  //double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  if(numGenes > 62)  
    throw std::range_error("Algo5H: numGenes cannot be larger than 62. Use GMP or non-id'ed version");
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // no GMP but long int
  std::vector<unsigned long> sp_id(1);
  sp_id.reserve(initSp);
  unsigned long new_id; 

  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsF tmpParam; 
  std::vector<spParamsF> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;
  nS_nonZero = 0; // ugly! FIXME: remove this variable?

  popParams[0].Flag = true;
  popParams[0].pos_in_NS = -99;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  
  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // no GMP, but unsigned long
  sp_id[0] = 0;
  

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n       sp_id =  " << sp_id[i] 
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
    

    // minNextMutationTime = std::numeric_limits<double>::infinity();
    // nextMutant = -99; 

    // for(int i = 0; i < popParams.size(); i++) {
    //   if(popParams[i].popSize > 0.0) {
    // 	if(popParams[i].nextMutationTime < minNextMutationTime) {
    // 	  nextMutant = i;
    // 	  minNextMutationTime = popParams[i].nextMutationTime;
    // 	}     
    //   }
    // }

    getMinNextMutationTime(nextMutant, minNextMutationTime, popParams);
 


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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      
      getMutatedPos(mutatedPos, numMutablePosParent, r, mutablePos,
		    nextMutant, Genotypes, numGenes);

      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      new_sp_longint1(sp, new_id, sp_id, numSpecies, nextMutant, 
		      mutatedPos);


      // numSpecies can decrease later if at sampling time 
      // a newly created species has popSize = 0.
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	sp_id.push_back(new_id);

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent sp_id = " << sp_id[nextMutant];
	}
#endif
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePosParent + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePosParent is found BEFORE the mutation
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
	tmpParam.pos_in_NS = -99;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species


#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (sp id = " << sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << sp_id[nextMutant] << ")";
	}

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

      sample_all_pop(nS_nonZero, currentTime, sp_to_remove, 
		     popParams, sp_id, tSample);

      numSpecies = nS_nonZero;
      timeNextPopSample += sampleEvery;
      
      if(sp_to_remove[0])
	remove_zero_sp_v2(sp_to_remove, sp_id, Genotypes, popParams);

      resize_outNS(outNS, outNS_i, numSpecies);
      
      totPopSize_and_fill_outNS(totPopSize, outNS_i, outNS, popParams, //iter
				sp_id, currentTime);

      
      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || (tSample >= finalTime)) 
	simulsDone = true;
	
      forceSample = false;
#ifdef DEBUGV
      std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  NumericVector popSizes(numSpecies);
  
  create_return(popSizes, returnGenotypes, outNS, outNS_i, 
		numSpecies, numGenes, popParams, Genotypes);
  
  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("FinalTime") = currentTime,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}



SEXP Algorithm5I(SEXP restrictTable_,
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
  // A copy of 5H but no use for sp_id at all

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
  int numSpecies = 0;
  int nS_nonZero = 0; // nS_nonZero <= numSpecies;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  if(numGenes > 62) {
    std::cout << "\n With more than 62 genes, expect problems in reporting id_code\n";
  } 

  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // no GMP nor long int
  // std::vector<unsigned long> sp_id(1);
  // sp_id.reserve(initSp);
  // unsigned long new_id; 

  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsF tmpParam; 
  std::vector<spParamsF> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;

  nS_nonZero = 0;
  popParams[0].pos_in_NS = -99;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // no GMP, but unsigned long
  //sp_id[0] = 0;
  

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n       sp_id =  " << id_code(Genotypes[i])  //sp_id[i] 
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
    

    // minNextMutationTime = std::numeric_limits<double>::infinity();
    // nextMutant = -99; 

    // for(int i = 0; i < popParams.size(); i++) {
    //   if(popParams[i].popSize > 0.0) {
    // 	if(popParams[i].nextMutationTime < minNextMutationTime) {
    // 	  nextMutant = i;
    // 	  minNextMutationTime = popParams[i].nextMutationTime;
    // 	}     
    //   }
    // }

    getMinNextMutationTime(nextMutant, minNextMutationTime, popParams);
 
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      getMutatedPos(mutatedPos, numMutablePosParent, r, mutablePos,
		    nextMutant, Genotypes, numGenes);

      // numMutablePos = 0;
      // for(int i=0; i < numGenes; ++i) {
      // 	if(!Genotypes[nextMutant][i]) { 
      // 	  mutablePos[numMutablePos] = i;
      // 	  ++numMutablePos;
      // 	}
      // }
      // if(numMutablePos > 1) {
      // 	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      // } else if (numMutablePos == 1) {
      // 	mutatedPos = mutablePos[0];
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5: run out of mutable places!!??");
      // }


#ifdef DEBUGW
      std::cout << "\n numMutablePosParent = " << numMutablePosParent;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

      // GMP
      // mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      // mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      // now in function
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      // sp = 0;
      
      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      // for(sp = 0; sp < numSpecies; ++sp) {
      // 	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      // 	  break;
      // }

      new_sp_genotypes(sp, numSpecies, numGenes, newGenotype, 
		       Genotypes);


      // numSpecies can decrease later if at sampling time 
      // a newly created species has popSize = 0.
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	//sp_id.push_back(new_id);

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << id_code(Genotypes[sp]); 
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent sp_id = " << id_code(Genotypes[nextMutant]);
	  //sp_id[nextMutant];
	}
#endif
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePosParent + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePosParent is found BEFORE the mutation
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	// if( newMutable[0] != (numGenes - tmpSumMutPos))
	//   throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// tmpParam.mutation = mu * newMutable[0];

	tmpParam.mutation = mu * (numMutablePosParent - 1);
	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));

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
	tmpParam.pos_in_NS = -99;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species


#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (sp id = " << id_code(Genotypes[sp]) << ")" 
	    //sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << id_code(Genotypes[nextMutant]) << ")"; 
	  //sp_id[nextMutant] << ")";
	}

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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;
      sp_to_remove[0] = 0;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	    std::cout << "\n\n     ********* 5.9 ******\n " 
		      << "     Species  = " << i 
		      << "\n      sp_id = " << id_code(Genotypes[i]) //sp_id[i]  
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
	  if( (popParams[i].pos_in_NS == -99) && (popParams[i].popSize <=  0.0) ) {
	    // this i has never been non-zero in any sampling time
	    sp_to_remove[0]++;
	    sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
	    std::cout << "\n\n     Removing species i = " << i 
		      << " with sp_id = " << id_code(Genotypes[i]); //sp_id[i];
#endif
	  } else {
	    if(popParams[i].pos_in_NS == -99) {
	      popParams[i].pos_in_NS = nS_nonZero;	      
	      ++nS_nonZero;
	    }
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
      
      

      numSpecies = nS_nonZero;
      timeNextPopSample += sampleEvery;

      // FIXME00: this function should take sp_to_remove, as
      // I am concerned about size.
      if(sp_to_remove[0])
	remove_zero_sp_v3(sp_to_remove, Genotypes, popParams);
      //FIXME00:
      // place an assert here to check pop sizes initially


      // FIXME00: all this resizing to a function.
      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      // FIXMErm
      // std::cout << "\n numSpecies = " << numSpecies;
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;


#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif

	// FIXME00: turn this into a function for profiling.
	for(int i = 0; i < popParams.size(); ++i) {
	  //	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  // I need one additional subscripting, unneeded if I keep all ever non-zero
	  // as popParams[i].pos_in_NS == i
	  outNS(popParams[i].pos_in_NS + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "\n       Species " << i 
		    << ", sp_id = " << id_code(Genotypes[i]) //sp_id[i] 
		    << ". Pop size = " << popParams[i].popSize ;
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}


SEXP Algorithm5H(SEXP restrictTable_,
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
  // A copy of 5G, but without GMP

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
  int numSpecies = 0;
  int nS_nonZero = 0; // nS_nonZero <= numSpecies;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  if(numGenes > 62)  
    throw std::range_error("Algo5H: numGenes cannot be larger than 62. Use GMP or non-id'ed version");
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // no GMP but long int
  std::vector<unsigned long> sp_id(1);
  sp_id.reserve(initSp);
  unsigned long new_id; 

  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsF tmpParam; 
  std::vector<spParamsF> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;

  nS_nonZero = 0;
  popParams[0].pos_in_NS = -99;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // no GMP, but unsigned long
  sp_id[0] = 0;
  

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n       sp_id =  " << sp_id[i] 
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
    

    // minNextMutationTime = std::numeric_limits<double>::infinity();
    // nextMutant = -99; 

    // for(int i = 0; i < popParams.size(); i++) {
    //   if(popParams[i].popSize > 0.0) {
    // 	if(popParams[i].nextMutationTime < minNextMutationTime) {
    // 	  nextMutant = i;
    // 	  minNextMutationTime = popParams[i].nextMutationTime;
    // 	}     
    //   }
    // }

    getMinNextMutationTime(nextMutant, minNextMutationTime, popParams);
 
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      getMutatedPos(mutatedPos, numMutablePosParent, r, mutablePos,
		    nextMutant, Genotypes, numGenes);

      // numMutablePos = 0;
      // for(int i=0; i < numGenes; ++i) {
      // 	if(!Genotypes[nextMutant][i]) { 
      // 	  mutablePos[numMutablePos] = i;
      // 	  ++numMutablePos;
      // 	}
      // }
      // if(numMutablePos > 1) {
      // 	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      // } else if (numMutablePos == 1) {
      // 	mutatedPos = mutablePos[0];
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5: run out of mutable places!!??");
      // }


#ifdef DEBUGW
      std::cout << "\n numMutablePosParent = " << numMutablePosParent;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

      // GMP
      // mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      // mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      // now in function
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      // sp = 0;
      
      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      // for(sp = 0; sp < numSpecies; ++sp) {
      // 	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      // 	  break;
      // }

      new_sp_longint1(sp, new_id, sp_id, numSpecies, nextMutant, 
		      mutatedPos);


      // numSpecies can decrease later if at sampling time 
      // a newly created species has popSize = 0.
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	sp_id.push_back(new_id);

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent sp_id = " << sp_id[nextMutant];
	}
#endif
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePosParent + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePosParent is found BEFORE the mutation
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	// if( newMutable[0] != (numGenes - tmpSumMutPos))
	//   throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// tmpParam.mutation = mu * newMutable[0];

	tmpParam.mutation = mu * (numMutablePosParent - 1);
	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));

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
	tmpParam.pos_in_NS = -99;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species


#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (sp id = " << sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << sp_id[nextMutant] << ")";
	}

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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;
      sp_to_remove[0] = 0;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	    std::cout << "\n\n     ********* 5.9 ******\n " 
		      << "     Species  = " << i 
		      << "\n      sp_id = " << sp_id[i]  
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
	  if( (popParams[i].pos_in_NS == -99) && (popParams[i].popSize <=  0.0) ) {
	    // this i has never been non-zero in any sampling time
	    sp_to_remove[0]++;
	    sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
	    std::cout << "\n\n     Removing species i = " << i 
		      << " with sp_id = " << sp_id[i];
#endif
	  } else {
	    if(popParams[i].pos_in_NS == -99) {
	      popParams[i].pos_in_NS = nS_nonZero;	      
	      ++nS_nonZero;
	    }
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
      
      

      numSpecies = nS_nonZero;
      timeNextPopSample += sampleEvery;

      if(sp_to_remove[0])
	remove_zero_sp_v2(sp_to_remove, sp_id, Genotypes, popParams);
      //FIXME00:
      // place an assert here to check pop sizes initially

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      // FIXMErm
      // std::cout << "\n numSpecies = " << numSpecies;
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;


#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	for(int i = 0; i < popParams.size(); ++i) {
	  //	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  // I need one additional subscripting, unneeded if I keep all ever non-zero
	  // as popParams[i].pos_in_NS == i
	  outNS(popParams[i].pos_in_NS + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "\n       Species " << i 
		    << ", sp_id = " << sp_id[i] 
		    << ". Pop size = " << popParams[i].popSize ;
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}



SEXP Algorithm5G(SEXP restrictTable_,
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
  // A copy of 5F, but with some parts turned into functions
  // for profiling and clearer code.
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
  int numSpecies = 0;
  int nS_nonZero = 0; // nS_nonZero <= numSpecies;
  int iter = 0;
  int numMutablePosParent = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsF tmpParam; 
  std::vector<spParamsF> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;

  nS_nonZero = 0;
  popParams[0].pos_in_NS = -99;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n       sp_id =  " << sp_id[i] 
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
    

    // minNextMutationTime = std::numeric_limits<double>::infinity();
    // nextMutant = -99; 

    // for(int i = 0; i < popParams.size(); i++) {
    //   if(popParams[i].popSize > 0.0) {
    // 	if(popParams[i].nextMutationTime < minNextMutationTime) {
    // 	  nextMutant = i;
    // 	  minNextMutationTime = popParams[i].nextMutationTime;
    // 	}     
    //   }
    // }

    getMinNextMutationTime(nextMutant, minNextMutationTime, popParams);
 
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      getMutatedPos(mutatedPos, numMutablePosParent, r, mutablePos,
		    nextMutant, Genotypes, numGenes);

      // numMutablePos = 0;
      // for(int i=0; i < numGenes; ++i) {
      // 	if(!Genotypes[nextMutant][i]) { 
      // 	  mutablePos[numMutablePos] = i;
      // 	  ++numMutablePos;
      // 	}
      // }
      // if(numMutablePos > 1) {
      // 	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      // } else if (numMutablePos == 1) {
      // 	mutatedPos = mutablePos[0];
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5: run out of mutable places!!??");
      // }


#ifdef DEBUGW
      std::cout << "\n numMutablePosParent = " << numMutablePosParent;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

      // GMP
      // mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      // mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      // now in function
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      // sp = 0;
      
      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      // for(sp = 0; sp < numSpecies; ++sp) {
      // 	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      // 	  break;
      // }

      new_sp_gmp(sp, add_to_id, new_id, sp_id, numSpecies, nextMutant, 
		 mutatedPos);

      // numSpecies can decrease later if at sampling time 
      // a newly created species has popSize = 0.
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}

#ifdef DEBUGV	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n parent sp_id = " << sp_id[nextMutant];
	}
#endif
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePosParent + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePosParent is found BEFORE the mutation
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	// if( newMutable[0] != (numGenes - tmpSumMutPos))
	//   throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// tmpParam.mutation = mu * newMutable[0];

	tmpParam.mutation = mu * (numMutablePosParent - 1);
	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));

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
	tmpParam.pos_in_NS = -99;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species


#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (sp id = " << sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << sp_id[nextMutant] << ")";
	}

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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;
      sp_to_remove[0] = 0;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	    std::cout << "\n\n     ********* 5.9 ******\n " 
		      << "     Species  = " << i 
		      << "\n      sp_id = " << sp_id[i]  
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
	  if( (popParams[i].pos_in_NS == -99) && (popParams[i].popSize <=  0.0) ) {
	    // this i has never been non-zero in any sampling time
	    sp_to_remove[0]++;
	    sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
	    std::cout << "\n\n     Removing species i = " << i 
		      << " with sp_id = " << sp_id[i];
#endif
	  } else {
	    if(popParams[i].pos_in_NS == -99) {
	      popParams[i].pos_in_NS = nS_nonZero;	      
	      ++nS_nonZero;
	    }
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
      
      

      numSpecies = nS_nonZero;
      timeNextPopSample += sampleEvery;

      if(sp_to_remove[0])
	remove_zero_sp_v1(sp_to_remove, sp_id, Genotypes, popParams);
      //FIXME00:
      // place an assert here to check pop sizes initially

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      // FIXMErm
      // std::cout << "\n numSpecies = " << numSpecies;
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;


#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	for(int i = 0; i < popParams.size(); ++i) {
	  //	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  // I need one additional subscripting, unneeded if I keep all ever non-zero
	  // as popParams[i].pos_in_NS == i
	  outNS(popParams[i].pos_in_NS + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "\n       Species " << i 
		    << ", sp_id = " << sp_id[i] 
		    << ". Pop size = " << popParams[i].popSize ;
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}





// As of 2013-03-04: Algorithm 5F, 5E, and 5D do what they are supposed to
// do. Note all use Bozic. Also, 5F and 5E do not give the same results,
// because removing species means that, say, sp. 2 can disappear in, say,
// iter 4, and be created again in iter 23. So the order in which ti is
// evaluated changes, and thus the sequence of random numbers.


// FIXME: this is thought, as of now, just for Bozic.
// This is VERY ugly, and just to get going.
// Later, do it well. Make w_f, r_f, algo2, algo3, ti_next, and some others
// take a pointer to struct, so we only care about bozic or linear
// when we return the birth rate.

SEXP Algorithm5F(SEXP restrictTable_,
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
  int numSpecies = 0;
  int nS_nonZero = 0; // nS_nonZero <= numSpecies;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  
  double tmpSize = 0.0;
  
  int max_remove = 1000000;
  std::vector<int>sp_to_remove(max_remove + 1);
  
  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsF tmpParam; 
  std::vector<spParamsF> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;

  nS_nonZero = 0;
  popParams[0].pos_in_NS = -99;

  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n       sp_id =  " << sp_id[i] 
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
    

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].popSize > 0.0) {
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }

 
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
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // Mutables
      // numMutablePos = Mutables[nextMutant][0];

      // if(numMutablePos > 1) {
      // 	indexMutatedPos =  1 + 
      // 	  static_cast<int>((gsl_rng_uniform_int(r, numMutablePos)));
      // } else if(numMutablePos == 1) {
      // 	indexMutatedPos = 1;
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5B: run out of mutable places!!??");
      // }

      // mutatedPos = Mutables[nextMutant][indexMutatedPos];
      // newMutable = Mutables[nextMutant];
      // newMutable.erase(newMutable.begin() + indexMutatedPos); 
      // newMutable[0] = numMutablePos - 1;
      
      // if(verbosity && (newMutable[0] == 0)) {
      // 	std::cout << "\n A species run out of mutable genes \n";
      // }

      // NOT using mutables
      // Uncomment if not using Mutables
      numMutablePos = 0;
      for(int i=0; i < numGenes; ++i) {
      	if(!Genotypes[nextMutant][i]) { 
      	  mutablePos[numMutablePos] = i;
      	  ++numMutablePos;
      	}
      }
      if(numMutablePos > 1) {
      	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      } else if (numMutablePos == 1) {
      	mutatedPos = mutablePos[0];
      } else {
      	// Should never happen, as mutation = 0 if no mutable positions.
      	throw std::out_of_range("Algo5: run out of mutable places!!??");
      }

#ifdef DEBUGW
      std::cout << "\n numMutablePos = " << numMutablePos;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

      // GMP
      mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      sp = 0;
      k = 0;

      
      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      for(sp = 0; sp < numSpecies; ++sp) {
      	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      	  break;
      }
      // numSpecies can decrease later if at sampling time 
      // a newly created species has popSize = 0.
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n patern sp_id = " << sp_id[nextMutant];
	}
	
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePos is found BEFORE the mutation
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	// if( newMutable[0] != (numGenes - tmpSumMutPos))
	//   throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// tmpParam.mutation = mu * newMutable[0];

	tmpParam.mutation = mu * (numMutablePos - 1);
	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));

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
	tmpParam.pos_in_NS = -99;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species


#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " (sp id = " << sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << sp_id[nextMutant] << ")";
	}

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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;
      sp_to_remove[0] = 0;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	    std::cout << "\n\n     ********* 5.9 ******\n " 
		      << "     Species  = " << i 
		      << "\n      sp_id = " << sp_id[i]  
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
	  if( (popParams[i].pos_in_NS == -99) && (popParams[i].popSize <=  0.0) ) {
	    // this i has never been non-zero in any sampling time
	    sp_to_remove[0]++;
	    sp_to_remove[sp_to_remove[0]] = i;
#ifdef DEBUGV
	    std::cout << "\n\n     Removing species i = " << i 
		      << " with sp_id = " << sp_id[i];
#endif
	  } else {
	    if(popParams[i].pos_in_NS == -99) {
	      popParams[i].pos_in_NS = nS_nonZero;	      
	      ++nS_nonZero;
	    }
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
      
      

      numSpecies = nS_nonZero;
      timeNextPopSample += sampleEvery;

      if(sp_to_remove[0])
	remove_zero_sp_v1(sp_to_remove, sp_id, Genotypes, popParams);
      //FIXME00:
      // place an assert here to check pop sizes initially

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      // FIXMErm
      // std::cout << "\n numSpecies = " << numSpecies;
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;


#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	for(int i = 0; i < popParams.size(); ++i) {
	  //	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  // I need one additional subscripting, unneeded if I keep all ever non-zero
	  // as popParams[i].pos_in_NS == i
	  outNS(popParams[i].pos_in_NS + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "\n       Species " << i 
		    << ", sp_id = " << sp_id[i] 
		    << ". Pop size = " << popParams[i].popSize ;
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
    }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}



// FIXME: this is thought, as of now, just for Bozic.
// This is VERY ugly, and just to get going.
// Later, do it well. Make w_f, r_f, algo2, algo3, ti_next, and some others
// take a pointer to struct, so we only care about bozic or linear
// when we return the birth rate.

SEXP Algorithm5E(SEXP restrictTable_,
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
  int numSpecies = 0;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  
  // needed if not using Mutables
  std::vector<int>mutablePos(numGenes);


  // We get rid of mutables
  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  // std::vector<myT> newMutable(numGenes + 1);
  // std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  // Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsD tmpParam; 
  std::vector<spParamsD> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;
  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  // Mutables[0][0] = numGenes;
  // for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL ))
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
		    << "\n      sp_id =  " << sp_id[i] 
		    << "\n      popSize = " << popParams[i].popSize 
		    << "\n      currentTime = " << currentTime 
		    << "\n      popParams[i].nextMutationTime = " 
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
    

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].popSize > 0.0) {
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }

 
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
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // Mutables
      // numMutablePos = Mutables[nextMutant][0];

      // if(numMutablePos > 1) {
      // 	indexMutatedPos =  1 + 
      // 	  static_cast<int>((gsl_rng_uniform_int(r, numMutablePos)));
      // } else if(numMutablePos == 1) {
      // 	indexMutatedPos = 1;
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5B: run out of mutable places!!??");
      // }

      // mutatedPos = Mutables[nextMutant][indexMutatedPos];
      // newMutable = Mutables[nextMutant];
      // newMutable.erase(newMutable.begin() + indexMutatedPos); 
      // newMutable[0] = numMutablePos - 1;
      
      // if(verbosity && (newMutable[0] == 0)) {
      // 	std::cout << "\n A species run out of mutable genes \n";
      // }

      // NOT using mutables
      // Uncomment if not using Mutables
      numMutablePos = 0;
      for(int i=0; i < numGenes; ++i) {
      	if(!Genotypes[nextMutant][i]) { 
      	  mutablePos[numMutablePos] = i;
      	  ++numMutablePos;
      	}
      }
      if(numMutablePos > 1) {
      	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      } else if (numMutablePos == 1) {
      	mutatedPos = mutablePos[0];
      } else {
      	// Should never happen, as mutation = 0 if no mutable positions.
      	throw std::out_of_range("Algo5: run out of mutable places!!??");
      }

#ifdef DEBUGW
      std::cout << "\n numMutablePos = " << numMutablePos;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
      
#endif

      // GMP
      mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      sp = 0;
      k = 0;

      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      for(sp = 0; sp < numSpecies; ++sp) {
      	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      	  break;
      }
      
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	//Mutables.push_back(newMutable);
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	  std::cout << "\n patern sp_id = " << sp_id[nextMutant];
	}
	
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// that expression IS correct: numMutablePos is found BEFORE the mutation
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	// if( newMutable[0] != (numGenes - tmpSumMutPos))
	//   throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// tmpParam.mutation = mu * newMutable[0];

	tmpParam.mutation = mu * (numMutablePos - 1);
	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));

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
		    << " (sp id = " << sp_id[sp] << ")"
		    << "\n from species "  <<   nextMutant
		    << " (sp_id = " << sp_id[nextMutant] << ")";
	}

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
	  if(verbosity >= 2) {
	    std::cout << "\n Mutation to an extinct species\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
	// popParams[sp].timeLastUpdate = currentTime; //not explicitly said in paper
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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
#ifdef DEBUGV
	    std::cout << "\n\n     ********* 5.9 ******\n " 
		      << "     Species  = " << i 
		      << "\n      sp_id = " << sp_id[i]  
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
	  // forceSampling, popSize was updated in previous loop.
	  if(tSample > popParams[i].timeLastUpdate) 
	    popParams[i].popSize = 
	      Algo2(popParams[i].popSize,
		    tSample - popParams[i].timeLastUpdate,
		    popParams[i].R,
		    popParams[i].W,
		    popParams[i].death,
		    popParams[i].birth);

	  popParams[i].Flag = true;

#ifdef DEBUGV
	    std::cout << "\n\n   post-update popSize = " 
		      << popParams[i].popSize << "\n";
#endif

	  
#ifdef DEBUGW	  
	  popParams[i].timeLastUpdate = -999999.99999;
#endif
	}
      }

      timeNextPopSample += sampleEvery;

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
totPopSize = 0.0;

#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	
	for(int i = 0; i < popParams.size(); ++i) {
	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "\n       Pop " << i 
		    << ", sp_id = " << sp_id[i] 
		    << ". Pop size = " << popParams[i].popSize ;
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
	}
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}


// FIXME: this is thought, as of now, just for Bozic.
// This is VERY ugly, and just to get going.
// Later, do it well. Make w_f, r_f, algo2, algo3, ti_next, and some others
// take a pointer to struct, so we only care about bozic or linear
// when we return the birth rate.

SEXP Algorithm5D(SEXP restrictTable_,
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
  int numSpecies = 0;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  

  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  std::vector<myT> newMutable(numGenes + 1);
  std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsD tmpParam; 
  std::vector<spParamsD> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;
  // FIXME??
  if(typeFitness == "bozic") {
    popParams[0].birth = 0.5;
    popParams[0].death = 0.5;
  } else {    
    popParams[0].birth = birthRate;
    popParams[0].death = death;
  }

  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(popParams[0].death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  Mutables[0][0] = numGenes;
  for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL ))
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
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
	  std::cout << "\n\n     ********* 5.2: call to ti_nextTime ******\n " <<
	    "     Species  = " << i <<
	    "\n   popSize = " << popParams[i].popSize << 
	    "\n   currentTime = " << currentTime <<
	    "; popParams[i].nextMutationTime = " << popParams[i].nextMutationTime;
	    std::cout << " \n     species R " << popParams[i].R;
	  std::cout << " \n     species W " << popParams[i].W;
	  std::cout << " \n     species death " << popParams[i].death;
	  std::cout << " \n     species birth " << popParams[i].birth;
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
    

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].popSize > 0.0) {
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }

 
    if(verbosity >= 2) {
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;


#ifdef DEBUGV
      std::cout << "\n DEBUGV In 5.3  currentTime " << currentTime <<"\n";

#endif
      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
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
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > -2) // we always warn about this
 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // Mutables
      numMutablePos = Mutables[nextMutant][0];

      if(numMutablePos > 1) {
      	indexMutatedPos =  1 + 
      	  static_cast<int>((gsl_rng_uniform_int(r, numMutablePos)));
      } else if(numMutablePos == 1) {
      	indexMutatedPos = 1;
      } else {
      	// Should never happen, as mutation = 0 if no mutable positions.
      	throw std::out_of_range("Algo5B: run out of mutable places!!??");
      }

      mutatedPos = Mutables[nextMutant][indexMutatedPos];
      newMutable = Mutables[nextMutant];
      newMutable.erase(newMutable.begin() + indexMutatedPos); 
      newMutable[0] = numMutablePos - 1;

#ifdef DEBUGW
      std::cout << "\n numMutablePos = " << numMutablePos;
      //      std::cout << "\n indexMutatedPos = " << indexMutatedPos;
      std::cout << "\n mutatedPos = " << mutatedPos  << "\n";
#endif

      
      if(verbosity && (newMutable[0] == 0)) {
	std::cout << "\n A species run out of mutable genes \n";
      }


      // GMP
      mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      sp = 0;
      k = 0;

      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      for(sp = 0; sp < numSpecies; ++sp) {
      	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      	  break;
      }
      
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	Mutables.push_back(newMutable);
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;

	if(typeFitness == "bozic") {
	  // if bozic, always death, so we pass death of parent
	  tmpParam.death = fitness_CBN(mutatedPos,
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
	  tmpParam.birth = fitness_CBN_std(mutatedPos,
					   restrictTable,
					   popParams[nextMutant].birth,
					   typeCBN,
					   newGenotype,
					   birthRate,
					   s,
					   numDrivers);
	  tmpParam.death = death;
	}
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birth of sp = " << tmpParam.birth;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	}
	
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	if( newMutable[0] != (numGenes - tmpSumMutPos))
	  throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	tmpParam.mutation = mu * newMutable[0];


#ifdef DEBUGV
	std::cout << "\n DEBUGV mu " << mu;	
	std::cout << "\n DEBUGV tmpParam.mutation " << tmpParam.mutation;
	std::cout << "\n DEBUGV newMutable[0] or numMutablePos  = " << newMutable[0];	
	
#endif


	
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

#ifdef DEBUGV
      std::cout << "\n DEBUGV Before crash currentTime " << currentTime <<"\n";
      std::cout << "\n DEBUGV Before crash timeLastUpdate " 
		<< popParams[sp].timeLastUpdate <<"\n";

#endif

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " from species "  <<   mutatedPos 
		    << ". New popSize = " << popParams[sp].popSize;
	}

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  popParams[sp].death,
		  popParams[sp].birth);
	} else {
	  if(verbosity >= 2) {
	    std::cout << "\n             Mutation " << 
	      "to an already extinct species, sp " <<
	      sp <<"\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
	// popParams[sp].timeLastUpdate = currentTime; //not explicitly said in paper
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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.9 ******\n " <<
	    "     Species  = " << i <<
	    "\n   pre-update popSize = " << popParams[i].popSize << 
	    "\n   time of sample = " << tSample <<
	    "; popParams[i].timeLastUpdate = " << popParams[i].timeLastUpdate <<
	    "; t for Algo2 = " << tSample - popParams[i].timeLastUpdate <<"\n";
	  std::cout << " \n     species R " << popParams[i].R;
	  std::cout << " \n     species W " << popParams[i].W;
	  std::cout << " \n     species death " << popParams[i].death;
	  std::cout << " \n     species birth " << popParams[i].birth;
	  std::cout << " \n     species nextMutationTIme " << popParams[i].nextMutationTime;
#endif
	  // Account for forceSampling. When 
	  // forceSampling, popSize was updated in previous loop.
	  if(tSample > popParams[i].timeLastUpdate) 
	    popParams[i].popSize = 
	      Algo2(popParams[i].popSize,
		    tSample - popParams[i].timeLastUpdate,
		    popParams[i].R,
		    popParams[i].W,
		    popParams[i].death,
		    popParams[i].birth);

	  popParams[i].Flag = true;

#ifdef DEBUGV
	  std::cout << "\n\n   post-update popSize = " << 
	    popParams[i].popSize << "\n";
#endif

	  
#ifdef DEBUGW	  
	  popParams[i].timeLastUpdate = -999999.99999;
	  // popParams[i].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	  }
	  //	}
      }

      timeNextPopSample += sampleEvery;

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
totPopSize = 0.0;

#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	
	for(int i = 0; i < popParams.size(); ++i) {
	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "       Pop size pop " << i << " = " << popParams[i].popSize << "\n";
#endif
	}
	
#ifdef DEBUGV
	std::cout << "\n\n       totPopSize   = " << totPopSize << "\n";
#endif
	
	if( !std::isfinite(totPopSize) ) {
	  throw std::range_error("totPopSize not finite");
	}
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling forceSampling is " << forceSample <<"\n";
#endif 
	  }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
    returnParams(i, 8) = popParams[i].death;

  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("iter") = iter,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}




SEXP Algorithm5C(SEXP restrictTable_,
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
		 SEXP ratioForce_) {
  BEGIN_RCPP
  using namespace Rcpp;
    
  const IntegerMatrix restrictTable(restrictTable_);
  const int numDrivers = as<int>(numDrivers_);
  const int numGenes = as<int>(numGenes_);
  const std::string typeCBN = as<std::string>(typeCBN_);
  const double birthRate = as<double>(birthRate_);
  const double s = as<double>(s_);
  const double death = as<double>(death_);
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
  int numSpecies = 0;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; 
  int timeL = 1000;
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  

  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 (or 127?)

  std::vector<myT> newMutable(numGenes + 1);
  std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParamsC tmpParam; 
  std::vector<spParamsC> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;
  popParams[0].birth = birthRate;
  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(death, popParams[0].birth, popParams[0].mutation);

  timeNextPopSample = currentTime + sampleEvery;

  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  Mutables[0][0] = numGenes;
  for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;


  while(!simulsDone) {
    iter++;
    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL ))
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
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

    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
	popParams[i].nextMutationTime = ti_nextTime_tmax(popParams[i].R,
							 popParams[i].W,
							 death,
							 popParams[i].birth,
							 popParams[i].popSize,
							 currentTime,
							 tSample);
	popParams[i].Flag = false;
	popParams[i].timeLastUpdate = currentTime;
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
    

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].popSize > 0.0) {
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }

 
    if(verbosity >= 2) {
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeNextPopSample = " << timeNextPopSample << "\n";
    }
    
    // Do we need to sample the population?
    if( minNextMutationTime <= tSample ) {// We are not sampling
      // ************   5.3   **************
      currentTime = minNextMutationTime;


#ifdef DEBUGV
      std::cout << "\n DEBUGV In 5.3  currentTime " << currentTime <<"\n";

#endif
      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	std::cout << "ERROR!! mutantTimeSinceLastUpdate " << 
	  mutantTimeSinceLastUpdate  << "  sampleEvery = " << sampleEvery << 
	  "  currentTime " << currentTime << " popParams[nextMutant].timeLastUpdate " <<
	  popParams[nextMutant].timeLastUpdate << 
	  " nextMutant = " << nextMutant << "\n";
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    death,
					    popParams[nextMutant].birth);
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	forceSample = true;
	ratioForce = std::min(1.0, 2 * ratioForce);
	if(verbosity > 0) {
	  std::cout << "\n Forced sampling triggered for next loop: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant <<
	    " at iteration " << iter << "\n";
	}
      }      
      
      // Check also for numSpecies, and force sampling if needed
      if(! (numSpecies % speciesFS )) {
      	forceSample = true;
	speciesFS *= 2; 
      	if(verbosity > 0) 
      	  std::cout << "\n Forced sampling triggered for next loop "
		    << " when numSpecies = " << 
      	    numSpecies << " at iteration " << iter << "\n";
      }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // Mutables
      numMutablePos = Mutables[nextMutant][0];

      if(numMutablePos > 1) {
      	indexMutatedPos =  1 + 
      	  static_cast<int>((gsl_rng_uniform_int(r, numMutablePos)));
      } else if(numMutablePos == 1) {
      	indexMutatedPos = 1;
      } else {
      	// Should never happen, as mutation = 0 if no mutable positions.
      	throw std::out_of_range("Algo5B: run out of mutable places!!??");
      }

      mutatedPos = Mutables[nextMutant][indexMutatedPos];
      newMutable = Mutables[nextMutant];
      newMutable.erase(newMutable.begin() + indexMutatedPos); 
      newMutable[0] = numMutablePos - 1;
      
      if(verbosity && (newMutable[0] == 0)) {
	std::cout << "\n A species run out of mutable genes \n";
      }


      // GMP
      mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());
      
      
      // ************   5.6   ***************
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      sp = 0;
      k = 0;

      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      // Did we create a new species?
      for(sp = 0; sp < numSpecies; ++sp) {
      	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      	  break;
      }
      
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	Mutables.push_back(newMutable);
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;
	tmpParam.birth = fitness_CBN_std(mutatedPos,
					 restrictTable,
					 popParams[nextMutant].birth,
					 typeCBN,
					 newGenotype,
					 birthRate,
					 s,
					 numDrivers);
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birthRate = " << birthRate;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	}
	
	
#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
	if( newMutable[0] != (numGenes - tmpSumMutPos))
	  throw std::out_of_range(" newMutable[0] != (numGenes - tmpSumMutPos)");
#endif
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	tmpParam.mutation = mu * newMutable[0];

	
	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(death, tmpParam.birth, tmpParam.mutation);
	}
#ifdef DEBUGW
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors
#endif
	tmpParam.Flag = true;
	popParams.push_back(tmpParam);
	// tis and nextMutationTime will be update above, as flagged
      } else {	// A mutation to pre-existing species

#ifdef DEBUGV
      std::cout << "\n DEBUGV Before crash currentTime " << currentTime <<"\n";
      std::cout << "\n DEBUGV Before crash timeLastUpdate " 
		<< popParams[sp].timeLastUpdate <<"\n";

#endif

#ifdef DEBUGW
	if( (currentTime - popParams[sp].timeLastUpdate) <= 0.0)
	  throw std::out_of_range("currentTime - timeLastUpdate out of range"); 
#endif
	
	if(verbosity >= 2) {
	  std::cout <<"\n     Mutated to existing species " << sp 
		    << " from species "  <<   mutatedPos 
		    << ". New popSize = " << popParams[sp].popSize;
	}

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  death,
		  popParams[sp].birth);
	} else {
	  if(verbosity >= 2) {
	    std::cout << "\n             Mutation " << 
	      "to an already extinct species, sp " <<
	      sp <<"\n";
	  }
	  popParams[sp].popSize = 1.0;
	}
	
	
#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
	// popParams[sp].timeLastUpdate = currentTime; //not explicitly said in paper
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
	  // STOPASSERT(tSample == timeNextPopSample);
	} else
	  std::cout <<". We reached finalTime " << finalTime << "\n";
      }
      
      currentTime = tSample;

      for(int i = 0; i < popParams.size(); i++) {
	if(popParams[i].popSize > 0.0) {
	  STOPASSERT(popParams[i].Flag == false);
	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(tSample - popParams[i].timeLastUpdate >= 0.0);
	  
	  
#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.9 ******\n " <<
	    "     Species  = " << i <<
	    "\n   pre-update popSize = " << popParams[i].popSize << 
	    "\n   time of sample = " << tSample <<
	    "; popParams[i].timeLastUpdate = " << popParams[i].timeLastUpdate <<
	    "; t for Algo2 = " << tSample - popParams[i].timeLastUpdate <<"\n";
	  std::cout << " \n     species R " << popParams[i].R;
	  std::cout << " \n     species W " << popParams[i].W;
	  std::cout << " \n     species death " << death;
	  std::cout << " \n     species birth " << popParams[i].birth;
	  std::cout << " \n     species nextMutationTIme " << popParams[i].nextMutationTime;
#endif
	  // Account for forceSampling. When 
	  // forceSampling, popSize was updated in previous loop.
	  if(tSample > popParams[i].timeLastUpdate) 
	    popParams[i].popSize = 
	      Algo2(popParams[i].popSize,
		    tSample - popParams[i].timeLastUpdate,
		    popParams[i].R,
		    popParams[i].W,
		    death,
		    popParams[i].birth);

	  popParams[i].Flag = true;

#ifdef DEBUGV
	  std::cout << "\n\n   post-update popSize = " << 
	    popParams[i].popSize << "\n";
#endif

	  
#ifdef DEBUGW	  
	  popParams[i].timeLastUpdate = -999999.99999;
	  // popParams[i].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	  }
	  //	}
      }

      timeNextPopSample += sampleEvery;

      type_resize = 0;
      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      if( outNS.n_cols <= (outNS_i + 1) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
totPopSize = 0.0;

#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	
	for(int i = 0; i < popParams.size(); ++i) {
	  outNS(i + 2, outNS_i) = popParams[i].popSize;
	  totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	  std::cout << "       Pop size pop " << i << " = " << popParams[i].popSize << "\n";
#endif
	}
	
	
	if( (totPopSize >= detectionSize) ||
	    (totPopSize <= 0.0) || (tSample >= finalTime)) 
	  simulsDone = true;
	
	forceSample = false;
#ifdef DEBUGV
	std::cout << "\n at     end of sampling froceSampliong is " << forceSample <<"\n";
#endif 
	  }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].nextMutationTime;
    returnParams(i, 7) = popParams[i].mutation;
  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}





SEXP Algorithm5B(SEXP restrictTable_,
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
		 SEXP initSize_species_,
		 SEXP initSize_iter_,
		 SEXP seed_gsl_,
		 SEXP verbose_) {
  BEGIN_RCPP
  using namespace Rcpp;
  
  std::cout << "\n WARNING: you have a good reason for using this version of Algo5, right?\n";


  //  test_gmp();

  //  Timer timer;
    
  const IntegerMatrix restrictTable(restrictTable_);
  const int numDrivers = as<int>(numDrivers_);
  const int numGenes = as<int>(numGenes_);
  const std::string typeCBN = as<std::string>(typeCBN_);
  const double birthRate = as<double>(birthRate_);
  const double s = as<double>(s_);
  const double death = as<double>(death_);
  const double mu = as<double>(mu_);
  const double initSize = as<double>(initSize_);
  const double sampleEvery = as<double>(sampleEvery_);
  const double detectionSize = as<double>(detectionSize_);
  const int initSp = as<int>(initSize_species_);
  const int initIt = as<int>(initSize_iter_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  const double ratioForce = 0.2; // If a single species this times
  // detectionSize, force a sampling to prevent going too far.
  const int seed = as<int>(seed_gsl_);


  
  // needed if not using Mutables
  //  std::vector<int>mutablePos(numGenes);

  bool Extinction = false;
  bool forceSample = false;
  bool simulsDone = false;

  double totPopSize = 0;
  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double currentTime = 0.0;
  double timeAllPopSample;

  int nextMutant;
  int numSpecies = 0;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int indexMutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  int iterL = 1000;
  int speciesL = 1000; // used also to force sampling
  int timeL = 1000;
  int speciesFS = 2000;
  

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // GMP
  // I need to use GMP because not a long enough type
  // for an id of > 64 genes.

  std::vector<mpz_class> sp_id(1);
  sp_id.reserve(initSp);
  mpz_t new_id; mpz_init (new_id);
  mpz_t add_to_id; mpz_init (add_to_id);
  

  // Mutables: in each row, first pos is the number of mutable genes left
  // remaing are the actual positions of the mutable genes.
  // myT cannot be a bool nor a char if more than 255 or 127

  std::vector<myT> newMutable(numGenes + 1);
  std::vector<std::vector<myT> > Mutables(1, std::vector<myT>(numGenes + 1));
  Mutables.reserve(initSp);
  

  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParams tmpParam; 
  std::vector<spParams> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.


  //timer.step("pre-init");

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;
  popParams[0].birth = birthRate;
  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(death, popParams[0].birth, popParams[0].mutation);

  timeAllPopSample = currentTime + sampleEvery;
  
  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  
  // GMP
  mpz_set_ui(sp_id[0].get_mpz_t(), 0);

  // Mutables
  Mutables[0][0] = numGenes;
  for(int i = 0; i < numGenes; ++i) Mutables[0][i + 1] = i;

  //timer.step("init");
  

  while(!simulsDone) {
    iter++;

    if(verbosity) {
      if(! (iter % iterL) ) {
	std::cout << "\n\n    ... iteration " << iter;
      // if(!(currentTime % timeL ))
	std::cout << "\n    ... currentTime " << currentTime <<"\n";
      }
      if(!(numSpecies % speciesL ))
	std::cout << "\n\n    ... numSpecies " << numSpecies << "\n";
    }


    //  ************   5.2   ***************

    if(verbosity >= 2)
      std::cout <<"\n\n\n***** Looping again through 5.2 \n";
#ifdef DEBUGV
    std::cout << " \n\n currentTime  " << currentTime;
#endif
    
    for(int i = 0; i < popParams.size(); i++) {
      if((popParams[i].Flag) && (popParams[i].popSize > 0.0))  {
#ifdef DEBUGV
	std::cout << " \n\n species i " << i;
	std::cout << " \n species R " << popParams[i].R;
	std::cout << " \n species W " << popParams[i].W;
	std::cout << " \n species death " << death;
	std::cout << " \n species birth " << popParams[i].birth;
	std::cout << " \n species popSize " << popParams[i].popSize;
#endif

	popParams[i].tis = ti_f(popParams[i].R,
				popParams[i].W,
				death,
				popParams[i].birth,
				popParams[i].popSize);
	popParams[i].nextMutationTime = popParams[i].tis + currentTime;
#ifdef DEBUGV
	 std::cout << " \n species tis " << popParams[i].tis;
	 std::cout << " \n species nextMutTIme " << popParams[i].nextMutationTime;
#endif
	popParams[i].Flag = false;
	popParams[i].timeLastUpdate = currentTime;
// #ifdef DEBUGW
// 	std::cout <<" in 5.2, species i " << i << "  currentTime = " << 
// 	  currentTime <<"\n";
// #endif
      }
    }

    #ifdef DEBUGW
    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].tis < 0.0) throw std::range_error("Algo5: tis < 0");
      else if(popParams[i].tis == 0.0) 
	std::cout << "Algo 5: BEWARE some tis == 0.0" << std::endl;
    }
    #endif
    
    // ******************** 5.3 and do we sample? *********** 

    // We need to find minimum to know if we need to sample the whole pop
    

    // Use a priority queue?

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    // DELETE
    // std::cout << " \n numSpecies = " << numSpecies << " numGenes = " << numGenes;
    // std::cout << " \n sizeGenotypes " << Genotypes.size();

    for(int i = 0; i < popParams.size(); i++) {
      // std::cout << " \n\n species = " << i;
      if(popParams[i].popSize > 0.0) {
	// std::cout << " \n    and it exists with nextmtime " << popParams[i].nextMutationTime;
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }

    if(nextMutant < -90) {
      if(verbosity)
	std::cout << "\n Algo5: population became extinct\n";
      Extinction = true;
      // Note: the algorithm says they become extinct before a mutation.
      // (as per eq. 12 in p. 1232 of paper) But we will sample once, and
      // then set final sizes (after sample) to 0.  However, a pop. with
      // no holes for mutation left would have a ti of infinity (see
      // Algo2) but could keep growing indefinitely. Signal if any pop has
      // not mut. holes left FIXME!!
    }

    if(verbosity >= 2) {
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeAllPopSample = " << timeAllPopSample << "\n";
    }

    // Do we need to sample the population?
    if((timeAllPopSample >= minNextMutationTime) && !forceSample) {// We are not sampling
      // ************   5.3   **************

      currentTime = minNextMutationTime;


      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    death,
					    popParams[nextMutant].birth);
      // Next line not in paper. It will be updated
      // in Algo2, as flagged.
      // popParams[nextMutant].timeLastUpdate = currentTime;
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	// What should we do here???
	// Allow bailing out even if tis become vanishingly small
	forceSample = true;
	if(verbosity) {
	  std::cout << "\n Forced sampling: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant << "\n";
	}
      }      
      
      // Check also for numSpecies, and force sampling if needed
      // if(! (numSpecies % speciesFS )) {
      // 	forceSample = true;
      // 	if(verbosity) 
      // 	  std::cout << "\n Forced sampling when numSpecies = " << 
      // 	    numSpecies << "\n";
      // }
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // Mutables
      numMutablePos = Mutables[nextMutant][0];

      if(numMutablePos > 1) {
      	indexMutatedPos =  1 + 
      	  static_cast<int>((gsl_rng_uniform_int(r, numMutablePos)));
      } else if(numMutablePos == 1) {
      	indexMutatedPos = 1;
      } else {
      	// Should never happen, as mutation = 0 if no mutable positions.
      	throw std::out_of_range("Algo5B: run out of mutable places!!??");
      }

      mutatedPos = Mutables[nextMutant][indexMutatedPos];
      newMutable = Mutables[nextMutant];
      newMutable.erase(newMutable.begin() + indexMutatedPos); 
      newMutable[0] = numMutablePos - 1;

      if(verbosity && (newMutable[0] == 0)) {
	std::cout << "\n A species run out of mutable genes \n";
      }

      // NOT using mutables
      // Uncomment if not using Mutables
      // numMutablePos = 0;
      // for(int i=0; i < numGenes; ++i) {
      // 	if(!Genotypes[nextMutant][i]) { 
      // 	  mutablePos[numMutablePos] = i;
      // 	  ++numMutablePos;
      // 	}
      // }
      // if(numMutablePos > 1) {
      // 	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      // } else if (numMutablePos == 1) {
      // 	mutatedPos = mutablePos[0];
      // } else {
      // 	// Should never happen, as mutation = 0 if no mutable positions.
      // 	throw std::out_of_range("Algo5: run out of mutable places!!??");
      // }

      // GMP
      mpz_ui_pow_ui(add_to_id, 2, mutatedPos);
      mpz_add(new_id, add_to_id, sp_id[nextMutant].get_mpz_t());

      
      // ************   5.6   ***************
      
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      
      // In lieu of the condition in the R code which is 
      // while((s <= num.species) && 
      //       !(identical(Genotype[s, ], newGenotype))) s <- s + 1
      // because if k >= numGenes in while it is because an 
      // existing genot. matched. 
      sp = 0;
      k = 0;

      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());
      
      //GMP
      for(sp = 0; sp < numSpecies; ++sp) {
      	if( mpz_cmp(new_id, sp_id[sp].get_mpz_t()) == 0)
      	  break;
      }

      // while ( (sp < numSpecies) && (k < numGenes) ) {
      // 	if (newGenotype[k] == Genotypes[sp][k]) k++;
      // 	else {
      // 	  sp++;
      // 	  k = 0;
      // 	}
      // }
      
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	// Mutables
	Mutables.push_back(newMutable);

	// GMP
	sp_id.push_back(mpz_class(new_id));

	tmpParam.popSize = 1;
	tmpParam.birth = fitness_CBN_std(mutatedPos,
					 restrictTable,
					 popParams[nextMutant].birth,
					 typeCBN,
					 newGenotype,
					 birthRate,
					 s,
					 numDrivers);
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << " at creation";
	  std::cout << "\n sp_id = " << sp_id[sp];
	  std::cout << "\n birthRate = " << birthRate;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	}


#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
#endif
	// number of mutates pos = those of the parent + 1!
	// Use assert for now, and later uncomment this
	// and delete accumulate expression underneath
	// tmpParam.mutation = mu * (numGenes - numMutablePos + 1); THIS EXPRESSION IS WRONG!!!
	// FIXME-mild
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	// Mutables 
	tmpParam.mutation = mu * newMutable[0];

	// tmpParam.mutation = mu * (numGenes - 
	//  			  std::accumulate(newGenotype.begin(),
	//  					  newGenotype.end(), 0));
	
	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(death, tmpParam.birth, tmpParam.mutation);
	}
#ifdef DEBUGW
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors
	// tmpParam.timeLastUpdate = currentTime;//not explicitly said in paper
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
		    << " from species "  <<   mutatedPos 
		    << ". New popSize = " << popParams[sp].popSize;
	}

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  death,
		  popParams[sp].birth);
	} else {
	  if(verbosity >= 2) {
	    std::cout << "\n             Mutation " << 
	      "to an already extinct species, sp " <<
	      sp <<"\n";
	  }
	  popParams[sp].popSize = 1.0;
	}


#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
	// popParams[sp].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      popParams[nextMutant].Flag = true;
      // pop of receiving mutant flagged above 
    } else { //       *********** We are sampling **********
      if(verbosity >= 2)
	std::cout <<"\n We are SAMPLING at time " << timeAllPopSample << "\n";
      
      currentTime = timeAllPopSample;
      
      //      for(int i = 0; i < numSpecies; i++) { // FIXME: make this change everywhere?
      for(int i = 0; i < popParams.size(); i++) {
	//	if(!popParams[i].Flag) {


	  if(popParams[i].popSize > 0.0) {
	    STOPASSERT(popParams[i].Flag == false);
	    STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	    STOPASSERT(timeAllPopSample - popParams[i].timeLastUpdate >= 0.0);


#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.9 ******\n " <<
	    "     Species  = " << i <<
	    "\n   pre-update popSize = " << popParams[i].popSize << 
	    "\n   timeAllPopSample = " << timeAllPopSample <<
	    "; popParams[i].timeLastUpdate = " << popParams[i].timeLastUpdate <<
	    "; t for Algo2 = " << timeAllPopSample - popParams[i].timeLastUpdate <<"\n";
	  
	  std::cout << " \n     species R " << popParams[i].R;
	  std::cout << " \n     species W " << popParams[i].W;
	  std::cout << " \n     species death " << death;
	  std::cout << " \n     species birth " << popParams[i].birth;
	  std::cout << " \n     species nextMutationTIme " << popParams[i].nextMutationTime;
#endif


	  popParams[i].popSize = 
	    Algo2(popParams[i].popSize,
		  timeAllPopSample - popParams[i].timeLastUpdate,
		  popParams[i].R,
		  popParams[i].W,
		  death,
		  popParams[i].birth);
	  popParams[i].Flag = true;

#ifdef DEBUGV
	  std::cout << "\n\n   post-update popSize = " << 
	    popParams[i].popSize << "\n";
#endif

	  
#ifdef DEBUGW	  
	  popParams[i].timeLastUpdate = -999999.99999;
	  // popParams[i].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	  }
	  //	}
      }

      timeAllPopSample += sampleEvery;

      type_resize = 0;

      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      // add one more, so as not check for resize if extinction
      if( outNS.n_cols <= (outNS_i + 2) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      
#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	
      for(int i = 0; i < popParams.size(); ++i) {
	outNS(i + 2, outNS_i) = popParams[i].popSize;
	totPopSize += popParams[i].popSize;
	// if(sps) {
	//   std::cout << "       Pop size pop " << i << " = " << popParams[i].popSize << "\n";
	// }

#ifdef DEBUGV
	std::cout << "       Pop size pop " << i << " = " << popParams[i].popSize << "\n";
#endif
      }


      // if(sps) {
      // 	std::cout << "\n  totPopSize = " << totPopSize << "\n";
      // }

#ifdef DEBUGV
      std::cout << "\n  totPopSize = " << totPopSize << "\n";
      if(Extinction)
	std::cout << "\n  But population became extinct as all ti are inf \n";

#endif

     


      if(Extinction) {
	// to make pretty output
	outNS_i++;
	outNS(0, outNS_i) = static_cast<double>(iter) + 1;
	outNS(1, outNS_i) = std::numeric_limits<double>::infinity();
	totPopSize = 0.0;
	for(int i = 0; i < numSpecies; ++i) {
	  outNS(i + 2, outNS_i) = 0.0;
	}
      }

      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || Extinction) simulsDone = true;
      
      forceSample = false;
    }
  }


  //timer.step("all big loop");

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].tis;
    returnParams(i, 7) = popParams[i].nextMutationTime;
    returnParams(i, 8) = popParams[i].mutation;
  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }

  //timer.step("finishing stuff");

  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("outi") = outNS_i + 1,
		      Named("TisInf") = Extinction);
  
  END_RCPP
}



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
		SEXP initSize_species_,
		SEXP initSize_iter_,
		SEXP seed_gsl_,
		SEXP verbose_) {
  BEGIN_RCPP
  using namespace Rcpp;

  std::cout << "\n WARNING: you have a good reason for using this version of Algo5, right?\n";

  //test_gmp();

    
  const IntegerMatrix restrictTable(restrictTable_);
  const int numDrivers = as<int>(numDrivers_);
  const int numGenes = as<int>(numGenes_);
  const std::string typeCBN = as<std::string>(typeCBN_);
  const double birthRate = as<double>(birthRate_);
  const double s = as<double>(s_);
  const double death = as<double>(death_);
  const double mu = as<double>(mu_);
  const double initSize = as<double>(initSize_);
  const double sampleEvery = as<double>(sampleEvery_);
  const double detectionSize = as<double>(detectionSize_);
  const int initSp = as<int>(initSize_species_);
  const int initIt = as<int>(initSize_iter_);
  const int seed = as<int>(seed_gsl_);
  const int verbosity = as<int>(verbose_);
  const double minNonZeroMut = mu * 0.01;  // to avoid == 0.0 comparisons
  const double ratioForce = 0.2; // If a single species this times
  // detectionSize, force a sampling to prevent going too far.

  std::vector<int>mutablePos(numGenes);

  bool Extinction = false;
  bool forceSample = false;
  bool simulsDone = false;

  double totPopSize = 0;
  double minNextMutationTime;
  double mutantTimeSinceLastUpdate;
  double currentTime = 0.0;
  double timeAllPopSample;

  int nextMutant;
  int numSpecies = 0;
  int iter = 0;
  int numMutablePos = 0;
  int mutatedPos = 0;
  int outNS_i = 0;
  int sp = 0;
  int k = 0;
  int type_resize = 0;

  //GSL rng
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, (unsigned long) seed);


  // The out.ns in R code; just an intermediate holder of popSizes
  // The first row is iter, the second is time
  // (in column major)
  // This should stay as arma, because of simple resizing of
  // both rows and cols (which initializes to 0)
  arma::mat outNS(initSp + 2, initIt);
  outNS.zeros();

  std::vector<myT> newGenotype(numGenes);
  //  std::vector<std::vector<int> > Genotypes(initSp, std::vector<int>(numGenes));
  std::vector<std::vector<myT> > Genotypes(1, std::vector<myT>(numGenes)); //int?
  Genotypes.reserve(initSp);

  spParams tmpParam; 
  std::vector<spParams> popParams(1);
  popParams.reserve(initSp);

  // I crucially assume that new mutations are placed in the next empty slot.
  // And if a species becomes extinct, its hole is left in all the data structures.
  // This could lead to large memory usage, but is the way to see species
  // disappearances.

  // 5.1 Initialize 
  numSpecies = 1;

  popParams[0].Flag = true;
  popParams[0].birth = birthRate;
  popParams[0].mutation = mu * numGenes;
  popParams[0].popSize = initSize;
  popParams[0].W = W_f(death, popParams[0].birth, popParams[0].mutation);
  popParams[0].R = R_f(death, popParams[0].birth, popParams[0].mutation);

  timeAllPopSample = currentTime + sampleEvery;
  
  outNS(0, 0) = 0.0;
  outNS(1, 0) = 0.0;
  outNS(2, 0) = initSize;
  


  while(!simulsDone) {
    iter++;

    //  ************   5.2   ***************

    if(verbosity >= 2)
      std::cout <<"\n\n\n***** Looping again through 5.2 \n";

    
    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].Flag) {
	// std::cout << " \n\n species i " << i;
	// std::cout << " \n species R " << popParams[i].R;
	// std::cout << " \n species W " << popParams[i].W;
	// std::cout << " \n species death " << death;
	// std::cout << " \n species birth " << popParams[i].birth;
	// std::cout << " \n species popSize " << popParams[i].popSize;
	popParams[i].tis = ti_f(popParams[i].R,
				popParams[i].W,
				death,
				popParams[i].birth,
				popParams[i].popSize);
	// std::cout << " \n species tis " << popParams[i].tis;
	popParams[i].nextMutationTime = popParams[i].tis + currentTime;
	// std::cout << " \n species nextMutTIme " << popParams[i].nextMutationTime;
	popParams[i].Flag = false;
	popParams[i].timeLastUpdate = currentTime;
// #ifdef DEBUGW
// 	std::cout <<" in 5.2, species i " << i << "  currentTime = " << 
// 	  currentTime <<"\n";
// #endif
      }
    }

    #ifdef DEBUGW
    for(int i = 0; i < popParams.size(); i++) {
      if(popParams[i].tis < 0.0) throw std::range_error("Algo5: tis < 0");
      else if(popParams[i].tis == 0.0) 
	std::cout << "Algo 5: BEWARE some tis == 0.0" << std::endl;
    }
    #endif
    
    // ******************** 5.3 and do we sample? *********** 

    // We need to find minimum to know if we need to sample the whole pop
    

    // Could some of this be avoided by keeping the index
    // of min nextMutationTime?
    // YES! 
    // add to previous loop and a first value.

    minNextMutationTime = std::numeric_limits<double>::infinity();
    nextMutant = -99; 

    // DELETE
    // std::cout << " \n numSpecies = " << numSpecies << " numGenes = " << numGenes;
    // std::cout << " \n sizeGenotypes " << Genotypes.size();

    for(int i = 0; i < popParams.size(); i++) {
      // std::cout << " \n\n species = " << i;
      if(popParams[i].popSize > 0.0) {
	// std::cout << " \n    and it exists with nextmtime " << popParams[i].nextMutationTime;
	if(popParams[i].nextMutationTime < minNextMutationTime) {
	  nextMutant = i;
	  minNextMutationTime = popParams[i].nextMutationTime;
	}     
      }
    }
    if(nextMutant < -90) {
      if(verbosity)
	std::cout << "Algo5: population became extinct";
      Extinction = true;
      // Note: the algorithm says they become extinct before a mutation.
      // (as per eq. 12 in p. 1232 of paper) But we will sample once, and
      // then set final sizes (after sample) to 0.  However, a pop. with
      // no holes for mutation left would have a ti of infinity (see
      // Algo2) but could keep growing indefinitely. Signal if any pop has
      // not mut. holes left FIXME!!
    }

    if(verbosity >= 2) {
      std::cout << "\n\n  iteration " << iter << "; minNextMutationTime = "
		<< minNextMutationTime 
		<< "; timeAllPopSample = " << timeAllPopSample << "\n";
    }

    // Do we need to sample the population?
    if((timeAllPopSample >= minNextMutationTime) && !forceSample) {// We are not sampling
      // ************   5.3   **************

      currentTime = minNextMutationTime;


      // ************   5.4   ***************

      STOPASSERT(popParams[nextMutant].timeLastUpdate >= 0.0);
      
      mutantTimeSinceLastUpdate = currentTime - 
	popParams[nextMutant].timeLastUpdate;
      
      ASSERT(mutantTimeSinceLastUpdate == (
					   popParams[nextMutant].nextMutationTime - 
					   popParams[nextMutant].timeLastUpdate));
      
      
      if(mutantTimeSinceLastUpdate > sampleEvery) {
	throw std::range_error("Algo5: mutantTimeSinceLastUpdate > sampleEvery");
      }
      
      popParams[nextMutant].popSize = Algo3(popParams[nextMutant].popSize,
					    mutantTimeSinceLastUpdate,
					    popParams[nextMutant].R,
					    popParams[nextMutant].W,
					    death,
					    popParams[nextMutant].birth);
      // Next line not in paper. It will be updated
      // in Algo2, as flagged.
      // popParams[nextMutant].timeLastUpdate = currentTime;
      
      if(popParams[nextMutant].popSize > (ratioForce * detectionSize)) {
	// What should we do here???
	// Allow bailing out even if tis become vanishingly small
	forceSample = true;
	if(verbosity) {
	  std::cout << "\n Forced sampling: \n    " << 
	    " popParams[nextMutant].popSize = " << 
	    popParams[nextMutant].popSize << " > ratioForce * detectionSize \n";
	  std::cout << " when nextMutant = " << nextMutant << "\n";
	}
      }      
      
      // ************   5.5   ***************
      // Note: impossible to have a second recorded mutation in
      // the same gene.  
      
      // This (the numMutablePos and the vector of mutablePos)
      // could be kept in the struct
      // Or in a different struct, since just accessed here
      // If we do that, below, when same species just recheck the fields
      // Actually, it is simpler; two species are the same if same mutable pos!
      numMutablePos = 0;
      for(int i=0; i < numGenes; ++i) {
      	if(!Genotypes[nextMutant][i]) { 
	  mutablePos[numMutablePos] = i;
      	  ++numMutablePos;
      	}
      }
      if(numMutablePos > 1) {
	mutatedPos = mutablePos[gsl_rng_uniform_int(r, numMutablePos)];
      } else if (numMutablePos == 1) {
	mutatedPos = mutablePos[0];
      } else {
	// Should never happen, as mutation = 0 if no mutable positions.
	throw std::out_of_range("Algo5: run out of mutable places!!??");
      }
      
      
      // ************   5.6   ***************
      
      newGenotype = Genotypes[nextMutant];
      newGenotype[mutatedPos] = 1;
      
      
      // In lieu of the condition in the R code which is 
      // while((s <= num.species) && 
      //       !(identical(Genotype[s, ], newGenotype))) s <- s + 1
      // because if k >= numGenes in while it is because an 
      // existing genot. matched. 
      sp = 0;
      k = 0;

      STOPASSERT(numSpecies == Genotypes.size());
      STOPASSERT(numSpecies == popParams.size());

      while ( (sp < numSpecies) && (k < numGenes) ) {
	if (newGenotype[k] == Genotypes[sp][k]) k++;
	else {
	  sp++;
	  k = 0;
	}
      }
      
      if(sp == numSpecies) {// New species
	++numSpecies;
	if(verbosity >= 2) {
	  std::cout <<"\n     Creating new species   " << (numSpecies - 1)
		    << "         from species "  <<   nextMutant;
	}
	
	Genotypes.push_back(newGenotype);
	tmpParam.popSize = 1;
	tmpParam.birth = fitness_CBN_std(mutatedPos,
					 restrictTable,
					 popParams[nextMutant].birth,
					 typeCBN,
					 newGenotype,
					 birthRate,
					 s,
					 numDrivers);
	
	if(verbosity >= 2) {
	  std::cout << " \n\n\n Looking at NEW species " << sp << "at creation\n";
	  std::cout << "\n birthRate = " << birthRate;
	  std::cout << "\n s = " << s;
	  std::cout << "\n parent fitness = " << popParams[nextMutant].birth;
	}


#ifdef DEBUGW
	int tmpSumMutPos =  std::accumulate(newGenotype.begin(),
					    newGenotype.end(),
					    0);
	if ((tmpSumMutPos < 0)  || (tmpSumMutPos > numGenes))
	  throw std::out_of_range("tmpSumMutPos out of range");
	if(! (tmpSumMutPos == (numGenes - numMutablePos + 1))) 
	  throw std::out_of_range("tmpSumMutPos != numMutPos expression");
	// tmpParam.mutation = mu * (numGenes - tmpSumMutPos);
#endif
	// number of mutates pos = those of the parent + 1!
	// Use assert for now, and later uncomment this
	// and delete accumulate expression underneath
	// tmpParam.mutation = mu * (numGenes - numMutablePos + 1);
	// FIXME-mild
	
	
	// Note: impossible to have a second recorded mutation in
	// the same gene.  See also 5.5
	
	tmpParam.mutation = mu * (numGenes - 
	 			  std::accumulate(newGenotype.begin(),
	 					  newGenotype.end(), 0));
	
	// Smallest mutations can be either 0 or mu * 1. Avoid equality comp.
	if(tmpParam.mutation < minNonZeroMut) {
	  // no mutation condition, in ti_f
	  tmpParam.W = -99.99;
	  tmpParam.R = -99.99; //unneeded, but do not set to 0
	} else{
	  tmpParam.W = W_f(death, tmpParam.birth, tmpParam.mutation);
	  tmpParam.R = R_f(death, tmpParam.birth, tmpParam.mutation);
	}
#ifdef DEBUGW
	tmpParam.timeLastUpdate = -99999.99999; // to catch errors
	// tmpParam.timeLastUpdate = currentTime;//not explicitly said in paper
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
		    << " from species "  <<   mutatedPos 
		    << ". New popSize = " << popParams[sp].popSize;
	}

	// Even if species became extint: Algo2 will return 0.
	if(popParams[sp].popSize > 0.0) {
	  popParams[sp].popSize = 1.0 + 
	    Algo2(popParams[sp].popSize,
		  currentTime - popParams[sp].timeLastUpdate,
		  popParams[sp].R,
		  popParams[sp].W,
		  death,
		  popParams[sp].birth);
	} else {
	  if(verbosity >= 2) {
	    std::cout << "\n             Mutation " << 
	      "to an already extinct species, sp " <<
	      sp <<"\n";
	  }
	  popParams[sp].popSize = 1.0;
	}


#ifdef DEBUGW
	popParams[sp].timeLastUpdate = -99999.99999; // to catch errors
	// popParams[sp].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	popParams[sp].Flag = true;
      }
      STOPASSERT( sp != nextMutant);
      
      //   ***************  5.7 ***************
      popParams[nextMutant].Flag = true;
      // pop of receiving mutant flagged above 
    } else { //       *********** We are sampling **********
      if(verbosity >= 2)
	std::cout <<"\n We are SAMPLING at time " << timeAllPopSample << "\n";
      
      currentTime = timeAllPopSample;
      
      //      for(int i = 0; i < numSpecies; i++) { // FIXME: make this change everywhere?
      for(int i = 0; i < popParams.size(); i++) {

#ifdef DEBUGV
	  std::cout << "\n\n     ********* 5.7 ******\n " <<
	    "     Species  = " << i <<
	    "\n   pre-update popSize = " << popParams[i].popSize << 
	    "\n   timeAllPopSample = " << timeAllPopSample <<
	    "; popParams[i].timeLastUpdate = " << popParams[i].timeLastUpdate <<
	    "; t for Algo2 = " << timeAllPopSample - popParams[i].timeLastUpdate <<"\n";
	  
	  std::cout << " \n     species R " << popParams[i].R;
	  std::cout << " \n     species W " << popParams[i].W;
	  std::cout << " \n     species death " << death;
	  std::cout << " \n     species birth " << popParams[i].birth;
	  std::cout << " \n     species s " << popParams[i].R;
#endif

	if(!popParams[i].Flag) {

	  STOPASSERT(popParams[i].timeLastUpdate >= 0.0);
	  STOPASSERT(timeAllPopSample - popParams[i].timeLastUpdate >= 0.0);


	  if(popParams[i].popSize > 0.0) {
	  popParams[i].popSize = 
	    Algo2(popParams[i].popSize,
		  timeAllPopSample - popParams[i].timeLastUpdate,
		  popParams[i].R,
		  popParams[i].W,
		  death,
		  popParams[i].birth);
	  popParams[i].Flag = true;
	  
#ifdef DEBUGW	  
	  popParams[i].timeLastUpdate = -999999.99999;
	  // popParams[i].timeLastUpdate = currentTime; //not explicitly said in paper
#endif
	  }
	}
#ifdef DEBUGV
	  std::cout << "\n\n   post-update popSize = " << 
	    popParams[i].popSize << "\n";
#endif
      }

      timeAllPopSample += sampleEvery;

      type_resize = 0;

      // resize outNS if needed
      if( outNS.n_rows <= (numSpecies + 2) ) type_resize += 1;
      // add one more, so as not check for resize if extinction
      if( outNS.n_cols <= (outNS_i + 2) ) type_resize += 2;

      // std::cout << "\n type resize = " << type_resize << "\n";
      // std::cout << "\n Before resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n outNS_i = " << outNS_i;


      if(type_resize == 1) 
	outNS.resize(2 * numSpecies + 2, outNS.n_cols);
      else if (type_resize == 2)
	outNS.resize(outNS.n_rows, 2 * outNS_i + 2);
      else if (type_resize == 3)
	outNS.resize(2 * numSpecies + 2, 2 * outNS_i + 2);

      // std::cout << "\n After resize";
      // std::cout << "\n outNS.n_rows = " << outNS.n_rows;
      // std::cout << "\n outNS.n_cols = " << outNS.n_cols;
      // std::cout << "\n  popParams.size() = " << popParams.size();
      // std::cout << "\n  numSpecies = " << numSpecies;

      outNS_i++;
      outNS(0, outNS_i) = static_cast<double>(iter);
      outNS(1, outNS_i) = currentTime;
      
      totPopSize = 0.0;
      
#ifdef DEBUGV
	std::cout << "\n Filling up outNS \n";
#endif
	
      for(int i = 0; i < popParams.size(); ++i) {
	outNS(i + 2, outNS_i) = popParams[i].popSize;
	totPopSize += popParams[i].popSize;
#ifdef DEBUGV
	std::cout << "       Pop size pop " << i << " = " << popParams[i].popSize << "\n";
#endif
      }

#ifdef DEBUGV
      std::cout << "\n  totPopSize = " << totPopSize << "\n";
#endif

     


      if(Extinction) {
	// to make pretty output
	outNS_i++;
	outNS(0, outNS_i) = static_cast<double>(iter) + 1;
	outNS(1, outNS_i) = std::numeric_limits<double>::infinity();
	totPopSize = 0.0;
	for(int i = 0; i < numSpecies; ++i) {
	  outNS(i + 2, outNS_i) = 0.0;
	}
      }

      if( (totPopSize >= detectionSize) ||
	  (totPopSize <= 0.0) || Extinction) simulsDone = true;
    }
  }

  // sanity checks
  if(Genotypes.size() != numSpecies) {
    std::cout << "\n ERROR: Genotypes.size() != numSpecies \n";
  }
  if(popParams.size() != numSpecies) {
    std::cout << "\n ERROR: popParams.size() != numSpecies \n";
  }


  // Resize before return
  outNS.resize(numSpecies + 2, outNS_i + 1);

  IntegerMatrix returnGenotypes(numSpecies, numGenes);
  
  for(int i = 0; i < numSpecies; ++i) {
    for(int j = 0; j < numGenes; ++j) {
      returnGenotypes(i, j) = Genotypes[i][j];
      }
  }
  

  // But this I do not need. Just initially, for checking.
  // I only need popSize. 
  #ifdef DEBUGW
  NumericMatrix returnParams(numSpecies, 9);
  for(int i = 0; i < numSpecies; ++i) {
    returnParams(i, 0) = popParams[i].Flag;
    returnParams(i, 1) = popParams[i].birth;
    returnParams(i, 2) = popParams[i].popSize;
    returnParams(i, 3) = popParams[i].timeLastUpdate;
    returnParams(i, 4) = popParams[i].W;
    returnParams(i, 5) = popParams[i].R;
    returnParams(i, 6) = popParams[i].tis;
    returnParams(i, 7) = popParams[i].nextMutationTime;
    returnParams(i, 8) = popParams[i].mutation;
  }
  #else
  std::string returnParams = "NA. Set DEBUGW to return it";
  #endif
  
  NumericVector popSizes(numSpecies);
  for(int i = 0; i < popParams.size(); ++i) {
    popSizes[i] = popParams[i].popSize;
  }
  
  return List::create(Named("NumSpecies") = numSpecies,
		      Named("TotalPopSize") = totPopSize,
		      Named("outNS") = wrap(outNS),
		      Named("Genotypes") = returnGenotypes,
		      Named("PopSizes") = popSizes,
		      Named("FinalTime") = currentTime,
		      Named("Params") = returnParams,
		      Named("outi") = outNS_i + 1);
  
  END_RCPP
}



// If I have a bidimension array of genotypes, then what? 
// I want something like slicing.

double fitness_linear_0(const int genotype[], const double birthRate, 
		      const double s, const int numDrivers) {
  // Crucial: drivers are always first in genotype
  // We could use valarrays, etc.
  int totalMut = 0;
  for (int i = 0; i < numDrivers; i++) totalMut += genotype[i];
  return birthRate + s * static_cast<double>(totalMut);
}


inline double fitness_linear_Rcpp(const Rcpp::IntegerMatrix::Column genotype, 
			     const double& birthRate, 
			     const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype
  // We could use sum of a sliced valarray
  // Or we could use a plain loop

  int totalMut = std::accumulate(genotype.begin(), 
				 genotype.begin() + numDrivers,
				 0);
  return birthRate + s * static_cast<double>(totalMut);
}



// this uses Armadillo
inline double fitness_linear_Arma(const arma::subview_col<unsigned int> Genotype,
			     const double& birthRate, 
			     const double& s, const int& numDrivers) {
  // Crucial: drivers are always first in genotype
  // We could use sum of a sliced valarray
  // Or we could use a plain loop

  int totalMut = arma::accu(Genotype.rows(0, numDrivers - 1));
  return birthRate + s * static_cast<double>(totalMut);
}



// how is this working with access to elements of subview_col via []??
// because that is OK too.
double fitness_CBN_Arma(const int& mutatedPos, 
		   Rcpp::IntegerMatrix restrictTable,
		   const double& fitnessParent, 
		   const std::string typeCBN,
		   const arma::subview_col<unsigned int> Genotype,
		   const double& birthRate, 
		   const double& s, 
		   const int& numDrivers,
			const int& DEBUG = 1,
			const int& DEBUG3 = 1) {
  //		   const int& DEBUG2 = 0) {
  
  using namespace Rcpp ;

  // FIXME: check numDrivers == ncol(restrictTable)??
  // or obtain the numDrivers from there??
  // Yes, numDrivers is known before hand. Take if from there
  // in some upper level function that calls this one.


  // FIXME: add checks of no negative numbers in deps: thisRestrict?

  // In fact, I could simplify, by having the restrictTable only
  // for drivers that DO depend. Would make it faster and exclude
  // a check below. But less clear.

  
  // will later become an argument
  double fitnessNo = 0.0;

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

    if(DEBUG) {
      if(thisRestrict[0] != mutatedPos ) {
	std::cout << std::endl << "thisRestrict[0] = " << thisRestrict[0] 
		<< "; mutatedPos  = " << mutatedPos  << std::endl;
      throw std::range_error("FitnessCBN: thisRestrict[0] != mutatedPos ");	
      }
      if(Genotype[mutatedPos] != 1){
	std::cout << " mutatedPos = " << mutatedPos << std::endl;
	throw std::range_error("FitnessCBN: genotype(mutatedPos) != 1");
      }
    }
    

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

    if(outFitnessYes) {
      return fitness_linear_Arma(Genotype, birthRate, s, numDrivers);
    } else {
      return fitnessNo;
    }

  }
}




// Does not use Armadillo
// FIXME: ojo, el restrictTable not const
double fitness_CBN_Rcpp(const int& mutatedPos, 
		   Rcpp::IntegerMatrix restrictTable,
		   const double& fitnessParent, 
		   const std::string typeCBN,
		   const Rcpp::IntegerMatrix::Column genotype,
		   const double& birthRate, 
		   const double& s, 
		   const int& numDrivers,
		   const int& DEBUG = 1) {
  //		   const int& DEBUG2 = 0) {
  
  using namespace Rcpp ;

  // FIXME: check numDrivers == ncol(restrictTable)??
  // or obtain the numDrivers from there??
  // Yes, numDrivers is known before hand. Take if from there
  // in some upper level function that calls this one.


  // FIXME: add checks of no negative numbers in deps: thisRestrict?

  // In fact, I could simplify, by having the restrictTable only
  // for drivers that DO depend. Would make it faster and exclude
  // a check below. But less clear.

  
  // will later become an argument
  double fitnessNo = 0.0;

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

    if(DEBUG) {
      if(thisRestrict[0] != mutatedPos ) {
	std::cout << std::endl << "thisRestrict[0] = " << thisRestrict[0] 
		<< "; mutatedPos  = " << mutatedPos  << std::endl;
      throw std::range_error("FitnessCBN: thisRestrict[0] != mutatedPos ");	
      }
      if(genotype[mutatedPos] != 1){
	std::cout << " mutatedPos = " << mutatedPos << std::endl;
	throw std::range_error("FitnessCBN: genotype(mutatedPos) != 1");
      }
    }
    

    if(!numDependencies) {
      // if(DEBUG2) {
      // 	std::cout << " FitnessCBN: exiting at !numDependencies " << std::endl;
      // }
      outFitnessYes = 1;
    } else {
      //outFitnessYes = 0;
      for(int i = 2; i < (2 + numDependencies); i++) {
	sumPresent += genotype[ thisRestrict[i] ];
      }
      if(typeCBN == "Multiple") {
        if(sumPresent) outFitnessYes = 1;
      } else{ // if(typeCBN == "CBN")
        if(sumPresent == numDependencies)
          outFitnessYes = 1; // return(fitnessYes(genotype, birth.rate, s, num.drivers))
      }
    }

    if(outFitnessYes) {
      return fitness_linear_Rcpp(genotype, birthRate, s, numDrivers);
    } else {
      return fitnessNo;
    }

  }
}




SEXP wrap_fitness_CBN_Arma(SEXP mutatedPos_,
			   SEXP genotypes_,
			   SEXP genNum_,
			   SEXP restrictTable_,
			   SEXP numDrivers_,
			   SEXP birthRate_, 
			   SEXP s_, 
			   SEXP fitnessParent_, 
			   SEXP typeCBN_,
			   SEXP retval_) {
  BEGIN_RCPP
    using namespace Rcpp;
  
  int mutatedPos = as<int>(mutatedPos_);
  IntegerMatrix restrictTable(restrictTable_);
  double fitnessParent = as<double>(fitnessParent_);
  std::string typeCBN = as<std::string>(typeCBN_);
  
  //  IntegerMatrix genotypes(genotypes_);
  
  arma::umat Genotypes = as<arma::umat>(genotypes_);
  
  
  int genNum = as<int>(genNum_);
  // IntegerMatrix::Column genotype = genotypes(_, genNum);
  
  double birthRate = as<double>(birthRate_);
  double s = as<double>(s_);
  int numDrivers = as<int>(numDrivers_);
  
  double retval = as<double>(retval_);
  
  retval = fitness_CBN_Arma(mutatedPos, restrictTable, fitnessParent,
			    typeCBN, Genotypes.col(genNum),
			    birthRate, s, numDrivers);
  
  return (wrap(retval));
  END_RCPP
    }


SEXP wrap_fitness_CBN_std(SEXP mutatedPos_,
			  SEXP genotypes_,
			  SEXP genNum_,
			  SEXP restrictTable_,
			  SEXP numDrivers_,
			  SEXP birthRate_, 
			  SEXP s_, 
			  SEXP fitnessParent_, 
			  SEXP typeCBN_,
			  SEXP retval_) {
  BEGIN_RCPP
    
    using namespace Rcpp;
  
  int mutatedPos = as<int>(mutatedPos_);
  IntegerMatrix restrictTable(restrictTable_);
  double fitnessParent = as<double>(fitnessParent_);
  std::string typeCBN = as<std::string>(typeCBN_);
  
  //arma::umat Genotypes = as<arma::umat>(genotypes_);

  // transpose and place in an std container
  IntegerMatrix genotypes(genotypes_);
  int ngenes = genotypes.nrow();
  int nsubjects = genotypes.ncol();
  int genNum = as<int>(genNum_);

  std::vector<std::vector<myT> > Genotypes(nsubjects, std::vector<myT>(ngenes));

  for(int i = 0; i < nsubjects; ++i){
    for(int j = 0; j < ngenes; ++j) {
      Genotypes[i][j] = genotypes(j, i);
    }
  }
    

  // IntegerMatrix::Column genotype = genotypes(_, genNum);
  
  double birthRate = as<double>(birthRate_);
  double s = as<double>(s_);
  int numDrivers = as<int>(numDrivers_);

  double retval = as<double>(retval_);
  
  retval = fitness_CBN_std(mutatedPos, restrictTable, fitnessParent,
			   typeCBN, Genotypes[genNum],
			   birthRate, s, numDrivers);

  return (wrap(retval));
  END_RCPP
}







// This is not using Armadillo,  but Rcpp matrices for genotypes
SEXP wrap_fitness_CBN_Rcpp(SEXP mutatedPos_,
		      SEXP genotypes_,
		      SEXP genNum_,
		      SEXP restrictTable_,
		      SEXP numDrivers_,
		      SEXP birthRate_, 
		      SEXP s_, 
		      SEXP fitnessParent_, 
		      SEXP typeCBN_,
		      SEXP retval_) {
  BEGIN_RCPP
  using namespace Rcpp;
  
  int mutatedPos = as<int>(mutatedPos_);
  IntegerMatrix restrictTable(restrictTable_);
  double fitnessParent = as<double>(fitnessParent_);
  std::string typeCBN = as<std::string>(typeCBN_);

  IntegerMatrix genotypes(genotypes_);
  int genNum = as<int>(genNum_);
  IntegerMatrix::Column genotype = genotypes(_, genNum);
  
  double birthRate = as<double>(birthRate_);
  double s = as<double>(s_);
  int numDrivers = as<int>(numDrivers_);

  double retval = as<double>(retval_);
  
  retval = fitness_CBN_Rcpp(mutatedPos, restrictTable, fitnessParent,
		       typeCBN, genotype, birthRate, s, numDrivers);

  return (wrap(retval));
  END_RCPP
}

// This is not using Armadillo,  but Rcpp matrices for genotypes
SEXP wrap_fitness_linear_Rcpp(SEXP allGenotypes_, 
			      SEXP genNum_, //this is the genotype number
			      SEXP birthRate_, 
			      SEXP s_, SEXP numDrivers_, 
			      SEXP v6_) {
  // we pass all the genotypes, and specify for which one we want fitness

  // Always use the RCPP macros to catch exceptions
  BEGIN_RCPP
  using namespace Rcpp ;
  IntegerMatrix allGenotypes(allGenotypes_);
  int genNum = as<int>(genNum_);
  double birthRate = as<double>(birthRate_);
  double s = as<double>(s_);
  int numDrivers = as<int>(numDrivers_);
  double v6 = as<double>(v6_);

  v6 = fitness_linear_Rcpp(allGenotypes(_, genNum), birthRate, s, numDrivers);

  return (wrap(v6));
  END_RCPP
}






SEXP wrap_ti(SEXP v1_, SEXP v2_, 
	     SEXP v3_, SEXP v4_, SEXP v5_,
	     SEXP v6_) {
  // Always use the RCPP macros to catch exceptions
  BEGIN_RCPP
  using namespace Rcpp ;
  double v1 = as<double>(v1_);
  double v2 = as<double>(v2_);
  double v3 = as<double>(v3_);
  double v4 = as<double>(v4_);
  double v5 = as<double>(v5_);
  double v6 = as<double>(v6_);

  v6 = ti_f(v1, v2, v3, v4, v5);
  return (wrap(v6));
  END_RCPP
}



SEXP wrap_Algo2(SEXP v1_, SEXP v2_, 
		SEXP v3_, SEXP v4_, SEXP v5_,
		SEXP v6_, SEXP vd_, SEXP v7_) {
  // Always use the RCPP macros to catch exceptions
  BEGIN_RCPP
  using namespace Rcpp ;
  double v1 = as<double>(v1_);
  double v2 = as<double>(v2_);
  double v3 = as<double>(v3_);
  double v4 = as<double>(v4_);
  double v5 = as<double>(v5_);
  double v6 = as<double>(v6_);
  double v7 = as<double>(v7_);
  int debug = as<int>(vd_);

  v7 = Algo2(v1, v2, v3, v4, v5, v6, debug);
  return (wrap(v7));
  END_RCPP
}


SEXP wrap_Algo3(SEXP v1_, SEXP v2_, 
		 SEXP v3_, SEXP v4_, SEXP v5_,
		 SEXP v6_, SEXP vd_, SEXP v7_) {
  // Always use the RCPP macros to catch exceptions
  BEGIN_RCPP
  using namespace Rcpp ;
  double v1 = as<double>(v1_);
  double v2 = as<double>(v2_);
  double v3 = as<double>(v3_);
  double v4 = as<double>(v4_);
  double v5 = as<double>(v5_);
  double v6 = as<double>(v6_);
  double v7 = as<double>(v7_);
  int debug = as<int>(vd_);

  v7 = Algo3(v1, v2, v3, v4, v5, v6, debug);
  return (wrap(v7));
  END_RCPP
}




// // An attempt to reduce size as spcies get popSize == 0.0. 
// // However, I doubt it will significantly impact speed, as
// // we would only reduce when sampling, which is a fraction of
// // number of iterations.
// // And in iterations, if popSize == 0.0 we do nothing, just loop
// // over and access an element.
// // In contrast, we need a cumbersome way of mapping things to outNS.


// void shrinkSinceZero(const int& zeroPos,
// 		  std::vector<spParamsE>& popParams,
// 		  std::vector<std::vector<myT> >& Genotypes,
// 		  std::vector<std::vector<myT> >& Mutables,		  
// 		  std::vector<mpz_class>& sp_id,
// 		  ) {
//   // We do not need to move all, just one, and shrink.
//   // Should be a lot faster than .erase(zeroPos)
  
//   lastSp = popParams.size() - 1;
//   if(zeroPos < lastSp) {
//     popParams[zeroPos] = popParams[lastSp];
//     Genotypes[zeroPos] = Genotypes[lastSp];
//     Mutables[zeroPos] = Mutables[lastSp];
//     sp_id[zeroPos] = sp_id[lastSp];
//   }
//   popParams.pop_back();
//   Genotypes.pop_back();
//   Mutables.pop_back();
//   sp_id.pop_back();
// }






void test_gmp() {

  std::cout << "\n *******  CHECKING  GMP  ********* \n";
  mpz_class a;
  a = 1234;

  std::vector<mpz_class> vm(1);

  vm.push_back(a);

  // mpz_t b; mpz_init (b);
  // mpz_pow_ui(b, a.get_mpz_t(), 2);
  
  // mpz_class c;
  // mpz_pow_ui(c.get_mpz_t(), a.get_mpz_t(), 2);

  mpz_t to_add;  mpz_init (to_add);
  mpz_ui_pow_ui(to_add, 2, 100); // instead of 100, the mutated pos
  
  mpz_t new_id; mpz_init (new_id);
  mpz_add(new_id, to_add, a.get_mpz_t()); // a would be the id of parent

  std::cout << "\n new id vs a " << mpz_cmp(new_id, a.get_mpz_t()) << "\n";

  mpz_class tmp1(to_add);
  vm.push_back(tmp1);

  vm.push_back(mpz_class(to_add));
    

  //mpz_class b;
  // mpz_t b;
  // mpz_init (b);
  // mpz_ui_pow_ui(b, 2, 100);
  

  // //  std::cout << "\n  b  = " << b << "\n";

  // mpz_t c;
  // mpz_init (c);
  // mpz_ui_pow_ui(c, 2, 100);
  // mpz_t d;
  // mpz_init (d);
  // mpz_ui_pow_ui(d, 2, 0);

  // mpz_add(c, b, d);

  // std::cout << " \n c > b " << mpz_cmp(c, b) << "\n"; 


  std::vector<mpz_class> myll(5);
  //std::vector<mpz_t> myl2(5);
  myll.reserve(150);

  for(int i = 0; i < 4; ++i) myll[i] = i + 10;
  myll[4] = 1236;
  
  for(int i = 0; i < 5; ++i) {
    std::cout << "\n myll[i] for i " << i << " =  "<< myll[i];
    std::cout << " comp [i] " << i << " " << cmp(myll[i], a);
    //std::cout << " comp2 [i] " << i << " " << cmp(myll[i], d);
  } 
  
  std::cout << "\n int " << std::numeric_limits<int>::max();
  std::cout << "\n long int " << std::numeric_limits<long int>::max();
  std::cout << "\n long long int " << std::numeric_limits<long long int>::max();
  std::cout << "\n uns long int " << std::numeric_limits<unsigned long long int>::max();
  long double pp = pow(2, 100);
  unsigned long long int mi = std::numeric_limits<unsigned long long int>::max();
  long double mm = static_cast<long double>(mi);
  bool the_comp = (mm > pp);
  bool the_comp2 = (mi > pp);
  std::cout << "\n larger than 2^100? mm " << the_comp <<
    ";  mi " << the_comp2;

  std::cout << "\n *******  DONE CHECKING GMP  ********* \n";

}






// Old stuff. When passing a column, I wanted to know if we were copying or not.
// So references are copied, but not content of columns. Look at this code, and uncomment
// includes for output_any and typeinfo, etc.


// double fitness_linear_verbose(Rcpp::IntegerMatrix::Column genotype, 
// 			      const double birthRate, 
// 			      const double s, const int numDrivers) {
//   // Crucial: drivers are always first in genotype
//   // We could use sum of a sliced valarray
//   // Or we could use a plain loop

//   // Verify we are only passing a reference
//   std::cout << "     Inside passed funct: genotype address " << &genotype << std::endl;
//   std::cout << "     Inside passed funct: what genotype points to? " << 
//     output_any(genotype) << std::endl;



//   genotype = genotype * 3;
//   for(int i = 0; i < 3; i++) {
//     std::cout << "        Inside passed funct, pointer genot[i]: " << &genotype[i] << std::endl;
//   }
  
//   int totalMut = std::accumulate(genotype.begin(), 
// 				 genotype.begin() + numDrivers,
// 				 0);
//   return birthRate + s * static_cast<double>(totalMut);
// }


// SEXP wrap_fitness_linear_verbose(SEXP allGenotypes_, 
// 				 SEXP genNum_,
// 				 SEXP birthRate_, 
// 				 SEXP s_, SEXP numDrivers_, 
// 				 SEXP v6_) {
//   // Always use the RCPP macros to catch exceptions
//   BEGIN_RCPP
//   using namespace Rcpp ;
//   IntegerMatrix allGenotypes(allGenotypes_);
//   int genNum = as<int>(genNum_);
//   double birthRate = as<double>(birthRate_);
//   double s = as<double>(s_);
//   int numDrivers = as<int>(numDrivers_);
//   double v6 = as<double>(v6_);

//   // Verify we are only passing a reference
//   Rcpp::IntegerMatrix::Column tmp = allGenotypes(_, genNum);
//   Rcpp::IntegerMatrix::Column* content_allGen;
//   content_allGen = &tmp;


//   std::cout << " Outside passed funct: the matrix address " << &allGenotypes << std::endl;
//   std::cout << " Outside passed funct: the column address " << &tmp << std::endl;

//   std::cout << " Outside passed funct: the column address, pos 1 " << &tmp[0] << std::endl;
  
//   tmp[0] = 55;
//   std::cout << " Outside passed funct: the column first pos value " << tmp[0] << std::endl;

//   std::cout << " Outside passed funct: the column first pos value, via all " << 
//     allGenotypes(0, genNum) << std::endl;

//   allGenotypes(0, 0) = 23;
//   std::cout << " Outside passed funct: allGenotypes(0,0) " << allGenotypes(0, 0) << std::endl;

  

//   std::string tipo = typeid(tmp).name();
//   std::cout << std::endl << " el tipo de la columna es " << tipo << std::endl;

//   // Something I am missing, as the next won't work  
//   // std::cout << " a lo que apunta? " << reinterpret_cast<int**>(tmp) <<std::endl;
//   // std::cout << " Outside passed funct: the column content? " << static_cast<Rcpp::MatrixColumn<13>*>(tmp) << std::endl;

//   // The next will not work. Oh well
//   // std::cout << " Outside passed funct: what matrix points to? " << 
//   //   output_any(allGenotypes) << std::endl;
//   std::cout << " Outside passed funct: what the column points to? " << 
//     output_any(tmp) << std::endl;
//   std::cout << " Outside passed funct: what the first pos contains? " << 
//     output_any(tmp[0]) << std::endl;

  
//   v6 = fitness_linear_verbose(allGenotypes(_, genNum), birthRate, s, numDrivers);

//   // Verify we are only passing a reference
//   for(int i = 0; i < 3; i++) {
//     std::cout << "   Outside passed funct, new value genot[i]: " << allGenotypes(i, genNum) << std::endl;
//     std::cout << "   Outside passed funct, pointer genot[i]: " << &allGenotypes(i, genNum) << std::endl;
//   }

//   return (wrap(v6));
//   END_RCPP
// }


// SEXP f1(SEXP inmat, SEXP dummy) {
//   BEGIN_RCPP
//   using namespace Rcpp ;
//   IntegerMatrix im(inmat);
//   int Dummy = as<int>(dummy);
//   std::cout << std::endl << " Matrix address " << &im << std::endl;
//   im(0, 0) = 99;
//   im(1, 1) = 33;
  
//   Dummy = 4;

//   return (wrap(Dummy));
//   END_RCPP
// }




// if using a matrix for popStats
  // likewise
  // rows are:
  //    birth, mutation, popSize, TimeLastUpdate, W, R, tis, NextMutationTime
  //    exists and flag are really boolean, but since using a matrix for this
  //    and for genotypes, I place them here.
  //    Alternatively, I could use a struct, and row major, but I'd use
  //    col major for genotypes and row major for stats. A mess.
  // int rowStats = 10;
  // NumericMatrix popParams(rowStats, maxSpecies);
  // // Initialize to catch errors
  // for(int j = 0; j < maxSpecies; j++) {
  //   popParams(0, j) = 0.0; //Exists
  //   popParams(1, j) = 0.0; //Flag
  //   for(int i = 2; i < rowStats; i++) 
  //     popParams(i, j) = -999999.99999;
  // }
  // popParams(0, 0) = 1.0;






// We have to loop over the struct for other reasons.
// Do not use this.
// popsizes is an element of a struct.
// inline bool areWeDone(const double popsizes, const double sizeForDetection) {
//   const double sumPopSizes = std::accumulate(popsizes.begin(), 
// 					     popsizes.end(),
// 					     0.0);
//   return (sumPopSizes >= sizeForDetection);
//   //   return true;
//   // else 
//   //   return false;
  
// }






// Some design decissions:


// Genotypes is from std. Unclear that arma gives me anything extra, and
// might be slower for just plain storage and access. 
// Moreover, of indiv. genotype as std object, a pain to add to arma matrix.
// Only con is having
// to convert to matrix at end and transpose (?)


// only popNS is arma object for easy resizing.

