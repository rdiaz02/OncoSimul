

#include "discrete_model.h"
#include "new_restrict.h"
#include <Rcpp.h>
#include <iomanip> 
#include <algorithm>
#include <random>
#include <chrono>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;
using namespace Rcpp;
  
// [[Rcpp::export]]
Rcpp::List discreteModel(Rcpp::List rFE,
	Rcpp::NumericVector mu_,
	double popIni,
	int tMax,
	double seed,
	double sampleEvery,
	double keepEvery, 
	int verbosity) {

	//initialization: population without mutations
	double popCurrent = popIni; // initial population
	int tPreset = 0; // initial time
	std::vector<Clon> listClones; 

	fitnessEffectsAll fE = convertFitnessEffects(rFE); // fitness especification

	std::vector<int> totalGenes = allGenesinFitness(fE); //total genes in genotype represented by int  

	// Vector with probabilities of mutation
	std::vector<double> mu;
	if(mu_.size()!=1){
		mu = Rcpp::as<std::vector<double> >(mu_);
	} else {
		std::vector<double> muAux = Rcpp::as<std::vector<double> >(mu_);
		mu = {};
		for(int i = 0; i < totalGenes.size(); i++){
			mu.push_back(muAux[0]);
		}
	}
	
	// It starts with a populaiton without mutations
	std::vector<int> orderEff = {};
	std::vector<int> epistRtEff = {};
	std::vector<int> rest = {};
	Genotype genotype = {orderEff, epistRtEff, rest};
	Clon clon = {genotype, popIni, calculateFitnessGenotype(genotype, fE)};
	listClones.push_back(clon);

	// For random operations
	unsigned int rseed = static_cast<unsigned int>(seed);
	if(seed == 0) {
	  rseed = std::random_device{}();
	}
	std::mt19937 ran_gen(rseed);

	// Initialization of DicreteModel
	std::map<int, std::string> intName = mapGenesIntToNames(fE);
	fitness_as_genes genesInFitness = fitnessAsGenes(fE);
	DiscreteModel dm = {popIni, popCurrent, tMax, tPreset, listClones, 0, totalGenes, mu, fE, intName, genesInFitness}; // 0 is average fitness
	dm.avefitness =  calculateFitnessAverageModelo(dm);

	//Initialization of Historical
	PhylogName phylog = { {}, {}, {} };
	Historical hist = { {}, {}, {}, {}, {}, {}, {}, {}, {}, -1, 2000, phylog};

	//Start (round zero)
	historical(hist, dm);
	dm.tPreset++;

	// ------- Graphic output --------
	if(verbosity >= 1){
    	Rcpp::Rcout << "\n Start --> Population size: " << dm.popCurrent << "---- average fitness: " << dm.avefitness << " ---- clones: " << dm.listClones.size();
	}

    if(verbosity >= 2){
      printDiscreteModel(dm);
    }
    // -------------------------------

	//Loop:
	for(dm.tPreset; dm.tPreset<=dm.tMax; dm.tPreset++){
		if(deathAndBirth(dm)==-1){
			return List::create(Named("finalTime") = dm.tPreset); //if population < 1 ---> end
		}
		dm = mutation(dm, hist, ran_gen);
		
		 // keep the evolution of simulation
		if(keepEvery == 0){
			if(dm.tPreset==dm.tMax){
				historical(hist, dm);
			}
		} else if ((((dm.tPreset)%(int)keepEvery)==0) || (dm.tPreset==dm.tMax)){
			historical(hist, dm);
		}

		// ------- Graphic output --------
		if(sampleEvery > 0 && verbosity >= 2){
			if(((dm.tPreset)%(int)sampleEvery)==0){
				Rcpp::Rcout << "\n After mutation " << dm.tPreset << "--> Population size: " << dm.popCurrent << "---- average fitness: " << dm.avefitness << " ---- clones: " << dm.listClones.size();
				if(verbosity >= 3){
					printDiscreteModel(dm);
				}
			}
		}
		// -------------------------------
	}

	// ------- Graphic output --------
	if(verbosity >= 1){
		Rcpp::Rcout << "\n End --> Population size: " << dm.popCurrent << "---- average fitness: " << dm.avefitness << " ---- clones: " << dm.listClones.size();
	}

	if(verbosity >= 2){
		printDiscreteModel(dm);
		Rcpp::Rcout <<"\n";
	}
	// -------------------------------

	return structDiscreteModel(dm, hist);	
}


int deathAndBirth(DiscreteModel &dm){
	#ifdef DEBUGdeathAndBirth
		Rcpp::Rcout << "\n\t Into deathAndBirth function: ";
	#endif
  
  	// Initialices a generator of random numbers
  	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	std::default_random_engine generator (seed);
  	
  	// Loop that simulates the growth and death of each clon
	for(int i=0; i < dm.listClones.size(); i++){

		// Chooses a random number proportional to the size of the clone and his fitness
		double lambda = (double)dm.popIni/(double)dm.popCurrent*(dm.listClones[i].absfitness/dm.avefitness)*dm.listClones[i].popSize;
		std::poisson_distribution<int> distribution((int)lambda);
		int number = distribution(generator);

		#ifdef DEBUGdeathAndBirth
			Rcpp::Rcout << "\n\t\t Clon: " << i << " ---- size: " << dm.listClones[i].popSize <<;
		#endif
		
		// Assigns the new population to clone
		dm.listClones[i].popSize = number;
		
		#ifdef DEBUGdeathAndBirth
			Rcpp::Rcout << "\n\t\t\t After variation: ---- size: " << dm.listClones[i].popSize << " ---- fitness: " << dm.listClones[i].absfitness << " ---- average fitness: " << dm.avefitness << " ---- lambda: " << lambda << " ---- number: " << number;
		#endif	
	}

	// Calculates and assigns the size of the total population
	dm.popCurrent = calculateSizePopulationModelo(dm);
	if(dm.popCurrent < 1){
				Rcpp::Rcout << "\n\n WARNING: Population = " << dm.popCurrent << " -->  extinct population";
				return -1; 
	}

	// Calculates and assigns the average fitness of the total population
	dm.avefitness = calculateFitnessAverageModelo(dm);
	return 0;
}


DiscreteModel mutation(const DiscreteModel &dm, Historical &hist, std::mt19937 &ran_gen){
	DiscreteModel dm2 = {dm.popIni, 0, dm.tMax, dm.tPreset, {}, 0, dm.totalGenes, dm.mu, dm.fE, dm.intName, dm.fg};

	// Initialices a generator of random numbers
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator (seed);
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	bool phylog_flag; //for keep phylog
	
	#ifdef DEBUGMutation
		Rcpp::Rcout << "\n\t Into Mutation function: ";
	#endif

	// Loop that simulates the process of mutation of each clon
	for(int i=0; i < dm.listClones.size(); i++){

		// Calculates the genes that are not mutated, so that genes can be mutated now
		std::vector<int> genesMut = allGenesinGenotype(dm.listClones[i].genotype); //mutated genes in that genotype
		std::sort (genesMut.begin(),genesMut.end());
		std::vector<int> candMut(dm.totalGenes.size()); //initalice a vector of size totalGenes 
		std::vector<int>::iterator itGenesMut;
		itGenesMut = std::set_difference (dm.totalGenes.begin(), dm.totalGenes.end(), genesMut.begin(), genesMut.end(), candMut.begin());
		candMut.resize(itGenesMut-candMut.begin()); // possible genes that can be mutated


		#ifdef DEBUGMutation
			Rcpp::Rcout << "\n\t\t Clon: " << i << " ---- size: " << dm.listClones[i].popSize;
			Rcpp::Rcout << "\n\t\t\t mutated genes : ";
			for(int w=0; w< genesMut.size(); w++){
		  		Rcpp::Rcout << genesMut[w] << " ";
		  	}
 
		  	Rcpp::Rcout << "\n\t\t\t non mutated genes: ";
			for(int w=0; w< candMut.size(); w++){
		  		Rcpp::Rcout << candMut[w] << " ";
		  	}
	  	#endif

		// Loop that simulates the process of mutation of each cell of the clon
		for(int j=0; j < dm.listClones[i].popSize; j++){
			phylog_flag = FALSE; // no mutations yet

			// Many random numbers as possible genes can be mutated are generated
		  	double number;
		  	std::vector<int> newGenesMut = {};
		  	for(int k=0; k < candMut.size(); k++){
				number = distribution(generator);
				if(number <= dm.mu[candMut[k]-1]){ // if random number is hicher than probabilities of mutation...
					#ifdef DEBUGMutation
						Rcpp::Rcout << "\n\t\t\t It has been mutated gene: " << candMut[k];
					#endif
					phylog_flag = TRUE; // Has happened a mutation
					newGenesMut.push_back(candMut[k]);
				}
			}

			// Creates a genotype by the actual genotype and the mutated genes
			Genotype g = createNewGenotype(dm.listClones[i].genotype, newGenesMut, dm.fE, ran_gen, false);

			// Keeps the phylog
			if(phylog_flag){
				addToPhylog(hist.phylog, dm.listClones[i].genotype, g, dm.tPreset, dm.intName, dm.fg);
			} 

			// Adds cell to the set of clones
			addCellToModel(dm2, g);			
		}
	}

	// Actualization of size population and average fitness
	dm2.popCurrent = calculateSizePopulationModelo(dm2);
	dm2.avefitness = calculateFitnessAverageModelo(dm2);

	return dm2;
}


void historical(Historical &hist, const DiscreteModel &dm){

	hist.outNS_i++;
	hist.time.push_back(dm.tPreset);
	int aux;
	int tmp_PopSizeLargestClone = 0;
	int tmp_MaxNumDrivers = -1;
	int tmp_NumDriversLargestClone = -1;
	for(int i=0; i < dm.listClones.size(); i++){
		hist.genot.push_back(dm.listClones[i].genotype);
		hist.popSizes.push_back(dm.listClones[i].popSize);
		hist.index.push_back(hist.outNS_i);
		
		if(dm.listClones[i].popSize > tmp_PopSizeLargestClone){
			tmp_PopSizeLargestClone = dm.listClones[i].popSize;
			tmp_NumDriversLargestClone = getGenotypeDrivers(dm.listClones[i].genotype, dm.fE.drv).size();
		}

		aux = getGenotypeDrivers(dm.listClones[i].genotype, dm.fE.drv).size();
		if(aux > tmp_MaxNumDrivers){
			tmp_MaxNumDrivers = aux;
		}
	}
	hist.totPopSize.push_back(dm.popCurrent);
	hist.popSizeLargestClone.push_back(tmp_PopSizeLargestClone);
	hist.propPopSizeLargestClone.push_back(tmp_PopSizeLargestClone/dm.popCurrent);
	hist.maxNumDrivers.push_back(tmp_MaxNumDrivers);
	hist.numDriversLargestClone.push_back(tmp_NumDriversLargestClone);
}


void addCellToModel(DiscreteModel &dm2, const Genotype &g){

	// Seeks clone with cells with the same genotype
	for(int i=0; i < dm2.listClones.size(); i++){
		if(eqGenotypes(dm2.listClones[i].genotype, g)){
			dm2.listClones[i].popSize++;
			return;
		}
	}

	// Else creates a new clon with that cell
	Clon clon = {g, 1, calculateFitnessGenotype(g, dm2.fE)};
	dm2.listClones.push_back(clon);
	return;
}


bool eqGenotypes(const Genotype &g1, const Genotype &g2){
	
	#ifdef DEBUGEqGenotypes
		Rcpp::Rcout << "\n\t\t Into eqGenotypes function: ";
		Rcpp::Rcout << "\n\t\t\t Genotype 1: ";
		for(int i=0; i < g1.orderEff.size(); i++){
			Rcpp::Rcout << g1.orderEff[i] << " ";
		}
		Rcpp::Rcout << "-- ";
		for(int i=0; i < g1.epistRtEff.size(); i++){
			Rcpp::Rcout << g1.epistRtEff[i] << " ";
		}
		Rcpp::Rcout << "-- ";
		for(int i=0; i < g1.rest.size(); i++){
			Rcpp::Rcout << g1.rest[i] << " ";
		}

		Rcpp::Rcout << "\n\t\t\t Genotype 2: ";
		for(int i=0; i < g2.orderEff.size(); i++){
			Rcpp::Rcout << g2.orderEff[i] << " ";
		}
		Rcpp::Rcout << "-- ";
		for(int i=0; i < g2.epistRtEff.size(); i++){
			Rcpp::Rcout << g2.epistRtEff[i] << " ";
		}
		Rcpp::Rcout << "-- ";
		for(int i=0; i < g2.rest.size(); i++){
			Rcpp::Rcout << g2.rest[i] << " ";
		}
	#endif


	if ((g1.orderEff.size() == g2.orderEff.size()) && (g1.epistRtEff.size() == g2.epistRtEff.size()) && (g1.rest.size() == g2.rest.size())){
		if (std::equal(g1.orderEff.begin(), g1.orderEff.begin() + g1.orderEff.size(), g2.orderEff.begin() ) && 
			std::equal(g1.epistRtEff.begin(), g1.epistRtEff.begin() + g1.epistRtEff.size(), g2.epistRtEff.begin() ) && 
			std::equal(g1.rest.begin(), g1.rest.begin() + g1.rest.size(), g2.rest.begin() ) ){
			#ifdef DEBUGEqGenotypes
				Rcpp::Rcout << " --> Are equals "<< std::endl;
			#endif
			return true;
		}
	}

	#ifdef DEBUGEqGenotypes
		Rcpp::Rcout << " --> Are differents "<< std::endl;
	#endif

	return false;
}


double calculateFitnessGenotype(const Genotype &g, const fitnessEffectsAll &fE){

	if((g.orderEff.size()==0) && (g.epistRtEff.size()==0) && (g.rest.size()==0)){
		return 1;
	}
	return prodFitness(evalGenotypeFitness(g, fE));
}


int calculateSizePopulationModelo(const DiscreteModel &dm){
	int popSizeTotal = 0;

	for(int i=0; i < dm.listClones.size(); i++){
		popSizeTotal = popSizeTotal + dm.listClones[i].popSize;
	}

	return popSizeTotal;
}


double calculateFitnessAverageModelo(const DiscreteModel &dm){
	double cumf = 0;

	for(int i=0; i < dm.listClones.size(); i++){
		cumf = cumf + dm.listClones[i].popSize*dm.listClones[i].absfitness;
	}

	return cumf/(double)dm.popCurrent;
}


void printDiscreteModel(const DiscreteModel &dm){

	for(int i=0; i < dm.listClones.size(); i++){
		std::vector<Genotype> genot_out;
		genot_out.push_back(dm.listClones[i].genotype);
		std::vector<std::vector<int> > genot_out_v = genot_to_vectorg(genot_out);
		std::vector<std::string> genotypesAsStrings = genotypesToNameString(genot_out_v, dm.fg, dm.intName);

		Rcpp::Rcout << "\n\t Clon " << i << " ---- size: " << dm.listClones[i].popSize;
		Rcpp::Rcout << " ---- genotype: " << genotypesAsStrings[0];
		Rcpp::Rcout << " ---- fitness: " << dm.listClones[i].absfitness;
		
		#ifdef DEBUGPrintDiscreteModel
			for(int j=0; j < dm.listClones[i].genotype.orderEff.size(); j++){
				Rcpp::Rcout << dm.listClones[i].genotype.orderEff[j] << " ";
			}
			Rcpp::Rcout << "-- ";
			for(int j=0; j < dm.listClones[i].genotype.epistRtEff.size(); j++){
				Rcpp::Rcout << dm.listClones[i].genotype.epistRtEff[j] << " ";
			}
			Rcpp::Rcout << "-- ";
			for(int j=0; j < dm.listClones[i].genotype.rest.size(); j++){
				Rcpp::Rcout << dm.listClones[i].genotype.rest[j] << " ";
			}
		#endif
		
	}
}


Rcpp::List structDiscreteModel(const DiscreteModel &dm, const Historical &hist){
	
	// General information and pops.by.time matrix
	std::vector<std::vector<int> > genot_out_v = genot_to_vectorg(hist.genot);
	std::vector<std::vector<int> > uniqueGenotypes_vector_nr  = uniqueGenot_vector(genot_out_v);
	IntegerMatrix returnGenotypes = nr_create_returnGenotypes(dm.fE.genomeSize, uniqueGenotypes_vector_nr);
	Rcpp::NumericMatrix outNS = create_outNS(uniqueGenotypes_vector_nr,
											   genot_out_v,
											   hist.popSizes,
											   hist.index, hist.time,
											   hist.outNS_i, hist.maxram);
	std::vector<std::string> genotypesAsStrings = genotypesToNameString(uniqueGenotypes_vector_nr, dm.fg, dm.intName);


	// Drivers information
	int maxNumDrivers = 0;
	int totalPresentDrivers = 0;
	std::vector<int> countByDriver(dm.fE.drv.size(), 0);
	std::vector<int> presentDrivers;
	driverCounts(maxNumDrivers, totalPresentDrivers,
	       countByDriver, presentDrivers,
	       returnGenotypes, dm.fE.drv);
	std::string driversAsString = driversToNameString(presentDrivers, dm.intName);

	// Stats
	NumericMatrix perSampleStats(hist.outNS_i + 1, 5);
  	fill_SStats(perSampleStats, hist.totPopSize, hist.popSizeLargestClone,
  	      hist.propPopSizeLargestClone, hist.maxNumDrivers, hist.numDriversLargestClone);

	return
    List::create(Named("pops.by.time") = outNS,
		 Named("NumClones") = uniqueGenotypes_vector_nr.size(), 
		 Named("TotalPopSize") = dm.popCurrent,
		 Named("Genotypes") = returnGenotypes,
		 Named("GenotypesWDistinctOrderEff") = Rcpp::wrap(uniqueGenotypes_vector_nr),
		 Named("GenotypesLabels") = Rcpp::wrap(genotypesAsStrings),
		 Named("MaxNumDrivers") = maxNumDrivers,
		 Named("MaxDriversLast") = hist.maxNumDrivers[hist.outNS_i],
		 Named("NumDriversLargestPop") =  hist.numDriversLargestClone[hist.outNS_i],
		 Named("LargestClone") = hist.popSizeLargestClone[hist.outNS_i],
		 Named("PropLargestPopLast") = hist.propPopSizeLargestClone[hist.outNS_i],
		 Named("FinalTime") = dm.tMax,
		 Named("NumIter") = dm.tMax, //no tiene mucho sentido, es el numero maximo de iteraciones en el BNB algoritmo
		 Named("HittedWallTime") = false, 
		 Named("HittedMaxTries") = false,
		 Named("TotalPresentDrivers") = totalPresentDrivers,
		 Named("CountByDriver") = countByDriver,
		 Named("OccurringDrivers") = driversAsString,
		 Named("PerSampleStats") = perSampleStats,
		 Named("other") = List::create(Named("attemptsUsed") = 1,
					       Named("errorMF") = -99,
					       //Named("errorMF_size") = e1,
					       //Named("errorMF_n_0") = n_0,    
					       Named("minDMratio") = -99,
					       Named("minBMratio") = -99,
					       //Named("errorMF_n_1") = n_1,
					       Named("PhylogDF") =  DataFrame::create(
										      Named("parent") = hist.phylog.parent,
										      Named("child") = hist.phylog.child,
										      Named("time") = hist.phylog.time
										      ),
					       Named("UnrecoverExcept") = false) 
	);
}