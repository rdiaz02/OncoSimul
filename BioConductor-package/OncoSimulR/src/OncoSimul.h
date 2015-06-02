#ifndef _OncoSimul_H
#define _OncoSimul_H

// #include <RcppGSL.h> 
#include <Rcpp.h> 


RcppExport SEXP Algorithm5(SEXP restrictTable_,
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
			   SEXP mutatorPhenotype_,
			   SEXP initMutant_,
			   SEXP maxWallTime_,
			   SEXP keepEvery_,
			   SEXP alpha_,
			   SEXP sh_,
			   SEXP K_,
			   SEXP endTimeEvery_,
			   SEXP finalDrivers_);


#endif

