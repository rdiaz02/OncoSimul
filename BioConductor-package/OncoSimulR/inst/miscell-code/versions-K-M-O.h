#ifndef _MatherR_H
#define _MatherR_H


// #include <Rcpp.h>
//#include <RcppArmadillo.h>
#include <RcppGSL.h> // is this here or in the .h?
//#include <Rcpp/Benchmark/Timer.h>



/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

/* RcppExport SEXP wrap_ti( SEXP v1_, SEXP v2_, 	 */
/* 			 SEXP v3_, SEXP v4_, SEXP v5_, */
/* 			 SEXP v6_); */


/* RcppExport SEXP wrap_Algo2( SEXP v1_, SEXP v2_, 	 */
/* 			    SEXP v3_, SEXP v4_, SEXP v5_, */
/* 			    SEXP v6_, SEXP vd_, SEXP v7_); */

/* RcppExport SEXP wrap_Algo3( SEXP v1_, SEXP v2_, 	 */
/* 			    SEXP v3_, SEXP v4_, SEXP v5_, */
/* 			    SEXP v6_, SEXP vd_, SEXP v7__); */

/* RcppExport SEXP wrap_fitness_linear_Rcpp(SEXP allGenotypes_,  */
/* 				    SEXP genNum_, */
/* 				    SEXP birthRate_,  */
/* 				    SEXP s_, SEXP numDrivers_,  */
/* 				    SEXP v6_); */


/* RcppExport SEXP wrap_fitness_CBN_Rcpp(SEXP mutatedPos_, */
/* 				 SEXP genotypes_, */
/* 				 SEXP genNum_, */
/* 				 SEXP restrictTable_, */
/* 				 SEXP numDrivers_, */
/* 				 SEXP birthRate_,  */
/* 				 SEXP s_,  */
/* 				 SEXP fitnessParent_,  */
/* 				 SEXP typeCBN_, */
/* 				 SEXP retval_); */

/* RcppExport SEXP wrap_fitness_CBN_Arma(SEXP mutatedPos_, */
/* 				 SEXP genotypes_, */
/* 				 SEXP genNum_, */
/* 				 SEXP restrictTable_, */
/* 				 SEXP numDrivers_, */
/* 				 SEXP birthRate_,  */
/* 				 SEXP s_,  */
/* 				 SEXP fitnessParent_,  */
/* 				 SEXP typeCBN_, */
/* 				 SEXP retval_); */

/* RcppExport SEXP wrap_fitness_CBN_std(SEXP mutatedPos_, */
/* 				 SEXP genotypes_, */
/* 				 SEXP genNum_, */
/* 				 SEXP restrictTable_, */
/* 				 SEXP numDrivers_, */
/* 				 SEXP birthRate_,  */
/* 				 SEXP s_,  */
/* 				 SEXP fitnessParent_,  */
/* 				 SEXP typeCBN_, */
/* 				 SEXP retval_); */


/* RcppExport SEXP Algorithm5(SEXP restrictTable_, */
/* 			   SEXP numDrivers_, */
/* 			   SEXP numGenes_, */
/* 			   SEXP typeCBN_, */
/* 			   SEXP birthRate_,  */
/* 			   SEXP s_,  */
/* 			   SEXP death_, */
/* 			   SEXP mu_, */
/* 			   SEXP initSize_, */
/* 			   SEXP sampleEvery_, */
/* 			   SEXP detectionSize_, */
/* 			   SEXP initSize_species_, */
/* 			   SEXP initSize_iter_, */
/* 			   SEXP seed_gsl_, */
/* 			   SEXP verbose_); */

/* RcppExport SEXP Algorithm5B(SEXP restrictTable_, */
/* 			   SEXP numDrivers_, */
/* 			   SEXP numGenes_, */
/* 			   SEXP typeCBN_, */
/* 			   SEXP birthRate_,  */
/* 			   SEXP s_,  */
/* 			   SEXP death_, */
/* 			   SEXP mu_, */
/* 			   SEXP initSize_, */
/* 			   SEXP sampleEvery_, */
/* 			   SEXP detectionSize_, */
/* 			   SEXP initSize_species_, */
/* 			   SEXP initSize_iter_, */
/* 			   SEXP seed_gsl_, */
/* 			   SEXP verbose_); */

/* RcppExport SEXP Algorithm5C(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_); */

/* RcppExport SEXP Algorithm5D(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

/* RcppExport SEXP Algorithm5E(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

/* RcppExport SEXP Algorithm5F(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

/* RcppExport SEXP Algorithm5G(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

/* RcppExport SEXP Algorithm5H(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */
/* RcppExport SEXP Algorithm5I(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

/* RcppExport SEXP Algorithm5J(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */

RcppExport SEXP Algorithm5K(SEXP restrictTable_,
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
			    SEXP typeFitness_);

/* RcppExport SEXP Algorithm5L(SEXP restrictTable_, */
/* 			    SEXP numDrivers_, */
/* 			    SEXP numGenes_, */
/* 			    SEXP typeCBN_, */
/* 			    SEXP birthRate_,  */
/* 			    SEXP s_,  */
/* 			    SEXP death_, */
/* 			    SEXP mu_, */
/* 			    SEXP initSize_, */
/* 			    SEXP sampleEvery_, */
/* 			    SEXP detectionSize_, */
/* 			    SEXP finalTime_, */
/* 			    SEXP initSize_species_, */
/* 			    SEXP initSize_iter_, */
/* 			    SEXP seed_gsl_, */
/* 			    SEXP verbose_, */
/* 			    SEXP speciesFS_, */
/* 			    SEXP ratioForce_, */
/* 			    SEXP typeFitness_); */


RcppExport SEXP Algorithm5M(SEXP restrictTable_,
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
			    SEXP typeFitness_);

RcppExport SEXP Algorithm5O(SEXP restrictTable_,
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
			    SEXP typeFitness_);


RcppExport SEXP Algorithm5P(SEXP restrictTable_,
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
			    SEXP typeFitness_);


/* RcppExport SEXP wrap_fitness_linear_verbose(SEXP allGenotypes_,  */
/* 					    SEXP genNum_, */
/* 					    SEXP birthRate_,  */
/* 					    SEXP s_, SEXP numDrivers_,  */
/* 					    SEXP v6_); */
/* RcppExport SEXP f1(SEXP inmat, SEXP dummy); */



// inline double W_f(const double& death, const double& growth, const double& mu);

#endif











