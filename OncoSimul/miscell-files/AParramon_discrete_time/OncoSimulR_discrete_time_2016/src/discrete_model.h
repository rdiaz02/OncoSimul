
#ifndef _DISCRETE_MODEL_H
#define _DISCRETE_MODEL_H

#include "new_restrict.h"
#include "bnb_common.h"


#include <Rcpp.h>

//Flags debug
//#define DEBUGMutation
//#define DEBUGEqGenotypes
//#define DEBUGdeathAndBirth
//#define DEBUGPrintDiscreteModel


struct Clon {
  Genotype genotype;
  double popSize;
  double absfitness; //absolute fitness
};


struct DiscreteModel {
  double popIni;
  double popCurrent;
  int tMax;
  int tPreset;
  std::vector<Clon> listClones;
  double avefitness; //average fitness
  std::vector<int> totalGenes;
  std::vector<double> mu;
  fitnessEffectsAll fE;
  std::map<int, std::string> intName;
  fitness_as_genes fg;
};


struct Historical {
  std::vector<Genotype> genot;
  std::vector<double> popSizes;
  std::vector<int> index;
  std::vector<double> time;
  std::vector<double> totPopSize;
  std::vector<double> popSizeLargestClone;
  std::vector<double> propPopSizeLargestClone;
  std::vector<int> maxNumDrivers;
  std::vector<int> numDriversLargestClone;
  int outNS_i;
  int maxram;
  PhylogName phylog;
};


// Simulates the growth and death of cells
int deathAndBirth(DiscreteModel &dm);

// Simulates the process of mutation of cells
DiscreteModel mutation(const DiscreteModel &dm, Historical &hist, std::mt19937 &ran_gen);

// Stores information about the evolutionary process in Historical struct
void historical(Historical &hist, const DiscreteModel &dm);

// Adding a cell with genotype g to the set of clones
void addCellToModel(DiscreteModel &dm2, const Genotype &g);

// Compares the two cells genotype
bool eqGenotypes(const Genotype &g1, const Genotype &g2);

// Calculates the fitness given a genotype
double calculateFitnessGenotype(const Genotype &g, const fitnessEffectsAll &fE);

// Calculates the size of the total population
int calculateSizePopulationModelo(const DiscreteModel &dm);

// Calculates the average fitness of the total population
double calculateFitnessAverageModelo(const DiscreteModel &dm);

// Prints information about the state of the simulation
void printDiscreteModel(const DiscreteModel &dm);

// Prepares the output of the simlation
Rcpp::List structDiscreteModel(const DiscreteModel &dm, const Historical &hist);


#endif