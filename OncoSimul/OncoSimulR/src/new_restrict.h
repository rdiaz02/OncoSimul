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



#ifndef _NEW_RESTRICT_H__
#define _NEW_RESTRICT_H__

#include "debug_common.h"
#include "common_classes.h"
// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include <Rcpp.h>
#include <limits>
#include <random>

// Yes, even if covr suggests epistasis, Poset_struct and
// Gene_Module_struct are not used they are used a lot.
// There are many vectors of these structs. This is just a problem
// of coverage testing of structs. Google for it.

enum class Dependency {monotone, semimonotone, xmpn, single, NA};
// enum class TypeModel {exp, bozic1, mcfarlandlog, mcfarland,
//     beerenwinkel, mcfarland0,  bozic2};
// enum class TypeModel {exp, bozic1, mcfarlandlog};

struct genesWithoutInt {
  int shift; // access the s as s[index of mutation or index of mutated
	     // gene in genome - shift]. shift is the min. of NumID, given
	     // how that is numbered from R. We assume mutations always
	     // indexed 1 to something. Not 0 to something.
  // If shift is -9, no elements The next first two are not really
  // needed. Will remove later. Nope! we use them to provide nice output.
  std::vector<int> NumID;
  std::vector<std::string> names;
  std::vector<double> s;
};


struct fitnessLandscape_struct {
  std::vector<int> NumID;
  std::vector<std::string> names;
  // zz: maybe not a char; hold on
  std::map<std::string, double> flmap;
  std::map<std::string, std::string> flFDFmap; //New line to define flFDFmap
  std::map<std::string, std::string> flfVarsmap; //New line to define flfVarsmap
};

struct evalFVars_struct {//structure to store the map fVars to fitness (double)
  std::map<std::string, double> evalFVarsmap;
};

struct Poset_struct {
  Dependency typeDep;
  int childNumID; //Not redundant
  double s;
  double sh;
  std::vector<int> parentsNumID;
    // The next two are clearly redundant but a triple check
  std::string child;
  std::vector<std::string> parents;
};

// We use same structure for epistasis and order effects. With order
// effects, NumID is NOT sorted, but reflects the order of the
// restriction. And checking is done using that fact.

struct epistasis {
  double s;
  std::vector<int> NumID; //a set instead? nope.using includes with epistasis
  std::vector<std::string> names; // will remove later
};



struct Gene_Module_struct {
  std::string GeneName;
  std::string ModuleName;
  int GeneNumID;
  int ModuleNumID;
};

struct fitnessEffectsAll {
  bool gMOneToOne;
  int genomeSize;
  // We use allOrderG or allEpistRTG to place new mutations in their
  // correct place (orderEff or epistRtEff). Only one is needed.  Use the
  // one that is presumably always shorter which is allOrderG. And this is
  // sorted.
  std::vector<int> allOrderG; // Modules or genes if one-to-one.
  // std::vector<int> allEpistRTG;

  // This makes it faster to run evalPosetConstraints
  std::vector<int> allPosetG; //Modules or genes if one-to-one. Only
			      //poset. Not epist.
  std::vector<Poset_struct> Poset;
  std::vector<epistasis> Epistasis;
  std::vector<epistasis> orderE;
  // std::vector<Gene_Module_struct> Gene_Module_tabl;
  std::vector<Gene_Module_struct> Gene_Module_tabl;
  std::vector<int> allGenes; //used whenever a mutation created. Genes,
			     //not modules. Sorted.
  std::vector<int> drv; // Sorted.
  genesWithoutInt genesNoInt;
  // zz:
  std::vector<std::string> fVars; //New line to store fVars
  fitnessLandscape_struct fitnessLandscape;
};

inline fitnessEffectsAll nullFitnessEffects() {
  // Make it explicit
  fitnessEffectsAll f;
  f.gMOneToOne = true;
  f.genomeSize = 0;
  f.allOrderG.resize(0);
  f.allPosetG.resize(0);
  f.Poset.resize(0);
  f.Epistasis.resize(0);
  f.orderE.resize(0);
  f.Gene_Module_tabl.resize(0);
  f.allGenes.resize(0);
  f.drv.resize(0);
  f.fVars.resize(0);//new line to initialize fVars
  f.genesNoInt.shift = -99L;
  f.genesNoInt.NumID.resize(0);
  f.genesNoInt.names.resize(0);
  f.genesNoInt.s.resize(0);
  f.fitnessLandscape.NumID.resize(0);
  f.fitnessLandscape.names.resize(0);
  f.fitnessLandscape.flmap.clear();
  f.fitnessLandscape.flFDFmap.clear();//new line to initialize flFDFmap
  f.fitnessLandscape.flfVarsmap.clear();//new line to initialize flFDFmap
  return f;
}

// FIXME: fitness_as_genes and Genotype are identical
// structures. Why not use the same thing?
// Because even if just four vectors of ints, have different meaning.
// Humm...
struct fitness_as_genes {
  // fitnessEffectsAll in terms of genes.  Useful for output
  // conversions. There could be genes that are both in orderG and
  // posetEpistG. In such a case, only in orderG.
  // We only use a small part for now.
  // All are ordered vectors.
  std::vector<int> orderG;
  std::vector<int> posetEpistG;
  std::vector<int> noInt;
  std::vector<int> flGenes;
};

inline fitness_as_genes zero_fitness_as_genes() {
  fitness_as_genes g;
  g.orderG.resize(0);
  g.posetEpistG.resize(0);
  g.noInt.resize(0);
  g.flGenes.resize(0);
  return g;
}
// There are no shared genes in order and epist.  Any gene in orderEff can
// also be in the posets or general epistasis, but orderEff is only for
// those that have order effects.

// For all genes for which there are no order effects, any permutation of
// the same mutations is the same genotype, and has the same fitness. That
// is why we separate orderEff, which is strictly in the order in which
// mutations accumulate, and thus usorted, from the other effects, that
// are always kept sorted.

// rest are those genes that have no interactions. Evaluating their
// fitness is simple, and there can be no modules here.
struct Genotype {
  std::vector<int> orderEff;
  std::vector<int> epistRtEff; //always sorted
  std::vector<int> rest; // always sorted
  std::vector<int> flGenes; // always sorted; the fitness landscape genes
};





inline Genotype wtGenotype() {
  // Not needed but to make it explicit
  Genotype g;
  g.orderEff.resize(0);
  g.epistRtEff.resize(0);
  g.rest.resize(0);
  g.flGenes.resize(0);
  return g;
}

// struct st_PhylogNum {
//   double time;
//   std::vector<int> parent;
//   std::vector<int> child;
// };

// This is all we need to then use igraph on the data frame.
struct PhylogName {
  std::vector<double> time;
  std::vector<std::string> parent;
  std::vector<std::string> child;
  std::vector<double> pop_size_child;
  // yes, implicit constructor clears
};


// This is all we need to then use igraph on the data frame.
// simplified for the LOD that is always stored
struct LOD {
  // std::vector<double> time;
  std::vector<std::string> parent;
  std::vector<std::string> child;
};

// We only need the string, but if we store the genotype as such
// we can avoid a costly conversion that often leads to storing nothing
// in
struct POM {
  // std::vector<double> time;
  std::vector<std::string> genotypesString;
  std::vector<Genotype> genotypes;
};




std::vector<int> genotypeSingleVector(const Genotype& ge);

bool operator==(const Genotype& lhs, const Genotype& rhs);

// bool operator<(const Genotype& lhs, const Genotype& rhs);


TypeModel stringToModel(const std::string& dep);

Dependency stringToDep(const std::string& dep);

// std::string depToString(const Dependency dep);

void obtainMutations(const Genotype& parent,
		     const fitnessEffectsAll& fe,
		     int& numMutablePosParent,
		     std::vector<int>& newMutations,
		     //randutils::mt19937_rng& ran_gen
		     std::mt19937& ran_gen,
		     std::vector<double> mu);

Genotype createNewGenotype(const Genotype& parent,
			   const std::vector<int>& mutations,
			   const fitnessEffectsAll& fe,
			   std::mt19937& ran_gen,
			   //randutils::mt19937_rng& ran_gen
			   bool random);

std::vector<double> evalGenotypeFitness(const Genotype& ge,
					const fitnessEffectsAll& F);


fitnessEffectsAll convertFitnessEffects(Rcpp::List rFE);
std::vector<int> getGenotypeDrivers(const Genotype& ge, const std::vector<int>& drv);
std::vector<int> allGenesinGenotype(const Genotype& ge);
void print_Genotype(const Genotype& ge);

fitness_as_genes fitnessAsGenes(const fitnessEffectsAll& fe);

std::map<int, std::string> mapGenesIntToNames(const fitnessEffectsAll& fe);

std::vector<int> getGenotypeDrivers(const Genotype& ge, const std::vector<int>& drv);

double prodFitness(const std::vector<double>& s);

double prodDeathFitness(const std::vector<double>& s);

double mutationFromScratch(const std::vector<double>& mu,
			   const spParamsP& spP,
			   const Genotype& g,
			   const fitnessEffectsAll& fe,
			   const int mutationPropGrowth,
			   const std::vector<int> full2mutator,
			   const fitnessEffectsAll& muEF);

// double mutationFromParent(const std::vector<double>& mu,
// 			  const spParamsP& newP,
// 			  const spParamsP& parentP,
// 			  const std::vector<int>& newMutations,
// 			  // const std::vector<int>& nonmutated,
// 			  const int mutationPropGrowth,
// 			  const Genotype& fullge,
// 			  const std::vector<int> full2mutator,
// 			  const fitnessEffectsAll& muEF);

double prodMuts(const std::vector<double>& s);


double set_cPDetect(const double n2, const double p2,
		    const double PDBaseline);

bool detectedSizeP(const double n, const double cPDetect,
		   const double PDBaseline, std::mt19937& ran_gen);

std::vector < std::vector<int> > list_to_vector_of_int_vectors(Rcpp::List vlist);

void addToPOM(POM& pom,
	      const Genotype& genotype,
	      const std::map<int, std::string>& intName,
	      const fitness_as_genes& fg);

void addToPOM(POM& pom,
	      const std::string string);

#endif
