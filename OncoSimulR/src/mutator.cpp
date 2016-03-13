// [[Rcpp::export]]
double evalRGenotypeMut(Rcpp::IntegerVector rG,
		       Rcpp::List rFE,
		       bool verbose) {
  // From evalRGenotype. Basically same thing except for a minor output.
  // So maybe just generalize previous?
  
  if(rG.size() == 0) {
    Rcpp::warning("WARNING: you have evaluated mutator genotype of a genotype of length zero.");
    return 1;
  }
    
  const Rcpp::List rF(rFE);
  fitnessEffectsAll F = convertFitnessEffects(rF);
  Genotype g = convertGenotypeFromR(rG, F);
  vector<double> s = evalGenotypeMut(g, F);
  if(verbose) {
    Rcpp::Rcout << "\n Individual mutator terms are :";
    for(auto const &i : s) Rcpp::Rcout << " " << i;
    Rcpp::Rcout << std::endl;
  }
  
  return prodMuts(s);
}

double prodMuts(const std::vector<double>& s) {
  // From prodFitness
  // return accumulate(s.begin(), s.end(), 1.0,
  // 		    [](double x, double y) {return (x * y);});
  return accumulate(s.begin(), s.end(), 1.0,
		    std::multiplies<double>());
}


std::vector<double> evalGenotypeMut(const Genotype& ge,
				    const fitnessEffectsAll& F){

  // From evalGenotypeFitness:
  // FIXME: ehhh? this is the exact same thing. Do we need a new function?


  // print_Genotype(ge);
  // check_disable_later
  checkLegitGenotype(ge, F);

  std::vector<double> s;
  if( (ge.orderEff.size() + ge.epistRtEff.size() + ge.rest.size()) == 0) {
    Rcpp::warning("WARNING: you have evaluated fitness of a genotype of length zero.");
    s.push_back(1.0);
    return s;
  }

  // Genes without any restriction or epistasis are just genes. No modules.
  // So simple we do it here.
  if(F.genesNoInt.shift > 0) {
    int shift = F.genesNoInt.shift;
    for(auto const & r : ge.rest ) {
      s.push_back(F.genesNoInt.s[r - shift]);
    }
  }

  // For the rest, there might be modules. Three different effects on
  // fitness possible: as encoded in Poset, general epistasis, order effects.
  
  // Epistatis and poset are checked against all mutations. Create single
  // sorted vector with all mutations and map to modules, if needed. Then
  // eval.

  // Why not use a modified genotypeSingleVector without the no ints? We
  // could, but not necessary. And you can place genes in any order you
  // want, since this is not for order restrictions. That goes below.
  // Why do I put the epist first? Se previous answer.
  // Why do I sort if one to one? binary searches. Not done below for order.
  std::vector<int> mutG (ge.epistRtEff);
  mutG.insert( mutG.end(), ge.orderEff.begin(), ge.orderEff.end());
  std::vector<int> mutatedModules;
  if(F.gMOneToOne) {
    sort(mutG.begin(), mutG.end()); 
    mutatedModules = mutG;
  } else {
    mutatedModules = GeneToModule(mutG, F.Gene_Module_tabl, true, true);
  }
  std::vector<double> srt =
    evalPosetConstraints(mutatedModules, F.Poset, F.allPosetG);
  std::vector<double> se =
    evalEpistasis(mutatedModules, F.Epistasis);

  // For order effects we need a new vector of mutatedModules:
  if(F.gMOneToOne) {
    mutatedModules = ge.orderEff;
  } else {
    mutatedModules = GeneToModule(ge.orderEff, F.Gene_Module_tabl, false, true);
  }
  
  std::vector<double> so =
    evalOrderEffects(mutatedModules, F.orderE);

  // I keep s, srt, se, so separate for now for debugging.
  s.insert(s.end(), srt.begin(), srt.end());
  s.insert(s.end(), se.begin(), se.end());
  s.insert(s.end(), so.begin(), so.end());
  
  return s;
}
