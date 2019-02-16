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



#ifndef _COMMON_CLS_H_
#define _COMMON_CLS_H_

#include<Rcpp.h>

// Simple custom exception for exceptions that lead to re-runs.
class rerunExcept: public std::runtime_error {
public:
  rerunExcept(const std::string &s) :
    std::runtime_error(s) {}
};

enum class TypeModel {exp, bozic1, mcfarlandlog};


struct spParamsP {
  double popSize;
  double birth;
  double death;
  double W;
  double R;
  double mutation; 
  double timeLastUpdate;
  std::multimap<double, int>::iterator pv;
  double absfitness; //convenient for Beerenwinkel
  int numMutablePos; //for mutator if need update of mutation
};



#endif

