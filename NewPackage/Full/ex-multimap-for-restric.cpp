#include <Rcpp.h>

using namespace Rcpp ;


// struct Deps {
//   std::vector<int> module;
//   double s;
//   double sh;
// };

// std::multimap<int, Deps> restrictTable;


// FIXME: check later penalty for using a string for typeDep
struct geneDeps {
  int typeDep; // smaller, predictable size. A lot less readable, though.
  double s;
  double sh;
  // std::vector<int> deps;
  std::vector< std::vector<int> > deps;
};


// // [[Rcpp::export]]
// void f2(){
//   geneDeps g1;
//   g1.s = 0.3;
//   g1.sh = 0.6;
//   g1.deps.push_back(3);
// }

// // [[Rcpp::export]]
// void f3(){
//   std::vector<geneDeps> g3;
//   g3.resize(4);
//   std::cout << "g3 size " << g3.size() << std::endl;
//   std::cout << "g3[0].sh " << g3[0].sh << std::endl;
//   g3[0].s = 0.3;
//   g3[0].sh = 0.6;
//   std::cout << "g3[0].sh " << g3[0].sh << std::endl;

//   std::cout << "g3[0].deps.size() " << g3[0].deps.size() << std::endl;
//   g3[0].deps.push_back(3);

//   std::cout << "g3[0].deps.size() " << g3[0].deps.size() << std::endl;
// }

// [[Rcpp::export]]
void f4(){
  std::vector<geneDeps> g3;
  g3.resize(4);
  g3[0].s = 0.3;
  g3[0].sh = 0.6;
  g3[0].deps.resize(2);
  // g3[0].deps[0].push_back(std::vector<int>(3));
  std::cout << "g3[0].deps.size() " << g3[0].deps.size() << std::endl;
  std::cout << "g3[0].deps[0].size() " << 
    g3[0].deps[0].size() << std::endl;
  g3[0].deps[0].push_back(3);
  std::cout << "g3[0].deps[0].size() " << 
    g3[0].deps[0].size() << std::endl;

  std::vector<int>vv(4, 3);
  g3[1].deps.push_back(vv);
  std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
  g3[1].deps.push_back(std::vector<int>(9, 5));
  std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
  
  for (auto c : g3[1].deps[0])
    std::cout << c << ' ';
  
  for (auto c : g3[1].deps[1])
    std::cout << c << ' ';
}



// [[Rcpp::export]]
void restrictTable_to_cpp(Rcpp::List rt,
			  std::vector<geneSeps>& restrictTable) { 
  int ndeps;

  if(restrictTable.size())
    restrictTable.clear(); //not needed later!! FIXME

  //std::vector<geneDeps> restrictTable;
  restrictTable.resize(rt.size());
  
  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;
  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    restrictTable[i].typeDep = rt_element["type"];
    restrictTable[i].s = rt_element["s"];
    restrictTable[i].sh = rt_element["sh"];
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      restrictTable[i].deps.push_back(Rcpp::as< std::vector<int> >(module));

      // std::cout <<"\n ** \n j = " << j << std::endl;
      // for (auto c : restrictTable[i].deps[j])
      // 	std::cout << c << ' ';
      // std::cout << std::endl;

    }
  }
}


void printRestrictTable(const std::vector<geneSeps>& restrictTable) {
  Rcpp::Rcout << "\n **********  Restriction table *******" << std::endl;
  Rcpp::Rcout << "\t Size = ", restrictTable.size() << std::endl;
  for(size_t i = 0; i != restrictTable.size(); ++i) {
    Rcpp::Rcout <<"\t\t Dependent node = " << i << std::endl;
    Rcpp::Rcout <<"\t\t\t typeDep = " << restrictTable[i].typeDep << " ";
    Rcpp::Rcout <<"\t s = " << restrictTable[i].s << " ";
    Rcpp::Rcout <<"\t sh = " << restrictTable[i].sh << std::endl;
    // here the code for parent modules


    
  }

}

// just a test line

// an a another in bufo


// // [[Rcpp::export]]
// void f1() {
//   //int n = 5;
//   std::vector<int> v1;
//   for(int i = 0; i != 5; ++i) {
//     for(int j = 0; j != (i + 1); ++j)
//       v1.push_back(j);
//     Deps d1;
//     d1.module = v1;
//     d1.s = 0.1 + i;
//     d1.sh = 0.6 + i;

//     restrictTable.insert(std::make_pair(i, d1));
//     v1.clear();
    
//   }
//   std::cout << "\n The multimap contains \n";
//   std::multimap<int, Deps>::iterator it;
//   for(it = restrictTable.begin(); it != restrictTable.end(); ++it ) {
//     std::cout << it->first << " => " << " Module = ";
//     for(auto vv : it->second.module) std::cout << vv << ' ';
//     std::cout << " s = " << it->second.s << "; sh = " << it->second.sh << std::endl;
//   }
//   restrictTable.clear();
// }




