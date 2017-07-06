using namespace std;
// #include <limits.h>
#include <iostream>
// #include <math.h>
#include <values.h>

#define DP2(x) {std::cout << "\n Value of " << #x << " = " << x << std::endl;}

int main() {
  double value(0.0);
  double ti(1e-308); //DBL_MIN in my laptop is 2.22e-308

  std::string ti_dbl_comp;                                                                                                                                                                                                                                                                                              
  if( ti == DBL_MIN) {                                                                                                                                                                                                                                                                                                  
    ti_dbl_comp = "ti_equal_DBL_MIN";                                                                                                                                                                                                                                                                                   
  } else if (ti == 0.0) {                                                                                                                                                                                                                                                                                               
    ti_dbl_comp = "ti_equal_0.0";                                                                                                                                                                                                                                                                                       
  } else if ( (ti < DBL_MIN) && (ti > 0.0) ) {                                                                                                                                                                                                                                                                          
    ti_dbl_comp = "ti_gt_0.0_lt_DBL_MIN";                                                                                                                                                                                                                                                                               
  }  else {                                                                                                                                                                                                                                                                                                             
    ti_dbl_comp = "IMPOSSIBLE!"; // for my use cases                                                                                                                                                                                                                                                                                       
  }                                                                                                                                                                                                                                                                                                                     

  double a(0.000001);
  
  DP2(DBL_MIN);
  DP2(ti);
  DP2(value);
  DP2(ti_dbl_comp);

  DP2( ((a + ti) - a ) ); 
  DP2( (((a + ti) - a) < DBL_MIN) );
  return 0;
}
