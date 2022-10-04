#include <Rcpp.h>
using namespace Rcpp;

// double DEG_EW_TO_KM = 111.320;
// double DEG_NS_TO_KM = 110.574;
// double RES = 4;
// double D = 100;
// double DT = 7*24*60*60;
// 
// // [[Rcpp::export]]
// double h_lon_cpp(double lat){
//   return(1000*DEG_EW_TO_KM*cos(lat*PI/180)/RES);
// }

// [[Rcpp::export]]
NumericMatrix conv_cpp(Function H, NumericMatrix Y, double p_length=0){
  int n = Y.nrow();
  int m = Y.ncol();
  NumericMatrix XconvY(n,m);
  
  for(int i=0; i<n; i++){
    // grab the R function
    NumericMatrix X = as<NumericMatrix>( H(i,p_length) );
    int mid_i = X.nrow()/2;
    int mid_j = X.ncol()/2;
    
    for(int j=0; j<m; j++){
	 
      double conv_sum = 0;
      XconvY(i,j) = 0;
      for(int s=0; s<X.nrow(); s++){
        for(int t=0; t<X.ncol(); t++){
           if(((mid_i+i-s) >= 0) &&
              ((mid_j+j-t) >= 0) &&
              ((mid_i+i-s) < n) && 
              ((mid_j+j-t) < m))
          {
            conv_sum = conv_sum + X(s,t) * Y(mid_i+i-s,mid_j+j-t);
          }
        }
      }
      XconvY(i,j) = conv_sum;
    }
  }
  return(XconvY);
}

