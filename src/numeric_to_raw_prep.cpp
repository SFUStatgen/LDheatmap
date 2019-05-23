#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix numeric_to_raw_prep(NumericMatrix mat) {
  
  int mat_size = mat.nrow() * mat.ncol();
  
  
  for (int i = 0; i < mat_size; i++) {
    
    if (mat[i] < 0 || mat[i] > 2) {
      mat[i] = 0;
      // Rcout << "values other than 0, 1 or 2 set to NA" << endl;
    } else {
      mat[i] = mat[i] + 1;
    }
    
  }
  
  return mat;
  
  
}

