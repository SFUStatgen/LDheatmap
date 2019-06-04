#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix raw_to_numeric(RawMatrix mat) {
  
  int mat_size = mat.nrow() * mat.ncol();
  NumericMatrix num_mat(mat.nrow(), mat.ncol());
  
  for (int i = 0; i < mat_size; i++) {
    
    if ((int)mat[i] == 0) {
      num_mat[i] = NA_REAL;
      // Rcout << "values other than 0, 1 or 2 set to NA" << endl;
    } else {
      num_mat[i] = (int)mat[i];
    }
    
  }
  
  return num_mat;
  
  
}
