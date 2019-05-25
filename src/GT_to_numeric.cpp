#include "Rcpp.h"
#include <iostream>
#include <string>

using namespace Rcpp;
using namespace std;



// Genotypes are in the form of 'a\b' or 'a|b'
// If genotypes are phased, make sum of a and b an entry of the returned matrix.
// If genotypes are not phased, make a and b sperate entries of the returned matrix.

// [[Rcpp::export]]
NumericMatrix GT_to_numeric (CharacterMatrix vcf_GT, bool phased = true) {
  
  int nrow = vcf_GT.nrow();
  int ncol = vcf_GT.ncol();
  List dimnames = (Rf_isNull(vcf_GT.attr("dimnames")))? List::create(R_NilValue, R_NilValue): vcf_GT.attr("dimnames");
  
  if (!phased) {
    
    NumericMatrix mat(nrow, ncol);
    int mat_size = mat.nrow() * mat.ncol();
    
    for (int i = 0; i < mat_size; i++) {
      
      mat[i] = (vcf_GT[i][0] - '0') + (vcf_GT[i][2] - '0');
      
    }
    
    mat.attr("dimnames") = dimnames;
    
    return mat;
    
  } else {
    

    NumericMatrix mat(2*nrow, ncol);
    List dimnamesNew = List::create(R_NilValue, R_NilValue);
    
    if (!Rf_isNull(dimnames[0])) {
      
      CharacterVector rowNames = dimnames[0];
      CharacterVector rowNamesNew(2*nrow);
      
      for (int i = 0; i < nrow; i++) {
        
        string str1(rowNames[i]);
        string str2(rowNames[i]);
        rowNamesNew[i*2] = (str1.append(".1")).c_str();
        rowNamesNew[i*2+1] = (str2.append(".2")).c_str();
        
        for (int j = 0; j < ncol; j++) { 
          
          mat(i*2, j) = vcf_GT(i, j)[0] - '0';
          mat(i*2+1, j) =  vcf_GT(i, j)[2] - '0';
          
        }
      }
      
      dimnamesNew[0] = rowNamesNew;

      
    } else {
    
      for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) { 
        
        mat(i, j*2) = vcf_GT(i, j)[0] - '0';
        mat(i, j*2+1) =  vcf_GT(i, j)[2] - '0';
        
        }
      }

    }
    
    
    if (!Rf_isNull(dimnames[1])) {
      
      dimnamesNew[1] = dimnames[1];
      
    }
    
    mat.attr("dimnames") = dimnamesNew;
    
    return mat;
    
    
  }
  
  
}