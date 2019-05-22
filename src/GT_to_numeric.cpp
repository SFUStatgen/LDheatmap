#include "Rcpp.h"
#include <iostream>
#include <string>

using namespace Rcpp;
using namespace std;




// void createUnphasedNumeric (NumericMatrix &mat, CharacterMatrix vcf_GT){}

// [[Rcpp::export]]
NumericMatrix GT_to_numeric (CharacterMatrix vcf_GT, bool phased = true) {
  
  int nrow = vcf_GT.nrow();
  int ncol = vcf_GT.ncol();
  List dimnames = vcf_GT.attr("dimnames");
  
  if (phased != true) {
    
    NumericMatrix mat(nrow, ncol);
    int mat_size = mat.nrow() * mat.ncol();
    
    for (int i = 0; i < mat_size; i++) {
      
      mat[i] = (vcf_GT[i][0] - '0') + (vcf_GT[i][2] - '0');
      
    }
    
    mat.attr("dimnames") = dimnames;
    
    return mat;
    
  } else {
    

    NumericMatrix mat(nrow, 2*ncol);
    List dimnamesNew = List::create(R_NilValue, R_NilValue);
    
    if (!Rf_isNull(dimnames[1])) {
      
      CharacterVector colNames = dimnames[1];
      CharacterVector colNamesNew(2*ncol);
      
      for (int j = 0; j < ncol; j++) {
        
        string str1(colNames[j]);
        string str2(colNames[j]);
        colNamesNew[j*2] = (str1.append(".1")).c_str();
        colNamesNew[j*2+1] = (str2.append(".2")).c_str();
        
        for (int i = 0; i < nrow; i++) { 
          
          mat(i, j*2) = vcf_GT(i, j)[0] - '0';
          mat(i, j*2+1) =  vcf_GT(i, j)[2] - '0';
          
        }
      }
      
      if (!Rf_isNull(dimnames[0])) {
        
        dimnamesNew[0] = dimnames[0];
        
      }
      
      dimnamesNew[1] = colNamesNew;
      mat.attr("dimnames") = dimnamesNew;
      
      return mat;
      
    } else {
    
      for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) { 
        
        mat(i, j*2) = vcf_GT(i, j)[0] - '0';
        mat(i, j*2+1) =  vcf_GT(i, j)[2] - '0';
        
        }
      }
      
      if (!Rf_isNull(dimnames[0])) {
        
        dimnamesNew[0] = dimnames[0];
        
      }
      
      mat.attr("dimnames") = dimnamesNew;
      return mat;
    }
  }
  
  
}