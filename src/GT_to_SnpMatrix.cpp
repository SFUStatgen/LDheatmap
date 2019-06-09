#include "Rcpp.h"
#include <iostream>
#include <string>

using namespace Rcpp;
using namespace std;


// testing, first asume phased = T
// [[Rcpp::export]]
S4 GT_to_SnpMatrix (CharacterMatrix GT, bool phased = true) {
  
  int nrow = GT.ncol();
  int ncol = GT.nrow();
  List GT_dimnames = (Rf_isNull(GT.attr("dimnames")))? List::create(R_NilValue, R_NilValue): GT.attr("dimnames");
  List dimnames = List::create(GT_dimnames[1], GT_dimnames[0]);

  if (!phased) {
    
    RawMatrix mat(nrow, ncol);
    
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) { 
      
        int num = (GT(j, i)[0] - '0') + (GT(j, i)[2] - '0') + 1;
        mat(i,j) = (num > 0 && num < 4)? num: 0; 
        
      }
    }
    
    S4 S("SnpMatrix");
    S.slot(".Data") = mat;
    S.attr("dimnames") = dimnames;
    
    return S;
    
  } else {  
  
  
    RawMatrix mat(2*nrow, ncol);
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
          
          int num1 = GT(j, i)[0] - '0';
          int num2 = GT(j, i)[2] - '0';
          
          if (num1 == 0) {
            mat(i*2, j) = 1;
          } else if (num1 == 1) {
            mat(i*2, j) = 3;
          } else {
            mat(i*2, j) = 0;
          }
          
          if (num2 == 0) {
            mat(i*2+1, j) = 1;
          } else if (num2 == 1) {
            mat(i*2+1, j) = 3;
          } else {
            mat(i*2+1, j) = 0;
          }

        }
      }

      dimnamesNew[0] = rowNamesNew;


    } else {
    
      for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
          
          int num1 = GT(j, i)[0] - '0';
          int num2 = GT(j, i)[2] - '0';
          
          if (num1 == 0) {
            mat(i*2, j) = 1;
          } else if (num1 == 1) {
            mat(i*2, j) = 3;
          } else {
            mat(i*2, j) = 0;
          }
          
          if (num2 == 0) {
            mat(i*2+1, j) = 1;
          } else if (num2 == 1) {
            mat(i*2+1, j) = 3;
          } else {
            mat(i*2+1, j) = 0;
          }

        }
      }
    }

    if (!Rf_isNull(dimnames[1])) {
      dimnamesNew[1] = dimnames[1];
    }
      
    S4 XS("XSnpMatrix");
    XS.slot("diploid") = LogicalVector(nrow*2, FALSE);
    XS.slot(".Data") = mat;
    XS.attr("dimnames") = dimnamesNew;
    
    return XS;
  }
}



    
  