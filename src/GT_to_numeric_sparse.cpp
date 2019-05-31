// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath>

#include <iostream>
#include <string>

using namespace Rcpp;
using namespace std;



// Genotypes are in the form of 'a\b' or 'a|b'
// If genotypes are phased, make sum of a and b an entry of the returned matrix.
// If genotypes are not phased, make a and b sperate entries of the returned matrix.
// The returned matrix has subjects in rows and mutations in columns (different from the input matrix) 

// [[Rcpp::export]]
S4 GT_to_S4 (CharacterMatrix vcf_GT, bool phased = true) {
  
  int nrow = vcf_GT.ncol();
  int ncol = vcf_GT.nrow();
  List vcf_GT_dimnames = (Rf_isNull(vcf_GT.attr("dimnames")))? List::create(R_NilValue, R_NilValue): vcf_GT.attr("dimnames");
  List dimnames = List::create(vcf_GT_dimnames[1], vcf_GT_dimnames[0]);
  
  if (!phased) {
    
    NumericMatrix mat(nrow, ncol);
    
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) { 
        
        mat(i, j) = (vcf_GT(j, i)[0] - '0') + (vcf_GT(j, i)[2] - '0');
        
        if ( mat(i, j) < 0 || mat(i, j) > 2) {
            mat(i, j) = 0;
         } else {
            mat(i, j) =  mat(i, j) + 1;
         }
      }
    }
    
    RawMatrix rawmat = as<RawMatrix>(mat);
    
    rawmat.attr("dimnames") = dimnames;
    
    S4 snpmat("SnpMatrix");
    snpmat.slot(".Data") = rawmat;
    
    return snpmat;
    
  } else {
    
    typedef Eigen::SparseMatrix<double> SpMat; 
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(0.05*ncol*nrow);
    List dimnamesNew = List::create(R_NilValue, R_NilValue);
    
      
    for (int j = 0; j < nrow; j ++) {
      for (int i = 0; i < ncol; i++) {
        
        if ((vcf_GT(i, j)[0] == '0') && (vcf_GT(i, j)[2] == '0')) {
          continue;
        } else {
        
          if (vcf_GT(i, j)[0] == '1') {
            tripletList.push_back(T(j*2,i,1));
          } else if (vcf_GT(i, j)[0] == '.') {
            tripletList.push_back(T(j*2,i,NAN));
          }
        
          if (vcf_GT(i, j)[2] == '1') {
            tripletList.push_back(T(j*2+1,i,1));
          } else if (vcf_GT(i, j)[2] == '.') {
            tripletList.push_back(T(j*2+1,i,NAN));
          }
        
        }
      }
    }

 
    if (!Rf_isNull(dimnames[0])) {

      CharacterVector rowNames = dimnames[0];
      CharacterVector rowNamesNew(2*nrow);
      
      for (int i = 0; i < nrow; i++) {

        string str1(rowNames[i]);
        string str2(rowNames[i]);
        rowNamesNew[i*2] = (str1.append(".1")).c_str();
        rowNamesNew[i*2+1] = (str2.append(".2")).c_str();
      }
      
      dimnamesNew[0] = rowNamesNew;
    }
    
    if (!Rf_isNull(dimnames[1])) {

      dimnamesNew[1] = dimnames[1];

    }
    
    SpMat mat(2*nrow,ncol);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
  
    S4 dgcmat(wrap(mat));
    dgcmat.slot("Dimnames") = dimnamesNew;
    
    return dgcmat;
  }
  
  
}


