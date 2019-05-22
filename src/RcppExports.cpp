// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// GT_to_numeric
NumericMatrix GT_to_numeric(CharacterMatrix vcf_GT, bool phased);
RcppExport SEXP _LDheatmap_GT_to_numeric(SEXP vcf_GTSEXP, SEXP phasedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type vcf_GT(vcf_GTSEXP);
    Rcpp::traits::input_parameter< bool >::type phased(phasedSEXP);
    rcpp_result_gen = Rcpp::wrap(GT_to_numeric(vcf_GT, phased));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LDheatmap_GT_to_numeric", (DL_FUNC) &_LDheatmap_GT_to_numeric, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_LDheatmap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
