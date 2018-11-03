// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// MI_categorical_worker_two
double MI_categorical_worker_two(const Rcpp::IntegerMatrix& d);
RcppExport SEXP _bigmi_MI_categorical_worker_two(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(MI_categorical_worker_two(d));
    return rcpp_result_gen;
END_RCPP
}
// MI_categorical_worker_all
Rcpp::NumericMatrix MI_categorical_worker_all(const Rcpp::IntegerMatrix& d);
RcppExport SEXP _bigmi_MI_categorical_worker_all(SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(MI_categorical_worker_all(d));
    return rcpp_result_gen;
END_RCPP
}
// MI_categorical_worker_sparse_all
Rcpp::List MI_categorical_worker_sparse_all(const Rcpp::IntegerMatrix& d, const bool progress);
RcppExport SEXP _bigmi_MI_categorical_worker_sparse_all(SEXP dSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(MI_categorical_worker_sparse_all(d, progress));
    return rcpp_result_gen;
END_RCPP
}
// pearson_correlation_sparse_all
Rcpp::List pearson_correlation_sparse_all(const Rcpp::NumericMatrix& d, const bool progress);
RcppExport SEXP _bigmi_pearson_correlation_sparse_all(SEXP dSEXP, SEXP progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const bool >::type progress(progressSEXP);
    rcpp_result_gen = Rcpp::wrap(pearson_correlation_sparse_all(d, progress));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bigmi_MI_categorical_worker_two", (DL_FUNC) &_bigmi_MI_categorical_worker_two, 1},
    {"_bigmi_MI_categorical_worker_all", (DL_FUNC) &_bigmi_MI_categorical_worker_all, 1},
    {"_bigmi_MI_categorical_worker_sparse_all", (DL_FUNC) &_bigmi_MI_categorical_worker_sparse_all, 2},
    {"_bigmi_pearson_correlation_sparse_all", (DL_FUNC) &_bigmi_pearson_correlation_sparse_all, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_bigmi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
