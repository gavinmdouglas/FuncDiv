// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// prep_all_sample_func_vec
List prep_all_sample_func_vec(NumericMatrix abun_tab, NumericMatrix func_tab);
RcppExport SEXP _FuncDiv_prep_all_sample_func_vec(SEXP abun_tabSEXP, SEXP func_tabSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type abun_tab(abun_tabSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type func_tab(func_tabSEXP);
    rcpp_result_gen = Rcpp::wrap(prep_all_sample_func_vec(abun_tab, func_tab));
    return rcpp_result_gen;
END_RCPP
}
// prep_all_sample_func_taxa_vec
List prep_all_sample_func_taxa_vec(NumericMatrix abun_tab, NumericMatrix func_tab);
RcppExport SEXP _FuncDiv_prep_all_sample_func_taxa_vec(SEXP abun_tabSEXP, SEXP func_tabSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type abun_tab(abun_tabSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type func_tab(func_tabSEXP);
    rcpp_result_gen = Rcpp::wrap(prep_all_sample_func_taxa_vec(abun_tab, func_tab));
    return rcpp_result_gen;
END_RCPP
}
// prep_func_contributor_dimnames
List prep_func_contributor_dimnames(arma::mat abun_tab, arma::mat func_tab);
RcppExport SEXP _FuncDiv_prep_func_contributor_dimnames(SEXP abun_tabSEXP, SEXP func_tabSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type abun_tab(abun_tabSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type func_tab(func_tabSEXP);
    rcpp_result_gen = Rcpp::wrap(prep_func_contributor_dimnames(abun_tab, func_tab));
    return rcpp_result_gen;
END_RCPP
}
// rbiom_par_unifrac
NumericVector rbiom_par_unifrac(List sparseMatrix, List tree, IntegerVector weighted);
RcppExport SEXP _FuncDiv_rbiom_par_unifrac(SEXP sparseMatrixSEXP, SEXP treeSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sparseMatrix(sparseMatrixSEXP);
    Rcpp::traits::input_parameter< List >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(rbiom_par_unifrac(sparseMatrix, tree, weighted));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FuncDiv_prep_all_sample_func_vec", (DL_FUNC) &_FuncDiv_prep_all_sample_func_vec, 2},
    {"_FuncDiv_prep_all_sample_func_taxa_vec", (DL_FUNC) &_FuncDiv_prep_all_sample_func_taxa_vec, 2},
    {"_FuncDiv_prep_func_contributor_dimnames", (DL_FUNC) &_FuncDiv_prep_func_contributor_dimnames, 2},
    {"_FuncDiv_rbiom_par_unifrac", (DL_FUNC) &_FuncDiv_rbiom_par_unifrac, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_FuncDiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
