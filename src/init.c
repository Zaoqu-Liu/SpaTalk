#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _SpaTalk_cpp_coexp_fast(SEXP, SEXP);
extern SEXP _SpaTalk_cpp_permutation_test(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpaTalk_cpp_random_walk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpaTalk_cpp_batch_coexp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _SpaTalk_cpp_fast_dist(SEXP, SEXP);
extern SEXP _SpaTalk_cpp_knn(SEXP, SEXP, SEXP);
extern SEXP _SpaTalk_cpp_batch_cor(SEXP, SEXP);
extern SEXP _SpaTalk_cpp_fast_sampling(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_SpaTalk_cpp_coexp_fast", (DL_FUNC) &_SpaTalk_cpp_coexp_fast, 2},
    {"_SpaTalk_cpp_permutation_test", (DL_FUNC) &_SpaTalk_cpp_permutation_test, 7},
    {"_SpaTalk_cpp_random_walk", (DL_FUNC) &_SpaTalk_cpp_random_walk, 7},
    {"_SpaTalk_cpp_batch_coexp", (DL_FUNC) &_SpaTalk_cpp_batch_coexp, 4},
    {"_SpaTalk_cpp_fast_dist", (DL_FUNC) &_SpaTalk_cpp_fast_dist, 2},
    {"_SpaTalk_cpp_knn", (DL_FUNC) &_SpaTalk_cpp_knn, 3},
    {"_SpaTalk_cpp_batch_cor", (DL_FUNC) &_SpaTalk_cpp_batch_cor, 2},
    {"_SpaTalk_cpp_fast_sampling", (DL_FUNC) &_SpaTalk_cpp_fast_sampling, 7},
    {NULL, NULL, 0}
};

void R_init_SpaTalk(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

