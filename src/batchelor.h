#ifndef BATCHELOR_H
#define BATCHELOR_H

#include "Rcpp.h"

extern "C" {

// MNN calculations.

SEXP find_mutual_nns(SEXP, SEXP);

SEXP smooth_gaussian_kernel(SEXP, SEXP, SEXP, SEXP);

SEXP adjust_shift_variance(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// Unlogging calculations

SEXP unlog_exprs_mean(SEXP, SEXP, SEXP, SEXP);

SEXP unlog_exprs_scaled(SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif 
