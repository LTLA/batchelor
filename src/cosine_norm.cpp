#include "batchelor.h"

#include "utils.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <algorithm>
#include <cmath>

/* Compute the L2 norm in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (SEXP input, SEXP subset_row) {
    auto mat=beachmat::create_matrix<M>(input);
    const size_t nrow=mat->get_nrow();
    const size_t ncol=mat->get_ncol();

    Rcpp::IntegerVector sub=check_subset_vector(subset_row, nrow);
    Rcpp::NumericVector l2norm(ncol), incol(nrow);

    for (size_t c=0; c<ncol; ++c) {
        mat->get_col(c, incol.begin());

        double& total=l2norm[c];
        for (auto s : sub) {
            const auto& i=incol[s];
            total+=i*i;
        }
        total=std::sqrt(total);
    }

    return l2norm;
}

SEXP cosine_norm(SEXP incoming, SEXP subset_row) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        return cosine_norm_internal<beachmat::integer_matrix>(incoming, subset_row);
    } else {
        return cosine_norm_internal<beachmat::numeric_matrix>(incoming, subset_row);
    }
    END_RCPP
}
