#include "batchelor.h"

#include "utils.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"

#include <algorithm>
#include <cmath>

/* Performs the cosine normalization in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (SEXP input, SEXP return_mat) {
    auto mat=beachmat::create_matrix<M>(input);
    const size_t& nrow=mat->get_nrow();
    const size_t& ncol=mat->get_ncol();

    // Deciding whether or not to return the matrix.
    const bool mat_return=check_logical_scalar(return_mat, "return matrix specification");
    std::unique_ptr<beachmat::numeric_output> optr=nullptr;
    if (mat_return) { 
        optr=beachmat::create_numeric_output(nrow, ncol, beachmat::output_param(mat.get()));
    }
   
    // Calculating the L2 norm of each vector and applying it. 
    beachmat::const_column<M> incoming(mat.get());
    Rcpp::NumericVector l2norm(ncol), outcol(nrow);

    for (size_t c=0; c<ncol; ++c) {
        incoming.fill(c);
        auto n=incoming.get_n();
        auto val=incoming.get_values();

        double& total=l2norm[c];
        for (size_t i=0; i<n; ++i, ++val) {
            total+=(*val)*(*val);
        }
        total=std::sqrt(total);

        if (mat_return) { 
            const double pos_total=std::max(total, 0.00000001); // avoid division by zero.
            auto val=incoming.get_values();
            auto oIt=outcol.begin();
            for (size_t i=0; i<n; ++i, ++val, ++oIt) {
                *(oIt)=*val/pos_total;
            }

            if (incoming.is_sparse()) {
                auto idx=incoming.get_indices();
                optr->set_col_indexed(c, n, idx, outcol.begin());
            } else {
                optr->set_col(c, outcol.begin());
            }
        }
    }

    // Figuring out what to return.
    if (mat_return) { 
        return Rcpp::List::create(optr->yield(), l2norm);
    } else {
        return Rcpp::List::create(R_NilValue, l2norm);
    }
}

SEXP cosine_norm(SEXP incoming, SEXP getmat) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(incoming);
    if (rtype==INTSXP) {
        return cosine_norm_internal<beachmat::integer_matrix>(incoming, getmat);
    } else {
        return cosine_norm_internal<beachmat::numeric_matrix>(incoming, getmat);
    }
    END_RCPP
}
