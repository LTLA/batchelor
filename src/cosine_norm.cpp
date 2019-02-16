#include "batchelor.h"

#include "utils.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <algorithm>
#include <cmath>

/* Performs the cosine normalization in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (M mat, SEXP original, SEXP return_mat) {
    const size_t& nrow=mat->get_nrow();
    const size_t& ncol=mat->get_ncol();

    // Deciding whether or not to return the matrix.
    bool mat_return=check_logical_scalar(return_mat, "return matrix specification");
    
    beachmat::numeric_output* optr=NULL;
    std::unique_ptr<beachmat::numeric_output> holder;
    if (mat_return) { 
        holder=beachmat::create_numeric_output(nrow, ncol, beachmat::output_param(original, false, true));
        optr=holder.get();
    }
   
    // Calculating the L2 norm of each vector and applying it. 
    Rcpp::NumericVector incoming(nrow);
    Rcpp::NumericVector l2norm(ncol);
    for (size_t c=0; c<ncol; ++c) {
        mat->get_col(c, incoming.begin());

        double& total=l2norm[c];
        for (const auto& val : incoming) { 
            total+=val*val;
        }
        total=std::sqrt(total);

        if (mat_return) { 
            const double pos_total=std::max(total, 0.00000001); // avoid division by zero.
            for (auto& val : incoming) { 
                val/=pos_total;
            }
            optr->set_col(c, incoming.begin());
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
        auto input=beachmat::create_integer_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming, getmat);
    } else {
        auto input=beachmat::create_numeric_matrix(incoming);
        return cosine_norm_internal(input.get(), incoming, getmat);
    }
    END_RCPP
}
