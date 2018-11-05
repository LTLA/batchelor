#include "batchelor.h"
#include "utils.h"

/* Performs the cosine normalization in a fairly efficient manner. */

template<class M>
SEXP cosine_norm_internal (M mat, SEXP original, SEXP return_mat) {
    const size_t& nrow=mat->get_nrow();
    const size_t& ncol=mat->get_ncol();

    // Deciding whether or not to return the matrix.
    bool mat_return=check_logical_scalar(return_mat, "return matrix specification");
    
    beachmat::numeric_output* optr=NULL;
    std::vector<std::unique_ptr<beachmat::numeric_output> > holder;
    if (mat_return) { 
        holder.push_back(beachmat::create_numeric_output(nrow, ncol, 
                         beachmat::output_param(original, false, true)));
        optr=holder.front().get();
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
        total=std::max(total, 0.00000001); // avoid division by zero.

        if (mat_return) { 
            for (auto& val : incoming) { 
                val/=total;
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
