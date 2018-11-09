#include "batchelor.h"
#include "utils.h"

SEXP unlog_exprs_mean (SEXP mat, SEXP base, SEXP pseudo, SEXP subset) {
    BEGIN_RCPP
    
    auto ptr=beachmat::create_numeric_matrix(mat);
    const double log_base=check_numeric_scalar(base, "log base");
    const double pseudo_count=check_numeric_scalar(pseudo, "pseudo-count");
    auto schosen=check_subset_vector(subset, ptr->get_nrow());
    const double log_scale=std::log(log_base);

    const size_t slen=schosen.size();
    Rcpp::NumericVector out_means(slen);
    const size_t ncells=ptr->get_ncol();
    Rcpp::NumericVector tmp(ncells);

    for (size_t s=0; s<slen; ++s) {
        ptr->get_row(schosen[s], tmp.begin());
        for (auto& x : tmp) {
            x=std::exp(x*log_scale);
        }

        double cur_mean=std::accumulate(tmp.begin(), tmp.end(), 0.0);
        cur_mean /= ncells;
        cur_mean -= pseudo_count;
        out_means[s]=cur_mean;
    }

    return out_means;
    END_RCPP
}

SEXP unlog_exprs_scaled (SEXP mat, SEXP base, SEXP pseudo, SEXP subset, SEXP rescale) {
    BEGIN_RCPP

    auto ptr=beachmat::create_numeric_matrix(mat);
    const double log_base=check_numeric_scalar(base, "log base");
    const double pseudo_count=check_numeric_scalar(pseudo, "pseudo-count");
    const double log_scale=std::log(log_base), log_pseudo=std::log(pseudo_count);
    
    auto schosen=check_subset_vector(subset, ptr->get_nrow());
    size_t slen=schosen.size();
    Rcpp::NumericVector rescaling(rescale);
    if (slen!=rescaling.size()) {
        throw std::runtime_error("rescaling vector should have length equal to the subsetting vector");
    }

    const size_t ncells=ptr->get_ncol();
    Rcpp::NumericVector tmp(ncells);
    auto optr=beachmat::create_numeric_output(slen, ncells, 
        beachmat::output_param(mat, true, pseudo_count==1) // preserve sparsity if pseudo count is unity.
    );

    for (size_t s=0; s<slen; ++s) {
        ptr->get_row(schosen[s], tmp.begin());
        const double curscale=rescaling[s];

        for (auto& x : tmp) {
            x=std::expm1(x*log_scale - log_pseudo);
            if (x <= 0) { 
                x=0; 
            } else {
                // Technically, we need to multiply by pseudo_count to adjust for expm1();
                // and then re-divide, to enable log1p() to work. But these two cancel out
                // so we won't do either of them.
                x=std::log1p(x * curscale); 
            }
        }
        if (log_pseudo!=0) {
            for (auto& x : tmp) { x+=log_pseudo; }
        }
        for (auto& x : tmp) { x/=log_scale; }

        optr->set_row(s, tmp.begin());
    }

    return optr->yield();
    END_RCPP
}
