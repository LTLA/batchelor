#include "batchelor.h"
#include "utils.h"

/* Perform smoothing with the Gaussian kernel */

SEXP smooth_gaussian_kernel(SEXP vect, SEXP index, SEXP data, SEXP sigma) {
    BEGIN_RCPP
    auto correction_vectors=beachmat::create_numeric_matrix(vect);
    const size_t npairs=correction_vectors->get_nrow();
    const size_t ngenes=correction_vectors->get_ncol();
    const Rcpp::IntegerVector _index(index);
    if (npairs!=_index.size()) { 
        throw std::runtime_error("number of rows in 'vect' should be equal to length of 'index'");
    }
    
    // Constructing the average vector for each MNN cell.
    std::deque<Rcpp::NumericVector> averages;
    std::set<int> mnncell;
    {
        std::deque<int> number;
        Rcpp::NumericVector currow(ngenes);
        int row=0;
        for (const auto& i : _index) {
            correction_vectors->get_row(row, currow.begin());
            ++row;
            
            if (i >= averages.size() || averages[i].size()==0) { 
                if (i >= averages.size()) { 
                    averages.resize(i+1);
                    number.resize(i+1);
                }
                averages[i]=currow;
                number[i]=1;
                mnncell.insert(i);
                currow=Rcpp::NumericVector(ngenes);
            } else {
                auto& target=averages[i];
                auto cIt=currow.begin();
                for (auto& t : target){
                    t += *cIt;
                    ++cIt;
                }
                ++(number[i]);
            }
        }
    
        for (const auto& i : mnncell) {
            auto& target=averages[i];
            const auto& num=number[i];
            for (auto& t : target) {
                t/=num;
            }
        }
    }

    // Setting up input constructs (including the expression matrix on which the distances are computed).
    auto mat=beachmat::create_numeric_matrix(data);
    const int ncells=mat->get_ncol();
    const int ngenes_for_dist=mat->get_nrow();
    const double s2=check_numeric_scalar(sigma, "sigma");

    // Setting up output constructs.
    Rcpp::NumericMatrix output(ngenes, ncells); // yes, this is 'ngenes' not 'ngenes_for_dist'.
    std::vector<double> distances2(ncells), totalprob(ncells);
    Rcpp::NumericVector mnn_incoming(ngenes_for_dist), other_incoming(ngenes_for_dist);

    // Using distances between cells and MNN-involved cells to smooth the correction vector per cell.
    for (const auto& mnn : mnncell) {
        auto mnn_iIt=mat->get_const_col(mnn, mnn_incoming.begin());

        for (int other=0; other<ncells; ++other) {
            double& curdist2=(distances2[other]=0);
            auto other_iIt=mat->get_const_col(other, other_incoming.begin());
            auto iIt_copy=mnn_iIt;

            for (int g=0; g<ngenes_for_dist; ++g) {
                const double tmp=(*iIt_copy  - *other_iIt);
                curdist2+=tmp*tmp;
                ++other_iIt;
                ++iIt_copy;
            }
        }

        // Compute log-probabilities using a Gaussian kernel based on the squared distances.
        // We keep things logged to avoid float underflow, and just ignore the constant bit at the front.
        for (auto& d2 : distances2) { 
            d2/=-s2;
        }
        
        // We sum the probabilities to get the relative MNN density. 
        // This requires some care as the probabilities are still logged at this point.
        double density=NA_REAL;
        for (const auto& other_mnn : mnncell) {
            if (ISNA(density)) {
                density=distances2[other_mnn];
            } else {
                const double& to_add=distances2[other_mnn];
                const double larger=std::max(to_add, density), diff=std::abs(to_add - density);
                density=larger + log1pexp(-diff);
            }
        }

        // Each correction vector is weighted by the Gaussian probability (to account for distance)
        // and density (to avoid being dominated by high-density regions).
        // Summation (and then division, see below) yields smoothed correction vectors.
        const auto& correction = averages[mnn];
        auto oIt=output.begin();
        for (int other=0; other<ncells; ++other) {
            const double mult=std::exp(distances2[other] - density);
            totalprob[other]+=mult;

            for (const auto& corval : correction) { 
                (*oIt)+=corval*mult;
                ++oIt;
            }
        }
    }

    // Dividing by the total probability.
    for (int other=0; other<ncells; ++other) {
        auto curcol=output.column(other);
        const double total=totalprob[other];

        for (auto& val : curcol) {
            val/=total;
        }
    }

    return output;
    END_RCPP
}
