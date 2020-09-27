#include "Rcpp.h"
#include "utils.h"

#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cmath>

double sq_distance_to_line(const double* ref, const double* grad, const double* point, std::vector<double>& working) {
    // Calculating the vector difference from "point" to "ref".
    for (auto& w : working) { 
        w=*ref - *point;
        ++point;
        ++ref;
    }

    // Calculating the vector difference from "point" to the line, and taking its norm.
    const double scale=std::inner_product(working.begin(), working.end(), grad, 0.0);
    double dist=0;
    for (auto& w : working) {
        w -= scale * (*grad);
        ++grad;
        dist+=w*w;
    }
   
    return dist;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector adjust_shift_variance(Rcpp::NumericMatrix data1, Rcpp::NumericMatrix data2, 
    Rcpp::NumericMatrix vect, double sigma2, Rcpp::IntegerVector restrict1, Rcpp::IntegerVector restrict2) 
{
    const size_t ngenes = data1.nrow();
    if (ngenes != data2.nrow() || ngenes != vect.ncol()) { 
        throw std::runtime_error("number of genes do not match up between matrices");
    }
    
    const size_t ncells = data2.ncol();
    if (ncells != vect.nrow()) {
        throw std::runtime_error("number of cells do not match up between matrices");
    }

    check_subset_vector(restrict1, data1.ncol());
    check_subset_vector(restrict2, ncells);

    std::vector<double> working(ngenes), grad(ngenes);
    std::vector<std::pair<double, double> > distance1(restrict1.size());
    Rcpp::NumericVector output(ncells);

    // Iterating through all cells.
    for (size_t cell=0; cell < ncells; ++cell) {
        auto tmpcell_current = data2.column(cell);
        auto curcell = tmpcell_current.begin();

        // Calculating the l2 norm and adjusting to a unit vector.
        double l2norm = 0;
        auto _grad = vect.row(cell);
        std::copy(_grad.begin(), _grad.end(), grad.begin());
        for (const auto& g : grad) {
            l2norm += g * g;
        }

        l2norm=std::sqrt(l2norm);
        if (l2norm) { // Avoid division by zero.
            for (auto& g : grad) { 
                g/=l2norm;
            }
        }
                
        const double curproj=std::inner_product(grad.begin(), grad.end(), curcell, 0.0);

        // Getting the cumulative probability of each cell in its own batch.
        // We keep values logged to avoid underflow at large distances.
        double prob2=0;
        {
            double totalprob2=0;
            bool starting_prob=true, starting_total=true;

            for (size_t same : restrict2) {
                bool add_prob=true;
                double log_prob=0;

                // If same==cell, probability of 1 is always added, so log_prob is simply 0.
                if (same!=cell) { 
                    auto tmpcell_same = data2.column(same);
                    auto samecell = tmpcell_same.begin();
                    const double sameproj=std::inner_product(grad.begin(), grad.end(), samecell, 0.0); // Projection
                    const double samedist=sq_distance_to_line(curcell, grad.data(), samecell, working); // Distance.
                    
                    log_prob=-samedist/sigma2;
                    if (sameproj > curproj) {
                        add_prob=false;
                    }
                }
    
                if (add_prob) {
                    if (starting_prob) {
                        prob2=log_prob;
                        starting_prob=false;
                    } else {
                        prob2=R::logspace_add(prob2, log_prob);
                    }
                }
                if (starting_total) {
                    totalprob2 =log_prob;
                    starting_total=false;
                } else {
                    totalprob2=R::logspace_add(totalprob2, log_prob);
                }
            }
            prob2-=totalprob2;
        }

        // Filling up the coordinates and weights for the reference batch.
        double totalprob1=0;
        {
            bool starting_total=true;
            for (size_t other=0; other<restrict1.size(); ++other) {
                auto tmpcell_other = data1.column(restrict1[other]);
                auto othercell = tmpcell_other.begin();
                distance1[other].first=std::inner_product(grad.begin(), grad.end(), othercell, 0.0); // Projection
                const double otherdist=sq_distance_to_line(curcell, grad.data(), othercell, working); // Distance.
                distance1[other].second=-otherdist/sigma2;

                if (starting_total) {
                    totalprob1=distance1[other].second;
                    starting_total=false;
                } else {
                    totalprob1=R::logspace_add(totalprob1, distance1[other].second);
                }
            }
        
            std::sort(distance1.begin(), distance1.end());
        }

        // Choosing the quantile in the projected reference coordinates that matches the cumulative probability in its own batch.
        double ref_quan=R_NaReal;

        if (distance1.size()) {
            const double target=prob2+totalprob1;
            double cumulative=0;
            bool starting_cum=true;
            ref_quan=distance1.back().first;

            for (const auto& val : distance1) {
                if (starting_cum) {
                    cumulative=val.second;
                    starting_cum=false;
                } else {
                    cumulative=R::logspace_add(cumulative, val.second);
                }
                if (cumulative >= target) { 
                    ref_quan=val.first;
                    break;
                }
            }
        }

        // Distance between quantiles represents the scaling of the original vector.        
        output[cell]=(ref_quan - curproj)/l2norm;
    }
    
    return output;
}

