#include "utils.h"

#include <stdexcept>
#include <sstream>

void check_subset_vector(Rcpp::IntegerVector subvec, int len) {
    for (const auto& s : subvec) {
        if (s==NA_INTEGER || s < 0 || s >= len) {
            throw std::runtime_error("subset indices out of range");
        }
    }
    return;
}
