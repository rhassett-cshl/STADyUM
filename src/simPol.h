#ifndef SIMPOL_H
#define SIMPOL_H

#include <Rcpp.h>
#include <vector>
#include <random>
#include <cmath>
#include <string>

// Function to simulate polymerase movement
Rcpp::List simulate_polymerase_cpp(
    int k, double ksd, int k_min, int k_max, 
    int gene_len, double alpha, double beta, 
    double zeta, double zeta_sd, double zeta_min, 
    double zeta_max, int cell_num, int pol_size, 
    int add_space, double time, Rcpp::Nullable<Rcpp::NumericVector> time_points_to_record,
    Rcpp::Nullable<Rcpp::NumericVector> zeta_vec = R_NilValue);

#endif // SIMPOL_H 