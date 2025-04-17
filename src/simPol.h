#ifndef SIMPOL_H
#define SIMPOL_H

#include <Rcpp.h>
#include <vector>
#include <random>
#include <cmath>
#include <string>

// Function to simulate polymerase movement
int simulate_polymerase(int k, int k_min, int k_max, double ksd, 
                       size_t gene_len, double alpha, double beta, 
                       double zeta, double zeta_sd, double zeta_max, 
                       double zeta_min, int total_cells, 
                       int s, int h, double time, 
                       double delta_t, int csv_steps_to_record);

#endif // SIMPOL_H 