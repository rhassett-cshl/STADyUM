#include <vector>
#include <list>
#include <chrono>
#include <algorithm>
#include <map>
#include <getopt.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sys/stat.h>
#include <Rcpp.h>
#include <Rinternals.h>  // For R_CheckUserInterrupt
#include <Rmath.h>      // For R's random number generators
#include "simPol.h"

using namespace std;
using namespace Rcpp;

void ConvertSiteDataToMatrix(vector<vector<int>> &input, vector<vector<int>> &output)
{
    for (size_t i = 0; i < output.size(); i++)
    {
        vector<int> *sites = &input[i];
        for (size_t j = 0; j < sites->size(); j++)
        {
            output[i][(*sites)[j]] = 1;
        }
    }
}

vector<double> NormalDistrubtionGenerator(double mean, double stddev, double min, double max, size_t length, bool round_result, double multiplication_factor=1)
{
    GetRNGstate();  // Initialize R's random number generator
    vector<double> random_values;
    while(random_values.size() < length)
    {
        double number = R::rnorm(mean, stddev);
        if (number >= min && number <= max) {
            if(round_result)
            {
                number = round(number);
            }
            random_values.push_back(number * multiplication_factor);
        }
    }
    PutRNGstate();  // Clean up R's random number generator
    return random_values;
}

// [[Rcpp::export]]
Rcpp::List simulate_polymerase_cpp(
    int k,
    double ksd,
    int k_min,
    int k_max,
    int gene_len,
    double alpha,
    double beta,
    double zeta,
    double zeta_sd,
    double zeta_min,
    double zeta_max,
    int cell_num,
    int pol_size,
    int add_space,
    double time,
    Rcpp::Nullable<Rcpp::NumericVector> time_points_to_record,
    Rcpp::Nullable<Rcpp::NumericVector> zeta_vec) {

    int steric_hindrance = pol_size + add_space;
    double delta_t = 1e-4; 
    const int total_sites = gene_len + 1;
    double steps = time / delta_t;

    // Convert time points to step indices
    std::vector<int> steps_to_record;
    std::vector<double> recorded_time_points;
    
    if (time_points_to_record.isNotNull()) {
        Rcpp::NumericVector time_vec = Rcpp::as<Rcpp::NumericVector>(time_points_to_record);
        for (int i = 0; i < time_vec.length(); i++) {
            double time_point = time_vec[i];
            if (time_point >= 0 && time_point <= time) {
                int step_idx = (int)(time_point / delta_t);
                steps_to_record.push_back(step_idx);
                recorded_time_points.push_back(time_point);
            }
        }
        
        // Sort steps to record in ascending order
        std::sort(steps_to_record.begin(), steps_to_record.end());
        std::sort(recorded_time_points.begin(), recorded_time_points.end());
    }

    /* Initialize an array to hold Pol II presence and absence*/
    std::vector<std::vector<int>> pos_matrix;

    for (int i = 0; i < cell_num; i++)
    {
        std::vector<int> sites;
        sites.push_back(0);
        pos_matrix.insert(pos_matrix.begin() + i, sites);
    }

    /* Construct a probability matrix to control RNAP movement
     * Generate pause sites located from kmin to kmax with sd = ksd
     */
    std::vector<double> y = NormalDistrubtionGenerator(k, ksd, k_min, k_max, cell_num, true);

    /* A matrix of probabilities to control transition from state to state
     * cols are cells, rows are positions
     */
    std::vector<double> zv;
    if (zeta_vec.isNull()) {
        zv = NormalDistrubtionGenerator(zeta, zeta_sd, zeta_min, zeta_max, 
                                      total_sites, false, delta_t);
    } else {
        // Convert Rcpp::NumericVector to std::vector<double>
        zv = Rcpp::as<std::vector<double>>(zeta_vec);
        
        // Ensure correct length
        if ((int)zv.size() > total_sites) {
            zv.resize(total_sites);
        }
        else if((int)zv.size() == total_sites - 1)
        {
            double mean = std::accumulate(zv.begin(), zv.end(), 0.0) / zv.size();
            zv.insert(zv.begin(), mean);
        }
        else {
            Rcpp::Rcout << "Vector for scaling zeta is too short, check total length of the vector!";
            return -1;
        }
        double transform_val = zeta * delta_t;
        std::transform(zv.begin(), zv.end(), zv.begin(), 
                      [transform_val](double c) -> double {
                          return c * transform_val;
                      });
    }

    GetRNGstate();  // Initialize R's random number generator

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<int>>> pos_matrices; 

    for (int step = 0; step < steps; step++)
    {
        // Check for user interruption every 1000 steps
        if (step % 1000 == 0) {
            R_CheckUserInterrupt();
        }
        
        for (int cell = 0; cell < cell_num; cell++)
        {
            std::vector<int> *sites = &pos_matrix[cell];
            for (size_t i = 0; i < sites->size(); i++)
            {
                /* Determine whether polymerase can move or not
                 * criteria 1, probability larger than random draw
                 * criteria 2, enough space ahead to let polymerase advance
                 */
                int site_idx = (*sites)[i];
                double prob = site_idx == 0            ? alpha * delta_t
                              : site_idx == y.at(cell) ? beta * delta_t
                                                  : zv.at(site_idx);
                double draw = R::runif(0.0, 1.0);  // Using R's uniform random number generator
                if (prob > draw)
                {
                    size_t last_polymerase = sites->size() - 1;
                    /* Check if space ahead is larger than polymerase size */
                    if (i != last_polymerase && (*sites)[i + 1] - (*sites)[i] > steric_hindrance)
                    {
                        (*sites)[i]++;
                    }
                    /* Always allow the polymerase at the end to move */
                    else if (i == last_polymerase)
                    {
                        if ((*sites)[i] + 1 < total_sites)
                        {
                            (*sites)[i]++;
                        }
                        else
                        {
                            /* Remove polymerase if past final site */
                            sites->pop_back();
                            break; // to prevent iterating past the end of the linked list
                        }
                    }
                }
            }

            /* Ensure there are always polymerases waiting to be initialized (i.e., first col is always 1) */
            if (sites->size() == 0 || (*sites)[0] != 0)
            {
                sites->insert(sites->begin(), 0);
            }
        }
        /* Record info for studying steric hindrance */
        bool record_to_pos_matrix = std::find(steps_to_record.begin(), steps_to_record.end(), step) != steps_to_record.end();
        if (record_to_pos_matrix)
        {
            pos_matrices.push_back(pos_matrix);
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    Rcpp::Rcout << "Total time used for the simulation is " << duration.count() / 60.0 << " mins.\n";

    /* Calculate the # of polymerase at each site across all cells and output the results in csv format */
    std::vector<int> res_all(total_sites, 0);
    for (int cell = 0; cell < cell_num; cell++)
    {
        std::vector<int> *sites = &pos_matrix[cell];
        for (size_t j = 0; j < sites->size(); j++)
        {
            res_all[(*sites)[j]]++;
        }
    }

    // Convert the recorded position matrices to a named list
    // Each element corresponds to a specific time point
    Rcpp::List pos_matrices_named;
    for (int step = 0; step < pos_matrices.size(); step++) {
        // Create a 2D matrix for this time point: [sites, cells]
        Rcpp::IntegerMatrix pos_matrix_2d(total_sites, cell_num);
        
        for (int cell = 0; cell < cell_num; cell++) {
            std::vector<int> *sites = &pos_matrices[step][cell];
            for (size_t j = 0; j < sites->size(); j++) {
                int site_idx = (*sites)[j];
                if (site_idx < total_sites) {
                    pos_matrix_2d(site_idx, cell) = 1;
                }
            }
        }
        
        // Name the matrix with its corresponding time point
        std::string time_name = "t_" + std::to_string(recorded_time_points[step]);
        pos_matrices_named[time_name] = pos_matrix_2d;
    }

    Rcpp::IntegerMatrix final_pos_matrix(total_sites, cell_num);
    for (int cell = 0; cell < cell_num; cell++) {
        std::vector<int> *sites = &pos_matrix[cell];
        for (size_t j = 0; j < sites->size(); j++) {
            int site_idx = (*sites)[j];
            if (site_idx < total_sites) {
                final_pos_matrix(site_idx, cell) = 1;
            }
        }
    }

    PutRNGstate();  // Clean up R's random number generator

    return Rcpp::List::create(
        Rcpp::Named("probabilityVector") = zv,
        Rcpp::Named("pauseSites") = y,
        Rcpp::Named("combinedCellsData") = res_all,
        Rcpp::Named("positionMatrices") = pos_matrices_named,
        Rcpp::Named("finalPositionMatrix") = final_pos_matrix
    );
}