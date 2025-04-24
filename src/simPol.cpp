#include <vector>
#include <random>
#include <list>
#include <chrono>
#include <algorithm>
#include <map>
#include <getopt.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sys/stat.h>
#include <omp.h>
#include <Rcpp.h>
#include "simPol.h"

using namespace std;
using namespace Rcpp;

template <typename T>
void PrintPositionMatrixToCSV(vector<T> &matrix, int nrows, int ncols, string file_name)
{
    ofstream out(file_name);

    for (int i = 0; i < nrows; i++)
    {
        vector<int> *sites = &matrix[i];
        int site_idx = 0;
        for (int j = 0; j < ncols; j++)
        {
            if (j == (*sites)[site_idx])
            {
                out << "1,";
                site_idx++;
            }
            else
            {
                out << "0,";
            }
        }
        out << '\n';
    }
}

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
    random_device rd; // Get seed for random number generator
    default_random_engine generator;
    normal_distribution<double> distribution(mean, stddev);
    generator.seed(rd());
    vector<double> random_values;
    while(random_values.size() < length)
    {
        double number = distribution(generator);
        if (number >= min && number <= max) {
            if(round_result)
            {
                number = round(number);
            }
            random_values.push_back(number * multiplication_factor);
        }
    }

    return random_values;
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List simulate_polymerase_cpp(int k, double ksd, int k_min, int k_max, int gene_len,
                                 double alpha, double beta, double zeta, double zeta_sd,
                                 double zeta_min, double zeta_max, int cell_num,
                                 int pol_size, int add_space, double time, int steps_to_record) {

    int steric_hindrance = pol_size + add_space;
    double delta_t = 1e-4; 
    const int total_sites = gene_len + 1;
    double steps = time / delta_t;

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
    if(true)//zeta_vec_str == "")
    {
        zv = NormalDistrubtionGenerator(zeta, zeta_sd, zeta_min, zeta_max, total_sites, false, delta_t);
    }
    else {
        /*std::ifstream data(zeta_vec_str);
        if(data.is_open())
        {
            std::string line;
            while(getline(data, line, '\n')) {
                zv.emplace_back(std::stod(line));
            }
        }
        if((int)zv.size() >= total_sites)
        {
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
        std::transform(zv.begin(), zv.end(), zv.begin(), [&transform_val](auto& c){return c*transform_val;});*/
    }

    std::random_device rd; // Get seed for random number generator
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distrib(0.0, 1.0);

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<int>>> pos_matrices; 

    for (int step = 0; step < steps; step++)
    {
#pragma omp parallel for
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
                double draw = distrib(gen);
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
        bool record_to_pos_matrix = step >= (steps - steps_to_record);
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

    // Convert position matrix to a format that can be returned to R
    // Now with sites as rows and cells as columns
    Rcpp::IntegerMatrix pos_matrix_r(total_sites, cell_num);
    for (int cell = 0; cell < cell_num; cell++) {
        std::vector<int> *sites = &pos_matrix[cell];
        for (size_t j = 0; j < sites->size(); j++) {
            pos_matrix_r((*sites)[j], cell) = 1;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("probability_vector") = zv,
        Rcpp::Named("pause_sites") = y,
        Rcpp::Named("combined_cells_data") = res_all,
        Rcpp::Named("position_matrix") = pos_matrix_r
    );
}