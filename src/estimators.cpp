// ----------------------------------------------------------------------------
// Spatial Covariance Estiamtors
// Author: Antonio Fioravanti
// Date: 2025-01-31
// ----------------------------------------------------------------------------

// Header for integration between R and C++
#include <Rcpp.h>
#include <algorithm> // for std::reverse
#include <cmath>      // for std::ceil

// Use the Rcpp namespace to avoid prefixing Rcpp classes
using namespace Rcpp;


int modulo(int a, int b) {
    int result = ((a % b) +b) % b;
    return result;
}
//-----------------------------------------------------------------------------
// M FUNCTIONS
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
int m0_fun_cpp(int k, int d = 1){
    // M0 function
    return modulo(k - 1, d) +1;
}

// [[Rcpp::export]]
int m_fun_cpp(int k, int s, int d, NumericVector nvec){
    // Parameters:
    // k - value
    // s - spatial dimension counter (1,...,g)
    // d - dimension of a potentially multivariate process (default should be 1)
    // nvec - vector containing the sizes of the spatial grid

    if (s == 1){
        // Calculate the ceiling of k/d and apply the modulus operation
        int result = modulo((int)std::ceil((double)k/d) - 1,
                            d * (int)nvec[0]) + 1;
        return result;
    } else {
        // Calculate the product of the previous dimensions nvec[1:(s-1)]
        int prod_previous_n = 1;
        for (int idx = 0; idx < (s - 1); ++idx) {
            prod_previous_n *= nvec[idx];
        }
        // Ceiling and modulus operations for the case s > 1
        int result = modulo((int)std::ceil((double)k /(d*prod_previous_n)) - 1, 
                     d * (int)(prod_previous_n * nvec[s - 1])) + 1;
        return result;
    }
}

// [[Rcpp::export]]
NumericVector compute_m0_values_cpp(int i, int j, int d = 1){
    // Call m0_fun_cpp to compute m0 values for i and j
    int m0_i = m0_fun_cpp(i,d);
    int m0_j = m0_fun_cpp(j, d);

    // Create NumericVector to store the results
    NumericVector m0_ij = NumericVector::create(m0_i, m0_j);

    // Return the vector
    return m0_ij;
}

// [[Rcpp::export]]
NumericVector compute_m_values_cpp(int i, int j, NumericVector nvec, 
                                   int d =1,
                                   bool flipsign = true, 
                                   bool flipposition = false){
    
    // Determine the number of spatial dimension
    int g = nvec.size();

    // Initialize out vector
    NumericVector m_ij(g);

    // Compute the differences using m_fun_cpp
    for (int s = 1; s <= g; s++){
        m_ij[s - 1] = m_fun_cpp(i, s, d, nvec) - m_fun_cpp(j, s, d, nvec);
    }
    // If flipsign is true, multiply each element by -1
    if (flipsign){
        for (int idx = 0; idx < g; ++idx){
            m_ij[idx] = -m_ij[idx];
        }
    }

    // If flip position is true, reverse the vector
    if (flipposition){
        std::reverse(m_ij.begin(), m_ij.end());
    }

    return m_ij;
}
