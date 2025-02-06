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
//-----------------------------------------------------------------------------
// GENERAL FUNCTIONS
int modulo(int a, int b) {
    int result = ((a % b) + b) % b;
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
                                   bool flipposition = true){
    
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

// [[Rcpp::export]]
NumericMatrix compute_M_matrix_cpp(int N, 
                                    NumericVector nvec,
                                    int d = 1,
                                    bool flipsign = true,
                                    bool flipposition = true){
    // Spatial dimension                                    
    int g = nvec.size();
    // Compute the M matrix
    NumericMatrix M_ij_matrix(N * N, g);

    int idx = 0;
    for (int i = 1; i <= N; i++){
        for (int j = 1; j <= N; j++){
            // Compute the m values for i and j
            M_ij_matrix(idx, _) = compute_m_values_cpp(i, j,
                                                nvec, d,
                                                flipsign, flipposition);
            idx += 1;

        }
    }   
    return M_ij_matrix;    
}


// [[Rcpp::export]]
List get_submatrices_cpp(NumericMatrix X, NumericVector hvec) {

    int h1 = hvec[0]; // lag1
    int h2 = hvec[1]; // lag2
    // Function to extract submatrices from a matrix X
    int n1 = X.nrow();
    int n2 = X.ncol();

    // Finding the regions on which to shift the lags for the autocorrelation
    // Horizontal lag
    int T11 = std::max(1, 1 - h1);
    int T12 = std::min(n1, n1 - h1);

    // Vertical lag
    int T21 = std::max(1, 1 - h2);
    int T22 = std::min(n2, n2 - h2);

    // Define the submatrices
    NumericMatrix X1 = X(Range(T21 - 1, T22 - 1), Range(T11 - 1, T12 - 1));
    NumericMatrix X2 = X(Range(T21 + h2 - 1, T22 + h2 - 1),
                         Range(T11 + h1 - 1, T12 + h1 - 1));

    // Return the submatrices as a list
    return List::create(Named("X1") = X1, Named("X2") = X2);
}

//-----------------------------------------------------------------------------
// FLAT TOP TAPERS
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
double flat_top_taper_cpp(NumericVector xvec,
                          double c = 1.0,
                          double l = 1.0,
                          double b = 0.5,
                          std::string type = "rectangular"){
    // Parameters:
    // xvec - vector of lags
    // c - taper parameter
    // l - taper bandwidth
    // b - differential taper parameter only
    // type - taper type

    // Rectangular Taper 
    auto rectangular_taper = [&](const NumericVector& xvec,
                                 double c) -> double {
        double kappa = 1.0;
        for (int idx = 0; idx < xvec.size(); ++idx){
            //Check each element of the vector xvec and see if
            // its absolute value is smaller or equal to c
            // if it is, then return 1, otherwise return 0
            if (std::abs(xvec[idx]) <= c){
                kappa *= 1.0;
            } else {
                kappa *= 0.0;
            }
        }
        // Return the product of the taper elements, such that if there
        // is at least one element that is not within the taper, the
        // product will be zero.
        return kappa;
    };

    // Trapezoid Taper
    auto trapezoid_taper = [&](const NumericVector& xvec,
                               double c) -> double {
        NumericVector vals(xvec.size());
        for (int idx = 0; idx < xvec.size(); ++idx){
           vals[idx] = std::max(0.0, std::min(1.0, c - std::abs(xvec[idx])));
        }
        double kappa = std::accumulate(vals.begin(), vals.end(),
                               1.0, std::multiplies<double>());
        return kappa;
    };


    // Differentiable Taper
    auto differentiable_taper = [&](const NumericVector& xvec,
                                    double c,
                                    double b) -> double {
        double kappa = 1.0;
        for (int i = 0; i < xvec.size(); i++){

            double abs_s = std::abs(xvec[i]);
            // If |s| <= c, then kappa = 1
            if (abs_s<= c){
                kappa *= 1.0;
            // If |s| >= 1, then kappa = 0
            } else if (abs_s>=1.0) {
                kappa *= 0.0;
            } else {
                // Compute: exp( (-b * exp(-b/((s-c)^2))) / ((s-1)^2) )
                double numerator = -b * std::exp(-b / std::pow(abs_s - c, 2));
                double denominator = std::pow(abs_s - 1, 2);
                kappa *= std::exp(numerator / denominator);
            }
        }
        return kappa;
    };

    NumericVector xvecl = xvec / l;
    double kappa;

    if (type == "rectangular"){
        kappa = rectangular_taper(xvecl, c);
    } else if (type == "trapezoid"){
        kappa = trapezoid_taper(xvecl, c);
    } else if (type == "differentiable"){
        kappa = differentiable_taper(xvecl, c, b);
    } else {
        Rcpp::stop("Invalid taper type");
    }

    return kappa;
}

// [[Rcpp::export]]
NumericVector compute_taper_vector_cpp(NumericMatrix M_ij_matrix, 
                                         double c, 
                                         double l,
                                         double b = 0.5,
                                         std::string type = "rectangular") {
  int n = M_ij_matrix.nrow();            // number of rows in the matrix
  NumericVector kappa_ij(n);             // output vector
  
  // Loop over each row of the matrix
  for (int i = 0; i < n; i++) {
    // Extract the i-th row as a NumericVector.
    // Note: M_ij_matrix(i, _) returns the entire row.
    NumericVector row = M_ij_matrix(i, _);
    
    // Compute the taper for this row using your flat_top_taper_cpp function.
    kappa_ij[i] = flat_top_taper_cpp(row, c = c, l = l, b = b, type = type);
  }
  
  return kappa_ij;
}

//-----------------------------------------------------------------------------
// AUTOCOVARIANCE ESTIMATORS
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
double SpatialAutoCov_cpp(NumericMatrix X, NumericVector hvec){

    // Estimator for the Autocovariance function (page 2 of the manuscript)
    // Inputs:
    //  X - matrix of size n1 x n2
    // hvec - vector of size 2 containing the spatial lags

    // hvec is the vector of spatial lags where 
    // h1 is the horizontal lag
    // h2 is the vertical lag

    int h1 = (int)hvec[0];
    int h2 = (int)hvec[1];

    // Spatial dimension sizes
    int n1 = X.nrow();
    int n2 = X.ncol();
    int N = n1 * n2;

    // Mean of the matrix X
    double X_mean = Rcpp::mean(X);

    // Determine the index bounds

    // Horizontal lag
    int T11 = std::max(1, 1 - h1);
    int T12 = std::min(n1, n1 - h1);
    // Vertical lag
    int T21 = std::max(1, 1 - h2);
    int T22 = std::min(n2, n2 - h2);

    // Define the submatrices
    NumericMatrix X1 = X(Range(T21 - 1, T22 - 1), Range(T11 - 1, T12 - 1));
    NumericMatrix X2 = X(Range(T21 + h2 - 1, T22 + h2 - 1),
                         Range(T11 + h1 - 1, T12 + h1 - 1));

   // Convert sub-matrices to vectors
    NumericVector sub_vector1 = as<NumericVector>(X1);
    NumericVector sub_vector2 = as<NumericVector>(X2);

    // Subtract the mean from the sub-vectors
    NumericVector subvector_1 = sub_vector1 - X_mean;
    NumericVector subvector_2 = sub_vector2 - X_mean;

    // Compute the autocovariance
    double autocov = sum(subvector_1 * subvector_2) / N;

    return autocov;
}


// [[Rcpp::export]]
double SpatialAutoCov_cpp_loop(NumericMatrix X, NumericVector hvec) {
    // Estimator for the Autocovariance function (page 2 of the manuscript)
    // Inputs:
    //  X - matrix of spatial / spatiotemporal observations
    //  hvec - vector containing the dimensional lags

    int h1 = hvec[0]; // lag1
    int h2 = hvec[1]; // lag2

    // Mean of the matrix X
    double X_mean = Rcpp::mean(X);

    int n1 = X.nrow();
    int n2 = X.ncol();
    int N = n1 * n2;

    // Finding the regions on which to shift the lags for the autocorrelation
    // Horizontal lag
    int T11 = std::max(1, 1 - h1);
    int T12 = std::min(n1, n1 - h1);

    // Vertical lag
    int T21 = std::max(1, 1 - h2);
    int T22 = std::min(n2, n2 - h2);

    // Define the submatrices
    NumericMatrix X1 = X(Range(T21 - 1, T22 - 1), Range(T11 - 1, T12 - 1));
    NumericMatrix X2 = X(Range(T21 + h2 - 1, T22 + h2 - 1),
                         Range(T11 + h1 - 1, T12 + h1 - 1));
    // Convert sub-matrices to vectors
    NumericVector sub_vector1 = as<NumericVector>(X1);
    NumericVector sub_vector2 = as<NumericVector>(X2);

    // Subtract the mean from the sub-vectors
    NumericVector subvector_1 = sub_vector1 - X_mean;
    NumericVector subvector_2 = sub_vector2 - X_mean;

    // Compute the autocovariance using a for loop
    double autocov = 0.0;
    for (int i = 0; i < subvector_1.size(); ++i) {
        autocov += subvector_1[i] * subvector_2[i];
    }
    autocov /= N;

    return autocov;
}

// TODO: Implement Kronecker Product Estimator