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
                                   int d = 1,
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
                          double l,
                          double c,
                          std::string type,
                          double b = 0.5){
    // Parameters:
    // xvec - vector of lags
    // l - taper bandwidth
    // c - taper parameter
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
                                        double l,
                                        double c,
                                        std::string type,
                                        double b = 0.5) {
  int n = M_ij_matrix.nrow();            // number of rows in the matrix
  NumericVector kappa_ij(n);             // output vector
  
  // Loop over each row of the matrix
  for (int i = 0; i < n; i++) {
    // Extract the i-th row as a NumericVector.
    // Note: M_ij_matrix(i, _) returns the entire row.
    NumericVector row = M_ij_matrix(i, _);
    
    // Compute the taper for this row using your flat_top_taper_cpp function.
    kappa_ij[i] = flat_top_taper_cpp(row, l, c, type, b);
  }
  
  return kappa_ij;
}


// [[Rcpp::export]]
double flat_top_taper_1d_cpp(double x_scalar,
                             double l,   
                             double c,
                             std::string type) {
    // Scale the lag
    double x_scaled = x_scalar / l;
    double kappa = 0.0;
    
    if (type == "rectangular") {
      // Rectangular kernel: 1 if |x_scaled| <= c, otherwise 0
      kappa = (std::abs(x_scaled) <= c) ? 1.0 : 0.0;
    } else if (type == "trapezoid") {
      // Compute kappa based on the sign and value of x_scaled
      if (x_scaled >= 0.0) {
          if (x_scaled <= 1.0) {
              kappa = 1.0;
          } else if (x_scaled <= 2.0) {
              kappa = 2.0 - x_scaled; // positive branch: 2 - x
          } else {
              kappa = 0.0;
          }
      } else { // x_scaled is negative
          if (x_scaled >= -1.0) {
              kappa = 1.0;
          } else if (x_scaled >= -2.0) {
              kappa = 2.0 + x_scaled; // negative branch: 2 + x
          } else {
              kappa = 0.0;
          }
        }
      } else {
        Rcpp::stop("type must be 'rectangular' or 'trapezoid'.");
    }
    
    return kappa;
  }

  // [[Rcpp::export]]
double flat_top_taper_multi_cpp(NumericVector x_vec,
                                NumericVector L,
                                double c,
                                std::string type) {
    int d = x_vec.size();
    
    if (L.size() != d) {
      Rcpp::stop("Length of L must match the length of x_vec.");
    }
    
    double product_kappa = 1.0;
    for (int i = 0; i < d; i++) {
      double kappa_i = flat_top_taper_1d_cpp(x_vec[i], L[i], c, type);
      product_kappa *= kappa_i;
    }
    
    return product_kappa;
  }


  // [[Rcpp::export]]
NumericVector compute_multi_taper_vector_cpp(NumericMatrix M_ij_matrix,
                                        NumericVector L,
                                        double c,
                                        std::string type) {
    int n = M_ij_matrix.nrow();
    NumericVector out(n);

    for (int i = 0; i < n; i++) {
    // Extract the i-th row as a vector.
    NumericVector x_vec = M_ij_matrix(i, _);
    // Compute the multi-dimensional taper using flat_top_taper_multi.
    out[i] = flat_top_taper_multi_cpp(x_vec, L, c, type);
    }

    return out;
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
// [[Rcpp::export]]
NumericVector compute_autocovariance_vector_cpp(NumericMatrix X,
                                    NumericMatrix M_ij_matrix) {
    int n = M_ij_matrix.nrow();            // number of rows in the matrix
    NumericVector C_ij(n);             // output vector

    // Loop over each row of the matrix
    for (int i = 0; i < n; i++) {
    // Extract the i-th row as a NumericVector.
    // Note: M_ij_matrix(i, _) returns the entire row.
    NumericVector hvec = M_ij_matrix(i, _);
    // Compute the autocovariance for this h vector
    C_ij[i] = SpatialAutoCov_cpp(X, hvec);
    }

    return C_ij;
}

//-----------------------------------------------------------------------------
// TODO: Implement Kronecker Product Estimator

// Kronecker Product of two matrices
// [[Rcpp::export]]
NumericMatrix kroneckerProduct_cpp(const NumericMatrix& A, 
                               const NumericMatrix& B) {

    // Get the dimensions of the input matrices
    int a_rows = A.nrow();
    int a_cols = A.ncol();
    int b_rows = B.nrow();
    int b_cols = B.ncol();

    NumericMatrix kronmat(a_rows * b_rows, a_cols * b_cols);
    for (int i = 0; i < a_rows; i++) {
        for (int j = 0; j < a_cols; j++) {
            double a_ij = A(i, j);
            for (int k = 0; k < b_rows; k++) {
                for (int l = 0; l < b_cols; l++) {
                    kronmat(i * b_rows + k, j * b_cols + l) = a_ij * B(k, l);
                }
            }
        }
    }
    return kronmat;
}

// Kronecker Estimator Part 1

// [[Rcpp::export]]
List tapered_sep_autocovariance_cpp(NumericMatrix X, 
                                    double c, 
                                    double l, 
                                    std::string type) {
  // Get dimensions. In R, X is assumed to have a dim attribute.
  IntegerVector dims = X.attr("dim");
  int g = dims.size();  // e.g., typically 2 for a matrix; can be >2 if X is an array.
  
  // In your R function, n1 = ncol(X) and n2 = nrow(X).
  int n1 = X.ncol();    // number of columns
  int n2 = X.nrow();    // number of rows

  // Compute C(0,0) = SpatialAutoCov(X, rep(0, g))
  NumericVector h0(g, 0.0);
  double C_00 = SpatialAutoCov_cpp(X, h0);

  // Prepare lists to hold the r-specific matrices
  List kappas_sep(g);       // taper matrices
  List cov_sep(g);          // normalized covariance matrices
  List tapered_cov_sep(g);  // elementwise product (tapered covariance)

  // Temporary matrices to be re-used for each spatial dimension r.
  // We use matrices of dimension (n2 rows, n1 columns) to mimic your R code.
  NumericMatrix kappa_sep_mat(n2, n1);
  NumericMatrix C_sep_mat(n2, n1);

  // Loop over each spatial coordinate (r = 0, …, g-1, corresponding to r=1,...,g in R)
  for (int r = 0; r < g; r++) {
    // Loop over the “grid.”  
    // In your original R code, the loops were over (i in 1:n1) and (j in 1:n2).
    // Here we let:
    //    i = 0,...,n1-1 (corresponding to column index)
    //    j = 0,...,n2-1 (corresponding to row index)
    // and we assign to element (j,i) in a matrix with n2 rows and n1 columns.
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        // Create a lag vector hvec of length g (all zeros)
        NumericVector hvec(g, 0.0);
        // In the original R code: hvec[r] = i - j (with R’s 1-indexing)
        // Here we mimic that by setting:
        hvec[r] = (i - j);
        // Compute the taper weight using your C++ taper function.
        kappa_sep_mat(j, i) = flat_top_taper_cpp(hvec, c, l, type, 0.5);
        // Compute the spatial autocovariance for lag vector hvec.
        C_sep_mat(j, i) = SpatialAutoCov_cpp(X, hvec);
      }
    }
    
    // Save the computed taper matrix.
    kappas_sep[r] = clone(kappa_sep_mat);
    
    // Compute cov_sep: normalize C_sep_mat by multiplying with C_00^(-(g-1))
    double factor = std::pow(C_00, -(g - 1));
    NumericMatrix cov_mat(n2, n1);
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < n1; j++) {
        cov_mat(i, j) = factor * C_sep_mat(i, j);
      }
    }
    cov_sep[r] = cov_mat;
    
    // Compute tapered_cov_sep: elementwise product of kappa_sep_mat and C_sep_mat.
    NumericMatrix tapered_mat(n2, n1);
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < n1; j++) {
        tapered_mat(i, j) = kappa_sep_mat(i, j) * C_sep_mat(i, j);
      }
    }
    tapered_cov_sep[r] = tapered_mat;
  }

  // Return the computed lists.
  // (The final Kronecker product will be computed in R.)
  return List::create(
    Named("kappas_sep") = kappas_sep,
    Named("cov_sep") = cov_sep,
    Named("tapered_cov_sep") = tapered_cov_sep
  );
}



// [[Rcpp::export]]
List tapered_sep_auto_multi_cpp(NumericMatrix X,
                                double c ,
                                NumericVector L,
                                std::string type) {

  // Check that L is a vector (it should have length equal to the number of spatial dimensions)
  if (L.size() < 1) {
    stop("L must be a vector containing bandwidths per spatial dimension.");
  }

  // Get dimensions. In R, X is assumed to have a dim attribute.
  IntegerVector dims = X.attr("dim");
  int g = dims.size(); // number of dimensions (for a matrix, g=2)

  int n1 = X.ncol(); // number of columns
  int n2 = X.nrow(); // number of rows

  // Compute C(0,0) = SpatialAutoCov(X, rep(0, g))
  NumericVector zeros(g, 0.0);
  double C_00 = SpatialAutoCov_cpp(X, zeros);

  // Prepare lists to hold the r-specific matrices
  List kappas_sep(g);       // taper matrices
  List cov_sep(g);          // covariance matrices
  List tapered_cov_sep(g);  // tapered covariance matrices

  // Temporary matrices to be re-used for each spatial dimension r.
  NumericMatrix kappa_sep_mat(n2, n1);
  NumericMatrix C_sep_mat(n2, n1);

  // Loop over each spatial coordinate 
  // (r = 0, …, g-1, corresponding to r=1,...,g in R)

  for(int r = 0; r < g; r++){
    for(int i = 0; i < n1; i++){ // In the R code, i runs from 1 to n1 (ncol)
      for(int j = 0; j < n2; j++){ // In the R code, i runs from 1 to n1 (nrow)
        // Create a lag vector hvec of length g (all zeros)
        NumericVector hvec(g, 0.0);
        // In the original R code: hvec[r] = i - j (with R’s 1-indexing)
        // Here we mimic that by setting:
        hvec[r] = (i - j);
        // Compute the taper weight 
        // NOTE: we are using the scalar version of flat_top_taper_cpp
        // because we need a differnt taper function for each
        // spatial dimension.
        kappa_sep_mat(j, i) = flat_top_taper_cpp(hvec, L[r], c, type);
        // Compute the spatial autocovariance for lag vector hvec.
        C_sep_mat(j, i) = SpatialAutoCov_cpp(X, hvec);
      }
    }
    // Save the computed taper matrix.
    kappas_sep[r] = clone(kappa_sep_mat);

    // Compute cov_sep as: C_00^( -(g-1) ) * C_sep_mat
    double factor = std::pow(C_00, -(g - 1));
    NumericMatrix cov_mat(n2, n1);
    for (int i = 0; i < n2; i++) {
      for (int j = 0; j < n1; j++) {
        cov_mat(i, j) = factor * C_sep_mat(i, j);
      }
    }
    cov_sep[r] = cov_mat;
  

  NumericMatrix tapered_mat(n2, n1);
  for(int i = 0; i < n2; i++){
    for(int j = 0; j < n1; j++){
      tapered_mat(i, j) = kappa_sep_mat(i, j) * C_sep_mat(i, j);
    }
  }
  tapered_cov_sep[r] = tapered_mat;

} 
  
// Return a list with the three components.
return List::create(
  Named("kappas_sep") = kappas_sep,
  Named("cov_sep") = cov_sep,
  Named("tapered_cov_sep") = tapered_cov_sep
  );
}

//-----------------------------------------------------------------------------
// BANDWIDTH SELECTION
//-----------------------------------------------------------------------------

// [[Rcpp::export]]
List bandwidth_selection_spatial_cpp_v1(NumericMatrix corrMat,
                                     int n1, int n2, 
                                     double C0 = 2, double c_ef = 1) {

// ----------------------------------------------------------------------------
// bandwidth_selection_spatial_cpp
//
// Inputs:
//   corrMat : an N x N correlation matrix (with N = n1*n2)
//   n1      : number of rows in the original lattice
//   n2      : number of columns in the original lattice
//   C0      : constant (default = 2)
//   c_ef    : effective constant (default = 1)
// 
// Returns a list with two elements, l1 and l2, the selected bandwidths.
// ----------------------------------------------------------------------------
  // Basic checks
  int N = corrMat.nrow();
  if (N != corrMat.ncol()) {
    stop("corrMat must be square");
  }
  if (N != n1 * n2) {
    stop("n1*n2 must match the dimension of corrMat.");
  }
  
  // Helper lambdas: convert 1-indexed k (from 1 to N) to (row, col)
  // using row-major ordering.
  auto get_row = [n2](int k) -> int {
    return ((k - 1) / n2) + 1;
  };
  auto get_col = [n2](int k) -> int {
    return ((k - 1) % n2) + 1;
  };
  
  // 2) Extract "pure row-lag" correlations: for lag h (0 to n1-1)
  std::vector<double> rowCorrVec(n1, 0.0);
  for (int h = 0; h < n1; h++) {
    double ss = 0.0;
    int count = 0;
    for (int i = 1; i <= N; i++) {
      int ri = get_row(i);
      int ci = get_col(i);
      int rj = ri - h;
      int cj = ci;
      if (rj >= 1 && rj <= n1) {
        // Convert (rj, cj) back to a 1-indexed vector index:
        int j = (rj - 1) * n2 + cj;
        ss += corrMat(i - 1, j - 1);
        count++;
      }
    }
    rowCorrVec[h] = (count > 0) ? ss / count : 0.0;
  }
  
  // 3) Extract "pure column-lag" correlations: for lag h (0 to n2-1)
  std::vector<double> colCorrVec(n2, 0.0);
  for (int h = 0; h < n2; h++) {
    double ss = 0.0;
    int count = 0;
    for (int i = 1; i <= N; i++) {
      int ri = get_row(i);
      int ci = get_col(i);
      int rj = ri;
      int cj = ci - h;
      if (cj >= 1 && cj <= n2) {
        int j = (rj - 1) * n2 + cj;
        ss += corrMat(i - 1, j - 1);
        count++;
      }
    }
    colCorrVec[h] = (count > 0) ? ss / count : 0.0;
  }
  
  // 4) Politis-style cutoff function: find smallest q such that for m = 0,..,K_T,
  //    |rho_vec[q+m]| < threshold.
  auto find_q_1d = [&](const std::vector<double>& rho_vec, double threshold, int K_T) -> int {
    int max_lag = rho_vec.size() - 1;
    for (int q = 0; q <= max_lag; q++) {
      bool ok = true;
      for (int m = 0; m <= K_T; m++) {
        int qm = q + m;
        if (qm > max_lag) break;
        if (std::abs(rho_vec[qm]) >= threshold) {
          ok = false;
          break;
        }
      }
      if (ok) return q;
    }
    return max_lag;
  };
  
  // 5) Apply threshold to find q1 and q2.
  double threshold_val = C0 * std::sqrt(std::log10(static_cast<double>(N)) / N);
  double temp = std::sqrt(std::log10(static_cast<double>(N)));
  int K_T_val = std::max(5, static_cast<int>(temp));
  
  int q1 = find_q_1d(rowCorrVec, threshold_val, K_T_val);
  int q2 = find_q_1d(colCorrVec, threshold_val, K_T_val);
  
  // 6) Convert q1 and q2 into final product-kernel bandwidths (l1, l2)
  int l1 = std::max(static_cast<int>(std::ceil(q1 / c_ef)), 1);
  int l2 = std::max(static_cast<int>(std::ceil(q2 / c_ef)), 1);
  
  return List::create(Named("l1") = l1, Named("l2") = l2);
}



// [[Rcpp::export]]
List bandwidth_selection_spatial_cpp_v2(NumericMatrix X,
                                     int n1, int n2, 
                                     double C0 = 2, double c_ef = 1) {

  // 1) Compute C_00 (the autocovariance at lag (0,0))
  NumericVector zeroLag = NumericVector::create(0.0, 0.0);
  double C_00 = SpatialAutoCov_cpp(X, zeroLag);

  // 2) Compute dim 1 autocovariances with h = (h1, 0)
  int offset_1 = n1 - 1;   // This maps lag = -(nrow-1) to index 0.
  std::vector<double> d1Rhos(2 * (n1 - 1) + 1, 0.0);
  for (int lag = -(n1 - 1); lag <= (n1 - 1); lag++) {
    int idx = lag + offset_1; // Convert lag to 0-based index.
    NumericVector lagVec = NumericVector::create(lag, 0.0);
    d1Rhos[idx] = SpatialAutoCov_cpp(X, lagVec) / C_00;
  }

  // 3) Compute dim 2 autocovariances with h = (0, h2)
  int offset_2 = n2 - 1;   // This maps lag = -(ncol-1) to index 0.
  std::vector<double> d2Rhos(2 * (n2 - 1) + 1, 0.0);
  for (int lag = -(n2 - 1); lag <= (n2 - 1); lag++) {
    int idx = lag + offset_2; // Convert lag to 0-based index.
    NumericVector lagVec = NumericVector::create(0.0, lag);
    d2Rhos[idx] = SpatialAutoCov_cpp(X, lagVec) / C_00;
  }

  
// 4) Remove duplicated entries and sort in decreasing order
  std::vector<double> d1Rhos_sorted = d1Rhos;
  std::sort(d1Rhos_sorted.begin(), d1Rhos_sorted.end(), std::greater<double>());
  d1Rhos_sorted.erase(std::unique(d1Rhos_sorted.begin(),
                      d1Rhos_sorted.end()),
                      d1Rhos_sorted.end());

  std::vector<double> d2Rhos_sorted = d2Rhos;
  std::sort(d2Rhos_sorted.begin(), d2Rhos_sorted.end(), std::greater<double>());
  d2Rhos_sorted.erase(std::unique(d2Rhos_sorted.begin(),
                      d2Rhos_sorted.end()),
                       d2Rhos_sorted.end());

  // 4) Politis-style cutoff function: find smallest q such that for m = 0,..,K_T,
  //    |rho_vec[q+m]| < threshold.
  auto find_q_1d = [&](const std::vector<double>& rho_vec,
                       double threshold, int K_T) -> int {
    int max_lag = rho_vec.size() - 1;
    for (int q = 0; q <= max_lag; q++) {
      bool ok = true;
      for (int m = 0; m <= K_T; m++) {
        int qm = q + m;
        if (qm > max_lag) break;
        if (std::abs(rho_vec[qm]) >= threshold) {
          ok = false;
          break;
        }
      }
      if (ok) return q;
    }
    return max_lag;
  };

  // 5) Find q1 and q2
  double threshold = C0 * std::sqrt(std::log10(static_cast<double>(n1 * n2)) / (n1 * n2));
  int K_T = std::max(5, static_cast<int>(std::sqrt(std::log10(static_cast<double>(n1 * n2)))));

  int q1 = find_q_1d(d1Rhos_sorted, threshold, K_T);
  int q2 = find_q_1d(d2Rhos_sorted, threshold, K_T);

  // 6) Compute l1 and l2
  int l1 = std::max(static_cast<int>(std::ceil(q1 / c_ef)), 1);
  int l2 = std::max(static_cast<int>(std::ceil(q2 / c_ef)), 1);

  return List::create(Named("l1") = l1, Named("l2") = l2);
}


