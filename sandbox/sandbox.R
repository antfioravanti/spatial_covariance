#-------------------------------------------------------------------------------
### SANDBOX ###
#-------------------------------------------------------------------------------
options(scipen = 5) # some decimals, set to 999 for max decimals
#options(scipen = 0) # scientific display
if (!require(sp)) install.packages("sp"); library(sp)
if (!require(tidyr)) install.packages("tidyr"); library(tidyr)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(gstat)) install.packages("gstat"); library(gstat)
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(fields)) install.packages("fields"); library(fields)
if (!require(rstudioapi)) install.packages("rstudioapi"); library(rstudioapi)
if (!require(pracma)) install.packages("pracma"); library(pracma)
if (!require(parallel)) install.packages("parallel"); library(parallel)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(ncf)) install.packages("ncf"); library(ncf)
if (!require(openxlsx)) install.packages("openxlsx"); library(openxlsx)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),".."))
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions
set.seed(42)
#-------------------------------------------------------------------------------
# MODIFIED EXPONENTIAL

# GENERATING GRID
# Define the dimensions of the 2D grid
n1 = 5 # Dimension 1 size 
n2 = 5 # Dimension 2 size
nvec = c(n1, n2)
N = prod(nvec) 
# Generate spatial locations/coordinates of integers
grid = expand.grid(t1 = 1:n1, t2 = 1:n2)

# Smoothness parameter:
#   high makes the spatial process smooth
#   low makes the spatial process rougher
alpha1 = 2
alpha2 = 2
# Range parameter:
#   high makes the covariance decay more slowly with distance
#   low makes the covariance decay more quickly and with distance
lambda1 = 5 
lambda2 = 5
# Variance parameter: 
#   controls the magnitude of the overall variance
sigma = 1

# Separability parameter:
#   0 separable
#   1 non-separable
beta = 0
res = ModifiedExponentialCovariance(grid,
                                    sigma = sigma,
                                    alpha1 = alpha1,
                                    alpha2 = alpha2,
                                    lambda1 = lambda1,
                                    lambda2 = lambda2,
                                    beta = beta,
                                    test_sep = T)
true_cov = res$covariance
res$norm_diff

mvn_sim = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_cov) 
# Checking that the simulated spatial data "falls" in the correct coordinates
#grid$sim = mvn_sim 

# NOTE: the process is structured from top left instead of from bottom left
X = matrix(mvn_sim, nrow = n1, ncol = n2, byrow = T)
plot_matrix(X)
#-------------------------------------------------------------------------------
# DEMONSTRATION AUTOCOVARIANCE FUNCTION

hvec = c(1,2)
# 
# SpatialAutoCov = function(X, hvec){
#   # Estimator for the Autocovariance function (page 2 of the manuscript)
#   # Input:
#   #   X     matrix of spatial / spatiotemporal observations
#   #   hvec  vector containing the dimensional lags
#   
  h1 = hvec[1] # lag1
  h2 = hvec[2] # lag2
#   
  X_mean = mean(as.vector(X))
  n1 = dim(X)[1]
  n2 = dim(X)[2]
  N = n1*n2
#   
#   #Finding the regions on which to shift the lags for the autocorrelation
#   # Horizontal lag
  T11 = max(1, 1-h1)
  T12 = min(n1, n1-h1)

  # Vertical lag
  T21 = max(1, 1-h2)
  T22 = min(n2, n2-h2)

  c(T11, T12, T21, T22)
  c(T11+h1, T12+h1, T21+h2, T22+h2)
  X
  X[T21:T22, T11:T12]
  X[(T21+h2):(T22+h2), (T11+h1):(T12+h1)]
#   
#   #Original way but we get vectors here
#   # cov_lags = 1/N* sum((X[T11:T12, T21:T22] - X_mean) %*%
#   #                       t(X[(T11+h1):(T12+h1), (T21+h2):(T22+h2)] - X_mean))
#   
#   sub_matrix1 = X[(T11+h1):(T12+h1), (T21+h2):(T22+h2)] - X_mean
#   sub_matrix2 = X[T11:T12, T21:T22] - X_mean
#   
#   sub_vector1 = as.vector(sub_matrix1) # column wise vectorization
#   sub_vector2 = as.vector(sub_matrix2) # column wise vectorization
#   
#   C_h_hat = 1/N* sum(t((sub_vector1) * t((sub_vector2))))
#   
#   return(C_h_hat)
# }
# 
# 
# sigma^2 * exp(-hvec[1]/decay) * exp(-hvec[2]/decay)


#-------------------------------------------------------------------------------
# DEMONSTRATION M FUNCTION

Xvec = as.vector(X) # column wise vectorization
i = 1 # Position i in the vectorized X
j = 6 # Position j in the vectorized X
c(m_fun(i, s=1, d=1, nvec) - m_fun(j, s=1, d=1, nvec),
  m_fun(i, s=2, d=1, nvec) - m_fun(j, s=2, d=1, nvec))
# OPINION: the correct way is to keep positive lags for movements that go
# to the right and below. This presumes origin at point 1,1 being top left
# corner.
# In addition, consider that we need to flip the position of the generated lags
# such that:
#   first number = horizontal shift (x-axis)
#   second number = vertical shift (y-axis)
compute_m_values(i,j,nvec, d=1, flipsign = T, flipposition = T)

#-------------------------------------------------------------------------------
# DEMONSTRATION FLAT TOP TAPERS
Xvec = as.vector(X) # column wise vectorization
i = 1 # Position i in the vectorized X
j = 4 # Position j in the vectorized X

# Taper parameters
c = 1
l = 1
type = "rectangular"
m_ij = compute_m_values(i,j,nvec, d=1, flipsign = T, flipposition = T)
flat_top_taper(m_ij, c=c, l=l,type=type)
# ----

# Computing M vector for correct lags
M_ij_list = vector("list", N*N)
for(i in 1:N){
  for(j in 1:N){
    M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec,
                                                flipsign = TRUE,
                                                flipposition = TRUE)
  }
}

# Vectorize m_ij computation by converting list to matrix
M_ij_matrix = do.call(rbind, M_ij_list)
# ------------------------------------------------------------------------------
# DEMONSTRATION OF THE EMPIRICAL BANDWIDTH SELCTION RULE
grid_size = 5
n1 = grid_size
n2 = grid_size
nvec = c(grid_size, grid_size)
N = prod(nvec)
alpha = 1
lambda = 5
beta = 0
# Taper parameters
type = "rectangular"
c = 1

params = list(sigma,
              alpha1 = alpha,
              alpha2 = alpha,
              lambda1 = lambda,
              lambda2 = lambda,
              beta = beta,
              test_sep= F)

spatial_process = simulate_spatial_process(
  covariance_function = ModifiedExponentialCovariance,
  grid_size = grid_size,
  params = params,
  seed = 42)
# Spatial Process
X = spatial_process$X
# True Covariance
true_cov = spatial_process$covariance

# Plotting the Simulated Matrix
plot_matrix(X)
# Plotting the respective True Covariance
plot_matrix(true_cov)


# Computing M vector for correct lags
M_ij_list = vector("list", N*N)
for(i in 1:N){
  for(j in 1:N){
    M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec,
                                                flipsign = TRUE,
                                                flipposition = TRUE)
  }
}

# Vectorize m_ij computation by converting list to matrix
M_ij_matrix = do.call(rbind, M_ij_list)

# Compute all autocovariances at once using vectorized form
# Use sapply to map through each row of the m_ij matrix
C_ij_vector = sapply(1:nrow(M_ij_matrix), 
                     function(idx) SpatialAutoCov(X, M_ij_matrix[idx, ]))

# Now reshape it back to a matrix form (filled by columns)
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute Autocovariance lag 0,0
C_00 = SpatialAutoCov(X, c(0,0))

# Compute regular GammaEst matrix normalized with C_00
GammaEst = C_ij_matrix

RhoEst = GammaEst/C_00
plot_matrix(GammaEst)
# Empirical rule for selecting bandwidth
select_bandwidth_empirical = function(GammaEst,
                                       N,
                                       C0 = 2,
                                       K_N = max(4, sqrt(log10(N))),
                                       c_eff = 0.5) {
  # Bandwidth selection based on empirical rule for flat-top kernel
  # GammaEst: estimated correlation values (matrix)
  # N: sample size
  # C0: constant threshold
  # K_T: function of sample size T
  # c_eff: effective constant for scaling
  
  n <- nrow(GammaEst)
  S <- matrix(0, n, n)
  
  for (j in 1:n) {
    for (k in j:n) {
      if (j == k) {
        # Case j = k
        q_hat = which.max(sapply(0:K_N, function(m) {
          if ((j + m) > n) {
            return(TRUE)
          }
          min(abs(GammaEst[j, j + m]) < C0 * sqrt(log10(N) / N),
            abs(GammaEst[j+m, j]) < C0 * sqrt(log10(N) / N)
          )
        }))
        S[j, k] <- max(ceiling(q_hat / c_eff), 1)
      } else {
        # Case j != k
        q_hat_jk <- which.max(sapply(0:K_N, function(m) {
          if ((k + m) > n) {
            return(TRUE)
          }
          abs(GammaEst[j, k + m]) < C0 * sqrt(log10(N) / N)
        })) - 1
        
        q_hat_kj <- which.max(sapply(0:K_N, function(m) {
          if ((j + m) > n) {
            return(TRUE)
          }
          abs(GammaEst[k, j + m]) < C0 * sqrt(log10(N) / N)
        })) - 1
        
        q_hat <- max(q_hat_jk, q_hat_kj)
        S[j, k] <- S[k, j] <- max(ceiling(q_hat / c_eff), 1)
      }
    }
  }
  return(S)
}


S = select_bandwidth_empirical(GammaEst, N)
S_vector <- as.vector(S)

# Expand S values to M_ij_matrix rows
M_ij_with_S <- cbind(M_ij_matrix, S_vector)
M_ij_with_S[, 1:2] / M_ij_with_S[, 3]

kappa_ij_vector = compute_taper_vector(M_ij_matrix, c=1, l=S, type = "rectangular")
# Compute taper matrix
kappa_ij_matrix = matrix(kappa_ij_vector, nrow=N, ncol=N)

#-------------------------------------------------------------------------------
bandwidth_selection_spatial = function(
    corrMat,   # an N x N correlation matrix, with N = n1*n2
    n1,        # number of "rows" in the original lattice
    n2,        # number of "columns" in the original lattice
    C0   = 2, 
    c_ef = 1
  ) {
  # ------------------------------------------------------------------
  # 1) Basic checks and definitions
  # ------------------------------------------------------------------
  N = nrow(corrMat)
  
  if(N != ncol(corrMat)) {
    stop("corrMat must be square")
  }
  if(N != n1*n2) {
    stop("n1*n2 must match the dimension of corrMat.")
  }
  
  # Helper to map vector index k in 1..N  ->  (row, col) in 1..n1 x 1..n2
  # using row-major ordering:
  #   row(i) = ((k-1) %/% n2) + 1
  #   col(i) = ((k-1) %%  n2) + 1
  get_row = function(k) ((k - 1) %/% n2) + 1
  get_col = function(k) ((k - 1) %%  n2) + 1
  
  # ------------------------------------------------------------------
  # 2) Extract "pure row-lag" correlations, i.e.  rho(h1, 0)
  #    We'll define a function rowCorr(h1) = average of corrMat[i,j]
  #    over all i,j s.t. row(i)-row(j)=h1, col(i)=col(j).
  #    We only need h1 in [0..(n1-1)] for stationarity reasons.
  # ------------------------------------------------------------------
  rowCorrVec = numeric(n1)  # rowCorrVec[h+1] = average correlation at lag h
  for(h in 0:(n1-1)) {
    ss = 0
    count = 0
    # loop over all pairs (i,j)
    for(i in 1:N) {
      ri = get_row(i)
      ci = get_col(i)
      # We require all the columns j such that 
      # row(i)-row(j)=h1 and col(i)=col(j) i.e. h2=0
      rj = ri - h
      cj = ci
      # check if valid
      if(rj >= 1 && rj <= n1) {
        # translate (rj, cj) back to index j
        j = (rj - 1)*n2 + cj
        # accumulate
        ss = ss + corrMat[i, j]
        count = count + 1
      }
    }
    rowCorrVec[h+1] = if(count > 0) ss / count else 0
  }
  
  # ------------------------------------------------------------------
  # 3) Extract "pure column-lag" correlations, i.e.  rho(0, h2)
  #    colCorrVec[h+1] = average correlation at lag h (0..n2-1)
  #    row(i)=row(j), col(i)-col(j)=h2
  # ------------------------------------------------------------------
  colCorrVec = numeric(n2)
  for(h in 0:(n2-1)) {
    ss = 0
    count = 0
    for(i in 1:N) {
      ri = get_row(i)
      ci = get_col(i)
      # We require all the rows j such that 
      # col(i)-col(j)=h2 and row(i)=row(j) i.e. h1=0
      rj = ri
      cj = ci - h
      if(cj >= 1 && cj <= n2) {
        j = (rj - 1)*n2 + cj
        ss = ss + corrMat[i, j]
        count = count + 1
      }
    }
    colCorrVec[h+1] = if(count > 0) ss / count else 0
  }
  
  # ------------------------------------------------------------------
  # 4) Politis-style cutoff function (1D)
  #    Finds smallest q >= 0 such that for m=0..K_T, |rho[q+m]| < threshold
  # ------------------------------------------------------------------
  find_q_1d = function(rho_vec, threshold, K_T) {
    max_lag = length(rho_vec) - 1
    for(q in 0:max_lag) {
      ok = TRUE
      for(m in 0:K_T) {
        qm = q + m
        if(qm > max_lag) break
        if(abs(rho_vec[qm + 1]) >= threshold) {
          ok = FALSE
          break
        }
      }
      if(ok) return(q)
    }
    return(max_lag)
  }
  
  # ------------------------------------------------------------------
  # 5) Apply threshold => find q1, q2
  # ------------------------------------------------------------------
  threshold_val = C0 * sqrt(log10(N) / N)
  K_T = max(5, sqrt(log10(N)))
  q1 = find_q_1d(rowCorrVec, threshold_val, K_T)
  q2 = find_q_1d(colCorrVec, threshold_val, K_T)
  
  # ------------------------------------------------------------------
  # 6) Convert q1, q2 into final product-kernel bandwidths (l1,l2)
  # ------------------------------------------------------------------
  l1 = max(ceiling(q1 / c_ef), 1)
  l2 = max(ceiling(q2 / c_ef), 1)
  
  # return
  return(list(l1 = l1, l2 = l2))
}

l_res = bandwidth_selection_spatial(GammaEst, n1 = n1, n2 = n2)
# Results:
# $l1
# [1] 9
# 
# $l2
# [1] 12

# !!!!!!!!!!!!!!!!!!!!!!!!!!!
#-------------------------------------------------------------------------------
# UNTIL HERE FOR NOW
#-------------------------------------------------------------------------------
# 2nd TRY

# automatic_bandwidth = function(M_ij_matrix,
#                                X,
#                                C0 = 2,
#                                K_N = max(5, sqrt(log10(N)))) {
#   # Bandwidth selection based on empirical rule for flat-top kernel
#   # GammaEst: estimated correlation values (matrix)
#   # N: sample size
#   # C0: constant threshold
#   # K_T: function of sample size T
#   
#   n = nrow(X)
#   
#   C_00 = SpatialAutoCov(X, c(0,0))
#   for(idx in 1:nrow(M_ij_matrix)){
#   vecotor= which.max(sapply(1:K_N, function(m) {
#              if ((M_ij_matrix[idx, 1] + m) > n |
#                   (M_ij_matrix[idx, 2] + m) > n) {
#                 return(TRUE)
#               }
#               abs(SpatialAutoCov(X, M_ij_matrix[idx, ] + m) <
#                                     C0 * sqrt(log10(N) / N)
#             })) - 1
#   }



#-------------------------------------------------------------------------------
# First try Kronecker Product
  for(i in 1:g){
    kappas_sep[[i]] = sapply(1:nrow(M_ij_matrix), function(idx){
      # Create a vector of zeros of length g
      hvec = rep(0, g)
      # Place the value of M_ij_matrix[idx, i] at the i-th position
      hvec[i] = M_ij_matrix[idx, i]
      flat_top_taper(hvec, c=c, l=l, type = type)
    })
    
    C_sep_vector = sapply(1:nrow(M_ij_matrix),function(idx){
      # Create a vector of zeros of length g
      hvec = rep(0, g)
      # Place the value of M_ij_matrix[idx, i] at the i-th position
      hvec[i] = M_ij_matrix[idx, i]
      SpatialAutoCov(X, hvec)
    })
    
    # Now reshape it back to a matrix form (filled by columns)
    cov_sep[[i]] = matrix(C_sep_vector, nrow=N, ncol=N)
    
    tapered_cov_sep[[i]] = kappas_sep[[i]] * cov_sep[[i]]
  }
  
  sep_gamma_est = tapered_cov_sep[[g]]
  for(i in (g-1):1){
    sep_gamma_est = kronecker(sep_gamma_est, tapered_cov_sep[[i]])
  }

  
  # 3rd try
  for(i in 1:g){
    kappas_sep[[i]] = sapply(1:length(m_ij_vector_n), function(idx){
      # Create a vector of zeros of length g
      hvec = rep(0, g)
      # Place the value of m_ij_vector_n[idx] at the i-th position
      hvec[i] = m_ij_vector_n[idx]
      flat_top_taper(hvec, c=c, l=l, type = type)
    })
    
    C_sep_vector = sapply(1:length(m_ij_vector_n),function(idx){
      # Create a vector of zeros of length g
      hvec = rep(0, g)
      # Place the value of m_ij_vector_n[idx] at the i-th position
      hvec[i] = m_ij_vector_n[idx]
      SpatialAutoCov(X, hvec)
    })
    
    # Now reshape it back to a matrix form (filled by columns)
    cov_sep[[i]] = matrix(C_sep_vector, nrow=n1, ncol=n1)
    
    tapered_cov_sep[[i]] = kappas_sep[[i]] * cov_sep[[i]]
  }
  
  sep_gamma_est = tapered_cov_sep[[g]]
  for(i in (g-1):1){
    sep_gamma_est = kronecker(sep_gamma_est, tapered_cov_sep[[i]])
  }
  #-----------------------------------------------------------------------------
  
  
  
  #-------------------------------------------------------------------------------
# 
# GammaEst = matrix(NA, nrow = N, ncol = N)
# GammaEstTapered = matrix(NA, nrow = N, ncol = N)
# C_00 = SpatialAutoCov_v2(X, c(0,0)) #  Autocovariance at lags 0, 0
# 
# 
# 
# # FASTER VERSION BELOW
# tic("Estimation Autocovariances")
# for(i in 1:N){
#   for(j in 1:N){
#     # Remember that we should vectorize theoretically row wise
#     # so that the direction in which i and j are running is along the rows
#     m_ij = compute_m_values(i,j, nvec, flipsign = T, flipposition = T)
# 
#     C_ij = SpatialAutoCov_v2(X, m_ij)
#     # Regular Autocovariance estimator
#     # (divided by Autocovariance at lags 0,0)
#     GammaEst[i,j] =  C_ij / C_00
# 
#     # Tapered Autocovariance estimator
#     # (divided by Autocovariance at lags 0,0)
#     # GammaEstTapered[i,j] = flat_top_taper(m_ij, c=c, l=l, type=type) *
#     #   C_ij / C_00
#   }
# }
# toc()
# norm(true_cov - GammaEst, type = "F")
# 
# # is_positive_semi_definite(GammaEst)
# # is_positive_semi_definite(GammaEstTapered)
# 
# #norm(GammaEst- GammaEstTapered, type = "2")
# 
# #norm(true_cov - GammaEstTapered, type = "2")
# 
# plot_matrix(true_cov, main = "True Covariance", labels = F)
# plot_matrix(GammaEst, main = "Estimated Covariance", labels = F)



# Define the non-separable covariance function
non_separable_cov <- function(h1, h2, sigma2 = 1, lambda1 = 1, lambda2 = 1, alpha1 = 1, alpha2 = 1, rho = 0.5) {
  return(sigma2 * exp(-(abs(h1)^alpha1 / lambda1 + abs(h2)^alpha2 / lambda2 + rho * abs(h1 - h2))))
}

# Define the dimensions of the grid (spatial locations)
n1 <- 10  # Number of grid points in the first dimension
n2 <- 10  # Number of grid points in the second dimension
n_total <- n1 * n2

# Generate spatial locations (grid)
grid <- expand.grid(x1 = 1:n1, x2 = 1:n2)

# Create an empty covariance matrix
cov_matrix <- matrix(0, n_total, n_total)

# Fill the covariance matrix based on the spatial distances
for (i in 1:n_total) {
  for (j in 1:n_total) {
    h1 <- grid[i, "x1"] - grid[j, "x1"]  # Difference in first dimension
    h2 <- grid[i, "x2"] - grid[j, "x2"]  # Difference in second dimension
    cov_matrix[i, j] <- non_separable_cov(h1, h2, sigma2 = 1, lambda1 = 5, lambda2 = 5, alpha1 = 1, alpha2 = 1, rho = 0.5)
  }
}


# Define the non-separable covariance function
non_separable_cov <- function(h1, h2, sigma2 = 1, lambda1 = 1, lambda2 = 1, alpha1 = 1, alpha2 = 1, rho = 0.5) {
  return(sigma2 * exp(-(abs(h1)^alpha1 / lambda1 + abs(h2)^alpha2 / lambda2 + rho * abs(h1 - h2))))
}

# Define the dimensions of the grid (spatial locations)
n1 <- 10  # Number of grid points in the first dimension
n2 <- 10  # Number of grid points in the second dimension
n_total <- n1 * n2

# Generate spatial locations (grid)
grid <- expand.grid(x1 = 1:n1, x2 = 1:n2)

# Vectorized computation of the covariance matrix
# Compute pairwise differences for each dimension
h1_diff <- outer(grid$x1, grid$x1, "-")  # Difference in first dimension
h2_diff <- outer(grid$x2, grid$x2, "-")  # Difference in second dimension

# Apply the non-separable covariance function to the differences in a vectorized way
cov_matrix_2 <- non_separable_cov(h1_diff, h2_diff, sigma2 = 1, lambda1 = 5, lambda2 = 5, alpha1 = 1, alpha2 = 1, rho = 0.5)
#-------------------------------------------------------------------------------

# Taper parameters
c = 2
l = 1
type = "trapezoid" # rectangular or trapezoid



tic()
# Vectorize the computation of m_ij
compute_all_m_values <- function(N, nvec, d = 1) {
  m_matrix <- matrix(0, nrow = N, ncol = N)
  
  grid <- expand.grid(1:N, 1:N)  # Create all pairs of (i, j)
  
  m_ij_list <- apply(grid, 1, function(idx) {
    i <- idx[1]
    j <- idx[2]
    compute_m_values(i, j, nvec, d = d, flipsign = TRUE, flipposition = FALSE)
  })
  
  return(matrix(unlist(m_ij_list), nrow = N, byrow = TRUE))  # Convert to matrix
}

compute_autocov_matrix <- function(X, hvec_list, N) {
  X_mean <- mean(as.vector(X))
  result <- matrix(NA, N, N)
  
  # Use mapply to apply autocovariance over hvec_list
  result <- mapply(function(hvec) {
    SpatialAutoCov_v2(X, hvec)
    
  }, hvec_list)
  
  return(result)
}

taper_autocov_matrix <- function(m_matrix, c, l, type = "trapezoid") {
  apply(m_matrix, 1, function(m_ij) {
    flat_top_taper(m_ij, c = c, l = l, type = type)
  })
}

GammaEst_2 <- matrix(NA, N, N)
GammaEstTapered_2 <- matrix(NA, N, N)

C_00 <- SpatialAutoCov_v2(X, c(0, 0))  # Autocovariance at lag (0, 0)

# Precompute m_ij for all i, j
m_matrix <- compute_all_m_values(N, nvec, d = 1)

# Vectorize the calculation of GammaEst and GammaEstTapered
GammaEst_2 <- compute_autocov_matrix(X, m_matrix, N) / C_00

# Vectorize the tapered estimator
tapered_values <- taper_autocov_matrix(m_matrix, c, l, type)
GammaEstTapered_2 <- tapered_values * GammaEst

toc()



m_values_list <- lapply(1:nrow(grid), function(idx) {
  i <- grid$t1[idx]
  j <- grid$t2[idx]
  compute_m_values(i, j, nvec, flipsign = TRUE, flipposition = TRUE)
})


C_ij <- sapply(m_values_list, function(m_ij) SpatialAutoCov_v2(X, m_ij))

est_mat = C_ij / C_00
est_mat = matrix(est_mat, nrow = N, ncol = N)
norm(est_mat - GammaEst, type = "2")



# Vectorize the computation of m_ij
compute_all_m_values <- function(N, grid, d = 1) {
  m_matrix <- matrix(0, nrow = N, ncol = N)
  
  
  m_ij_list <- apply(grid, 1, function(idx) {
    i <- idx[1]
    j <- idx[2]
    compute_m_values(i, j, nvec, d = d, flipsign = TRUE, flipposition = T)
  })
  
  return(matrix(unlist(m_ij_list), nrow = N, byrow = TRUE))  # Convert to matrix
}

m_matrix  = compute_all_m_values(N, grid)
compute_autocov_matrix <- function(X, hvec_list, N) {
  X_mean <- mean(as.vector(X))
  result <- matrix(NA, N, N)
  
  # Use mapply to apply autocovariance over hvec_list
  result <- mapply(function(hvec) {
    SpatialAutoCov_v2(X, hvec)
  }, hvec_list)
  
  return(result)
}
GammaEst <- compute_autocov_matrix(X, m_matrix, N) / C_00

outer(1:N, 1:N, Vectorize(function(i, j) compute_m_values(i, j, nvec,
                                                          flipsign = TRUE,
                                                          flipposition = TRUE)[1]))







norm(true_cov - GammaEst_2, type = "2")



tic()
# Precompute m_ij values for all pairs (i, j)
M_ij_list = vector("list", N*N)
for(i in 1:N){
  for(j in 1:N){
    M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec, flipsign = TRUE, flipposition = TRUE)
  }
}

# Vectorize m_ij computation by converting list to matrix
M_ij_matrix = do.call(rbind, M_ij_list)

# Compute all autocovariances at once using vectorized form
# Use lapply to map through each row of the m_ij matrix
C_ij_vector = sapply(1:nrow(M_ij_matrix), function(idx) SpatialAutoCov_v2(X,
                                                                          M_ij_matrix[idx, ]))

# Now reshape it back to a matrix form
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute regular GammaEst matrix
GammaEst_2 = C_ij_matrix / C_00
toc()

norm(GammaEst - GammaEst_2, type = "2")
norm(true_cov - GammaEst_2, type = "2")


# Apply tapering (vectorized over the whole matrix)
Taper_matrix = apply(M_ij_matrix, 1, function(m_ij) flat_top_taper(m_ij, c=c, l=l, type=type))

# Reshape Taper_matrix and compute GammaEstTapered
Taper_matrix = matrix(Taper_matrix, nrow=N, ncol=N)
GammaEstTapered = Taper_matrix * (C_ij_matrix / C_00)









# Vectorized computation of m_ij values with correct order
compute_m_values_vectorized_faster = function(grid, nvec, d=1, flipsign=TRUE, flipposition=FALSE) {
  
  # Number of spatial dimensions
  g <- length(nvec)
  
  # Function to compute m_fun for all indices at once
  m_fun_vectorized <- function(k, s, d, nvec) {
    if (s == 1) {
      m <- ((ceiling(k / d) - 1) %% (d * nvec[1])) + 1
    } else {
      prod_prev_n <- prod(nvec[1:(s-1)])
      m <- ((ceiling(k / (d * prod_prev_n)) - 1) %% (d * prod_prev_n * nvec[s])) + 1
    }
    return(m)
  }
  
  # Precompute m_i and m_j for all spatial dimensions and all pairs
  m_ij_matrix <- matrix(NA, nrow = nrow(grid), ncol = g)
  
  for (s in 1:g) {
    # Calculate m_i and m_j for each dimension `s`
    m_i_s <- m_fun_vectorized(grid$t1, s, d, nvec)
    m_j_s <- m_fun_vectorized(grid$t2, s, d, nvec)
    
    # Compute the difference for the current dimension `s`
    m_ij_matrix[, s] <- m_i_s - m_j_s
  }
  
  # Apply the flipsign and flipposition if needed
  if (flipsign) {
    m_ij_matrix <- m_ij_matrix * -1
  }
  if (flipposition) {
    m_ij_matrix <- t(apply(m_ij_matrix, 1, rev))
  }
  
  return(m_ij_matrix)
}

tic()
# Call the vectorized function to compute the M_ij_matrix
M_ij_matrix <- compute_m_values_vectorized_faster(grid,
                                                  nvec, flipsign = F, flipposition = TRUE)

# Now proceed with the remaining computation

# Compute all autocovariances at once using vectorized form
C_ij_vector = sapply(1:nrow(M_ij_matrix), function(idx) SpatialAutoCov_v2(X, M_ij_matrix[idx, ]))

# Reshape it back to a matrix form
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute regular GammaEst matrix
GammaEst_2 = C_ij_matrix / C_00
toc()


# Apply tapering (vectorized over the whole matrix)
Taper_matrix = apply(M_ij_matrix, 1, function(m_ij) flat_top_taper(m_ij, c=c, l=l, type=type))

# Reshape Taper_matrix and compute GammaEstTapered
Taper_matrix = matrix(Taper_matrix, nrow=N, ncol=N)
GammaEstTapered = Taper_matrix * (C_ij_matrix / C_00)


#-------------------------------------------------------------------------------
sigma = 1
grid_size = 5
n1 = grid_size
n2 = grid_size
nvec = c(grid_size, grid_size)
N = prod(nvec)
alpha = 1
lambda = 5
beta = 0
# Taper parameters
type = "rectangular"
c = 1

params = list(sigma,
              alpha1 = alpha,
              alpha2 = alpha,
              lambda1 = lambda,
              lambda2 = lambda,
              beta = beta,
              test_sep= F)

spatial_process = simulate_spatial_process(
  covariance_function = ModifiedExponentialCovariance,
  grid_size = grid_size,
  params = params,
  seed = 42)
# Spatial Process
X = spatial_process$X
# True Covariance
true_cov = spatial_process$covariance

# Plotting the Simulated Matrix
plot_matrix(X, main = "Spatial Process")
# Plotting the respective True Covariance
plot_matrix(true_cov, main = "True Covariance")


# 2) Computing M Matrix for lags corrected with the m functions
M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)

# 3) Compute all autocovariances
C_ij_vector = compute_autocovariance_vector_cpp(X, M_ij_matrix)

# Reshape back to a Matrix
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N, byrow=T)

# Compute Autocovariance lag 0,0
C_00 = SpatialAutoCov_cpp(X, c(0,0))

# Naive Autocovariance Estimator
GammaEst = C_ij_matrix

plot_matrix(GammaEst, main = "Estimated Covariance")
# Naive Autocorrelation Estimator
RhoEst = C_ij_matrix / C_00
plot_matrix(RhoEst, main = "Estimated Correlation")


C_00 = SpatialAutoCov_cpp(X, c(0,0))
n1 = dim(X)[2]
n2 = dim(X)[1]

rowRhos <- numeric(2*(n2-1) + 1) 
for(lag in -(n2-1):(n2-1)){
  rowRhos[lag + n2] = SpatialAutoCov_cpp(X, c(lag,0) ) / C_00

}


colRhos <- numeric(2*(n1-1) + 1) 
for(lag in -(n1-1):(n1-1)){
  colRhos[lag + n1] = SpatialAutoCov_cpp(X, c(0,lag) ) / C_00
  
}


rowRhos_sorted = sort(unique(rowRhos), decreasing = T)
colRhos_sorted = sort(unique(colRhos), decreasing = T)


find_q_1d = function(rho_vec, threshold, K_T) {
  max_lag = length(rho_vec) - 1
  for(q in 0:max_lag) {
    ok = TRUE
    for(m in 0:K_T) {
      qm = q + m
      if(qm > max_lag) break
      if(abs(rho_vec[qm + 1]) >= threshold) {
        ok = FALSE
        break
      }
    }
    if(ok) return(q)
  }
  return(max_lag)
}




# ------------------------------------------------------------------
# 5) Apply threshold => find q1, q2
# ------------------------------------------------------------------
threshold_val = C0 * sqrt(log10(N) / N)
K_T = max(5, sqrt(log10(N)))
q1 = find_q_1d(rowRhos_sorted, threshold_val, K_T)
q2 = find_q_1d(colRhos_sorted, threshold_val, K_T)

# ------------------------------------------------------------------
# 6) Convert q1, q2 into final product-kernel bandwidths (l1,l2)
# ------------------------------------------------------------------
l1 = max(ceiling(q1 / c_ef), 1)
l2 = max(ceiling(q2 / c_ef), 1)

















corrMat = RhoEst
get_row = function(k) ((k - 1) %/% n2) + 1
get_col = function(k) ((k - 1) %%  n2) + 1



rescorrVec <- data.frame(h = numeric(0),
                         ri = numeric(0),
                         ci = numeric(0),
                         rj = numeric(0),
                         cj = numeric(0),
                         i = numeric(0),
                         j = numeric(0),
                         corr = numeric(0))


rowCorrVec = numeric(n1)  # rowCorrVec[h+1] = average correlation at lag h
for(h in 0:(n1-1)) {
  ss = 0
  count = 0
  # loop over all pairs (i,j)
  for(i in 1:N) {
    ri = get_row(i)
    ci = get_col(i)
    # We require all the columns j such that 
    # row(i)-row(j)=h1 and col(i)=col(j) i.e. h2=0
    rj = ri - h
    cj = ci
    # check if valid
    if(rj >= 1 && rj <= n1) {
      # translate (rj, cj) back to index j
      j = (rj - 1)*n2 + cj
      # accumulate
      ss = ss + corrMat[i, j]
      count = count + 1
      
      iterdf = data.frame(h = h,
                          ri = ri,
                          ci = ci,
                          rj = rj,
                          cj = cj,
                          i = i,
                          j = j,
                          corr = corrMat[i, j])
      
      rescorrVec = rbind(rescorrVec, iterdf)
    }
  }
  rowCorrVec[h+1] = if(count > 0) ss / count else 0
}

# ------------------------------------------------------------------
# 3) Extract "pure column-lag" correlations, i.e.  rho(0, h2)
#    colCorrVec[h+1] = average correlation at lag h (0..n2-1)
#    row(i)=row(j), col(i)-col(j)=h2
# ------------------------------------------------------------------
colCorrVec = numeric(n2)
for(h in 0:(n2-1)) {
  ss = 0
  count = 0
  for(i in 1:N) {
    ri = get_row(i)
    ci = get_col(i)
    # We require all the rows j such that 
    # col(i)-col(j)=h2 and row(i)=row(j) i.e. h1=0
    rj = ri
    cj = ci - h
    if(cj >= 1 && cj <= n2) {
      j = (rj - 1)*n2 + cj
      ss = ss + corrMat[i, j]
      count = count + 1
    }
  }
  colCorrVec[h+1] = if(count > 0) ss / count else 0
}
