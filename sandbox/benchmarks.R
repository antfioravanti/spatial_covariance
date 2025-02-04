#-------------------------------------------------------------------------------
### SPATIAL PROCESS SIMULATION ###
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
if (!require(Rcpp)) install.packages("Rcpp"); library(Rcpp)
if (!require(microbenchmark)) install.packages("microbenchmark"); library(microbenchmark)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),".."))
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions

sourceCpp("src/estimators.cpp")
set.seed(42)
#-------------------------------------------------------------------------------
# COMPARE M VALUES FUNCTIONS

# Test input
i = 21
j = 1
nvec = c(20,20)
N = prod(nvec)
d = 1
flipsign = TRUE
flipposition = FALSE

# Run the R function
result_r = compute_m_values(i, j, nvec, d, flipsign, flipposition)

# Run the C++ function
result_cpp = compute_m_values_cpp(i, j, nvec, d, flipsign, flipposition)

result_r
result_cpp

# Benchmark both versions
benchmark_results <- microbenchmark(
  R_version = compute_m_values(5, 3, nvec, d = 1, flipsign = TRUE, flipposition = FALSE),
  Cpp_version = compute_m_values_cpp(5, 3, nvec, d = 1, flipsign = TRUE, flipposition = FALSE),
  times = 1000
)
print(benchmark_results)

#-------------------------------------------------------------------------------
# COMPARE M MATRIX FUNCTIONS
M_ij_matrix = compute_M_matrix(N, nvec = nvec)
M_ij_matrix_cpp = compute_M_matrix_cpp(N, nvec = nvec)

all.equal(M_ij_matrix, M_ij_matrix_cpp)

# Benchmark both versions
benchmark_results = microbenchmark(
  R_version = compute_M_matrix(N, nvec, flipsign = TRUE, flipposition = FALSE),
  Cpp_version = compute_M_matrix_cpp(N, nvec=nvec, flipsign = TRUE, flipposition = FALSE),
  times = 100
)
print(benchmark_results)
#-------------------------------------------------------------------------------
# COMPARE AUTOCOVARIANCE FUNCTIONS

n1 = 5
n2 = 5
hvec = c(4,0)

h1 = hvec[1]
h2 = hvec[2]


T11 = max(1, 1-h1)
T12 = min(n1, n1-h1)

# Vertical lag
T21 = max(1, 1-h2)
T22 = min(n2, n2-h2)

c(T11, T12, T21, T22)

sigma = 1
alpha = 1
beta = 0
grid_size = 5
lambda = 3
params = list(sigma = sigma,
              alpha1 = alpha,
              alpha2 = alpha,
              lambda1 = lambda,
              lambda2 = lambda,
              beta = beta,
              test_sep = F)


spatial_process = simulate_spatial_process(
  covariance_function = ModifiedExponentialCovariance,
  grid_size = grid_size,
  params = params,
  seed = 42)



# Get Submatrices
hvec = c(0,4)

res = GetSubmatrices(spatial_process$X, hvec)
res_cpp = get_submatrices_cpp(spatial_process$X, hvec)
spatial_process$X
res$X1
res_cpp$X1
res$X2
res_cpp$X2

estCov = SpatialAutoCov(spatial_process$X, hvec)
estCov_v2 = SpatialAutoCov_v2(spatial_process$X, hvec)
estCov_cpp = SpatialAutoCov_cpp(spatial_process$X, hvec)
estCov_cpp_loop = SpatialAutoCov_cpp_loop(spatial_process$X, hvec)
estCov
estCov_v2
estCov_cpp
estCov_cpp_loop

# Benchmark both versions
benchmark_results = microbenchmark(
  R_version = SpatialAutoCov_v2(spatial_process$X, hvec),
  Cpp_version = SpatialAutoCov_cpp(spatial_process$X, hvec=hvec),
  Cpp_version_loop = SpatialAutoCov_cpp_loop(spatial_process$X,hvec=hvec),
  times = 1000
)
