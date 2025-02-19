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
sigma = 1
alpha = 1
beta = 0
grid_size = 20
n1 = grid_size
n2 = grid_size
nvec = c(grid_size, grid_size)
g = length(nvec)
N = prod(nvec)
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


X = spatial_process$X
#-------------------------------------------------------------------------------
# COMPARE M VALUES FUNCTIONS

# Test input
i = 1
j = 15
nvec = c(20,20)
N = prod(nvec)
d = 1
flipsign = TRUE
flipposition = TRUE

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

M_ij_matrix[(N*N-10):(N*N),1:2]
M_ij_matrix_cpp[(N*N-10):(N*N),1:2]

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
hvec = c(1,3)

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
estCov_v2 = SpatialAutoCov_v0(spatial_process$X, hvec)
estCov_cpp = SpatialAutoCov_cpp(spatial_process$X, hvec)
estCov_cpp_loop = SpatialAutoCov_cpp_loop(spatial_process$X, hvec)
estCov
estCov_v2
estCov_cpp
estCov_cpp_loop

# Benchmark both versions
benchmark_results = microbenchmark(
  R_version = SpatialAutoCov(spatial_process$X, hvec),
  Cpp_version = SpatialAutoCov_cpp(spatial_process$X, hvec=hvec),
  Cpp_version_loop = SpatialAutoCov_cpp_loop(spatial_process$X,hvec=hvec),
  times = 10000
)
print(benchmark_results)
#-------------------------------------------------------------------------------
# Compare the autocovariance computation for the whole M matrix

M_ij_matrix = compute_M_matrix(N, nvec = nvec)
sapplycij = function(X, M_ij_matrix){
  
  C_ij_vector = sapply(1:nrow(M_ij_matrix), 
                       function(idx) SpatialAutoCov(X, M_ij_matrix[idx, ]))
  
  return(C_ij_vector)
}

r_cij = sapplycij(X, M_ij_matrix)
cpp_ij = compute_autocovariance_vector_cpp(X, M_ij_matrix)

all.equal(r_cij, cpp_ij)

benchmark_results <- microbenchmark(
  R_version   = sapplycij(X, M_ij_matrix),
  CPP_version = compute_autocovariance_vector_cpp(X, M_ij_matrix),
  times = 100
)
print(benchmark_results)
#-----------------------------------------------------------------------
# TESTING TAPER FUNCTION
# COMPARE M MATRIX FUNCTIONS

M_ij_matrix = compute_M_matrix(N, nvec = nvec)
M_ij_matrix_cpp = compute_M_matrix_cpp(N, nvec = nvec)

type = "rectangular"
l = 1
c = 1

# Compute Taper vector
kappa_ij_vector = sapply(1:nrow(M_ij_matrix), 
                         function(idx) flat_top_taper(M_ij_matrix[idx, ],
                                                      c=c, l=l, type = type))
# Compute Taper vector
kappa_ij_vector_cpp = sapply(1:nrow(M_ij_matrix), 
                         function(idx) flat_top_taper_cpp(M_ij_matrix[idx, ],
                                                        c=c, l=l, type = type))


benchmark_results = microbenchmark(
  R_version = flat_top_taper(M_ij_matrix[2, ],
                             c=c, l=l, type = type),
  Cpp_version = flat_top_taper_cpp(M_ij_matrix[2, ],
                                   c=c, l=l, type = type),
  times = 10000
)
summary(benchmark_results)
all.equal(kappa_ij_vector, kappa_ij_vector_cpp)

#-------------------------------------------------------------------------------
# Compare the sapply taper vs C++ taper loop


sapply_taper = function(M_ij_matrix, c, l, type){
  kappa_ij_vector = sapply(1:nrow(M_ij_matrix), 
                           function(idx) flat_top_taper(M_ij_matrix[idx, ],
                                                    c=c, l=l, type = type))
  return(kappa_ij_vector)
}

type = "rectangular"
l = 1
c = 1

kappa_ij_vector = sapply_taper(M_ij_matrix, c=c, l=l, type = type)
kappa_ij_vector_cpp = compute_taper_vector_cpp(M_ij_matrix, c=c, l=l, type = type)

all.equal(kappa_ij_vector, kappa_ij_vector_cpp)


benchmark_results = microbenchmark(
  R_version = sapply_taper(M_ij_matrix, c=c, l=l, type = type),
  Cpp_version = compute_taper_vector_cpp(M_ij_matrix, c=c, l=l, type = type),
  times = 100
)
print(benchmark_results)
#-------------------------------------------------------------------------------

# Compare tapered Autocovariance Matrix


# Compare Kronecker Product

set.seed(123)
A <- matrix(runif(20 * 20), nrow = 20, ncol = 20)
B <- matrix(runif(20 * 20), nrow = 20, ncol = 20)

# Compute the Kronecker product using R's built-in function
kronecker_R <- kronecker(A, B)

# Compute the Kronecker product using the C++ function
kronecker_cpp <- kroneckerProduct(A, B)

all.equal(kronecker_R, kronecker_cpp)
benchmark_results <- microbenchmark(
  R_kronecker = kronecker(A, B),
  Cpp_kronecker = kroneckerProduct(A, B),
  times = 100  # number of repetitions; adjust as needed
)
print(benchmark_results)
#-------------------------------------------------------------------------------

compute_multi_taper_vector_R <- function(M_ij_matrix, L, c = 1, type = "rectangular") {
  sapply(seq_len(nrow(M_ij_matrix)), function(idx) {
    x_vec <- M_ij_matrix[idx, ]
    flat_top_taper_multi(x_vec, c = c, L = L, type = type)
  })
}

L <- c(1, 1)
taper_vector_R <- compute_multi_taper_vector_R(M_ij_matrix, L, c = 1, type = "rectangular")
taper_vector_CPP <- compute_multi_taper_vector_cpp(M_ij_matrix, L, c = 1, type = "rectangular")
all.equal(taper_vector_R, taper_vector_CPP)

benchmark_results <- microbenchmark(
  R_version = compute_multi_taper_vector_R(M_ij_matrix, L, c = 1, type = "rectangular"),
  CPP_version = compute_multi_taper_vector_cpp(M_ij_matrix, L, c = 1, type = "rectangular"),
  times = 100
)
print(benchmark_results)
#-------------------------------------------------------------------------------

# Compare C++ and R Tapered preparation
set.seed(123)
X <- matrix(rnorm(100), nrow = 10, ncol = 10)

# Compute results using the R-only version.
result_R <- Tapered_Sep_Autocovariance_Kron(X, c = 1, l = 1, type = "rectangular")

# Compute results using the C++-assisted version.
result_cpp <- tapered_sep_autocovariance_cpp(X, c = 1, l = 1, type = "rectangular")

# Function to compare two lists element-by-element.
compare_list_elements <- function(list1, list2, name) {
  for (i in seq_along(list1)) {
    eq <- all.equal(list1[[i]], list2[[i]])
    if (!isTRUE(eq)) {
      cat(sprintf("Mismatch in %s at index %d:\n", name, i))
      print(eq)
    } else {
      cat(sprintf("%s at index %d match.\n", name, i))
    }
  }
}

compare_list_elements(result_R$kappas_sep, result_cpp$kappas_sep, "kappas_sep")
compare_list_elements(result_R$cov_sep, result_cpp$cov_sep, "cov_sep")
compare_list_elements(result_R$tapered_cov_sep, result_cpp$tapered_cov_sep, "tapered_cov_sep")




compare_component <- function(comp_name, r_list, cpp_list) {
  # If the component is a list (e.g., "kappas_sep") then compare each element.
  if (is.list(r_list)) {
    for (i in seq_along(r_list)) {
      equal <- all.equal(r_list[[i]], cpp_list[[i]])
      if (isTRUE(equal)) {
        message(sprintf("%s element %d: MATCH", comp_name, i))
      } else {
        message(sprintf("%s element %d: MISMATCH", comp_name, i))
        print(equal)
      }
    }
  } else {
    equal <- all.equal(r_list, cpp_list)
    if (isTRUE(equal)) {
      message(sprintf("%s: MATCH", comp_name))
    } else {
      message(sprintf("%s: MISMATCH", comp_name))
      print(equal)
    }
  }
}


X <- matrix(rnorm(100), nrow = 10, ncol = 10)

result_R <- Tapered_Sep_Autocovariance_Kron(X, c = 1, l = 1, type = "rectangular")
result_v3 <- Tapered_Sep_Autocovariance_Kron_v3(X, c = 1, l = 1, type = "rectangular")


compare_component("KronTaperCov", result_R$KronTaperCov, result_v3$KronTaperCov)
compare_component("kappas_sep", result_R$kappas_sep, result_v3$kappas_sep)
compare_component("cov_sep", result_R$cov_sep, result_v3$cov_sep)
compare_component("tapered_cov_sep", result_R$tapered_cov_sep, result_v3$tapered_cov_sep)

benchmark_results <- microbenchmark(
  R_version   = Tapered_Sep_Autocovariance_Kron(X, c = 1, l = 1, type = "rectangular"),
  Cpp_version = Tapered_Sep_Autocovariance_Kron_v3(X, c = 1, l = 1, type = "rectangular"),
  times = 100
)
print(benchmark_results)
#-------------------------------------------------------------------------------
# BANDWIDTH SELECTION BENCHMARK


# Computing M vector for correct lags
M_ij_matrix = compute_M_matrix(N, nvec = nvec)

# Compute all autocovariances at once using vectorized form
# Use sapply to map through each row of the m_ij matrix
C_ij_vector = sapply(1:nrow(M_ij_matrix), 
                     function(idx) SpatialAutoCov_v0(X, M_ij_matrix[idx, ]))

# Now reshape it back to a matrix form (filled by columns)
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute Autocovariance lag 0,0
C_00 = SpatialAutoCov_v0(X, c(0,0))

# Compute regular GammaEst matrix normalized with C_00
RhoEst = C_ij_matrix / C_00

result_R <- bandwidth_selection_spatial(RhoEst, n1, n2, C0 = 2, c_ef = 1)
result_CPP <- bandwidth_selection_spatial_cpp(RhoEst, n1, n2, C0 = 2, c_ef = 1)

all.equal(result_R, result_CPP)

benchmark_results <- microbenchmark(
  R_version   = bandwidth_selection_spatial(RhoEst, n1, n2, C0 = 2, c_ef = 1),
  CPP_version = bandwidth_selection_spatial_cpp(RhoEst, n1, n2, C0 = 2, c_ef = 1),
  times = 100
)
print(benchmark_results)


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# COMPARE TAPERED SEP AUTOCOVARIANCE MULLTIPLE KRONECKER

resR = Tapered_Sep_Autocovariance_Kron_multi(X, c=1, L=c(1,1),
                type = "rectangular")

rescpp = Tapered_Sep_Autocovariance_Kron_multi_cpp(X, c=1, L=c(1,1),
                                                   type = "rectangular")

Kron_CPP <- rescpp$KronTaperCov
Kron_R <- resR$KronTaperCov
all.equal(Kron_R, Kron_CPP)

benchmark_results <- microbenchmark(
  R_version   = Tapered_Sep_Autocovariance_Kron_multi(X, c=1, L=c(1,1),
                                                      type = "rectangular"),
  
  CPP_version = Tapered_Sep_Autocovariance_Kron_multi_cpp(X, c=1, L=c(1,1),
                                                          type = "rectangular"),
  
  times = 50
)
print(benchmark_results)

#-------------------------------------------------------------------------------
simulation_og = function(grid_size, params, c, type){
  nvec = c(grid_size, grid_size)
  g = length(nvec)
  N = prod(nvec)
  spatial_process = simulate_spatial_process(
    covariance_function = ModifiedExponentialCovariance,
    grid_size = grid_size,
    params = params,
    seed = 42)
  # Spatial Process
  X = spatial_process$X
  # True Covariance
  true_cov = spatial_process$covariance
  
  # Computing M vector for correct lags
  M_ij_matrix = compute_M_matrix(N, nvec = nvec)
  
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
  
  RhoEst = C_ij_matrix / C_00
  
  # 5) Select empirical bandwidth for each dimension
  Lvals = bandwidth_selection_spatial(
    corrMat = RhoEst,
    n1 = grid_size,
    n2 = grid_size,
    C0 = 2,
    c_ef = 1)
  
  # Turn it into a vector c(l1, l2)
  L = unlist(Lvals, use.names = FALSE)
  
  # 6) Compute Taper matrix
  kappa_ij_vector = sapply(seq_len(nrow(M_ij_matrix)), function(idx) {
    x_vec = M_ij_matrix[idx, ]  
    flat_top_taper_multi(
      x_vec, 
      c    = c,
      L    = L,   # c(l1, l2)
      type = type)
  })
  
  # Turn into matrix
  kappa_ij_matrix = matrix(kappa_ij_vector, nrow = N, ncol = N)
  
  # Compute Tapered Covariance Matrix
  GammaEstTaper = kappa_ij_matrix * GammaEst
  
  # 8) Compute Separable Taper Estimator
  SepResults = Tapered_Sep_Autocovariance_Kron_multi(
    X,
    L   = L,
    c   = c,
    type= type)
  
  GammaEstTaperSep = SepResults$KronTaperCov
  
  spectral_norm = norm(true_cov - GammaEst, type = "2")
  spectral_norm_tap = norm(true_cov - GammaEstTaper, type = "2")
  spectral_norm_sepkron = norm(true_cov - GammaEstTaperSep, type = "2")
  
  return(list(true_cov = true_cov,
              GammaEst = GammaEst,
              RhoEst = RhoEst,
              taper_matrices= kappa_ij_matrix,
              GammaEstTaper = GammaEstTaper,
              GammaEstTaperSep = GammaEstTaperSep,
              spectral_norm = spectral_norm,
              spectral_norm_tap = spectral_norm_tap,
              spectral_norm_sepkron = spectral_norm_sepkron,
              L1Selected         = L[1],
              L2Selected         = L[2]))
}


simulation_cpp = function(grid_size, params, c, type ){
  nvec = c(grid_size, grid_size)
  g = length(nvec)
  N = prod(nvec)
  spatial_process = simulate_spatial_process(
    covariance_function = ModifiedExponentialCovariance,
    grid_size = grid_size,
    params = params,
    seed = 42)
  # Spatial Process
  X = spatial_process$X
  # True Covariance
  true_cov = spatial_process$covariance
  
  # Computing M vector for correct lags
  M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)
  
  # Compute all autocovariances at once using vectorized form
  # Use sapply to map through each row of the m_ij matrix
  C_ij_vector = compute_autocovariance_vector_cpp(X, M_ij_matrix)
  
  # Now reshape it back to a matrix form (filled by columns)
  C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)
  
  # Compute Autocovariance lag 0,0
  C_00 = SpatialAutoCov_cpp(X, c(0,0))
  
  # Compute regular GammaEst matrix normalized with C_00
  GammaEst = C_ij_matrix
  
  RhoEst = C_ij_matrix / C_00
  
  
  # 5) Select empirical bandwidth for each dimension
  Lvals = bandwidth_selection_spatial_cpp_v1(
    corrMat = RhoEst,
    n1 = grid_size,
    n2 = grid_size,
    C0 = 2,
    c_ef = 1)
  
  # Turn it into a vector c(l1, l2)
  L = unlist(Lvals, use.names = FALSE)
  
  # 6) Compute Taper matrix
  kappa_ij_vector = compute_multi_taper_vector_cpp(M_ij_matrix,
                                                   L = L,
                                                   c ,
                                                   type)
  
  # Turn into matrix
  kappa_ij_matrix = matrix(kappa_ij_vector, nrow = N, ncol = N)
  
  # Compute Tapered Covariance Matrix
  GammaEstTaper = kappa_ij_matrix * GammaEst
  
  # 8) Compute Separable Taper Estimator
  SepResults = Tapered_Sep_Autocovariance_Kron_multi_cpp(
    X,
    c   = c,
    L   = L,
    type= type)
  
  
  GammaEstTaperSep = SepResults$KronTaperCov
  
  spectral_norm = norm(true_cov - GammaEst, type = "2")
  spectral_norm_tap = norm(true_cov - GammaEstTaper, type = "2")
  spectral_norm_sepkron = norm(true_cov - GammaEstTaperSep, type = "2")
  
  return(list(true_cov = true_cov,
              GammaEst = GammaEst,
              RhoEst = RhoEst,
              taper_matrices= kappa_ij_matrix,
              GammaEstTaper = GammaEstTaper,
              GammaEstTaperSep = GammaEstTaperSep,
              spectral_norm = spectral_norm,
              spectral_norm_tap = spectral_norm_tap,
              spectral_norm_sepkron = spectral_norm_sepkron,
              L1Selected         = L[1],
              L2Selected         = L[2]))
}


type = "rectangular"
c = 1
grid_size = 20

# Run the original R version
tic()
result_og <- simulation_og(grid_size, params, c = c, type = type)
toc()

# Run the C++-assisted version
tic()
result_cpp <- simulation_cpp(grid_size, params, c = c, type = type)
toc()

# Run the C++-assisted version
tic()
result_cpp_2 <- simulation_cpp(grid_size, params, c = c, type = type)
toc()


components <- c("truecov_matrices", "cov_matrices", "taper_matrices",
                "taper_covariances", "RhoEst",
                "spectral_norm", "spectral_norm_tap", "spectral_norm_sepkron",
                "L1Selected", "L2Selected")

for (comp in components) {
  cat(sprintf("Comparing component '%s':\n", comp))
  comp_og <- result_og[[comp]]
  comp_cpp <- result_cpp[[comp]]
  cmp <- all.equal(comp_og, comp_cpp)
  print(cmp)
  cat("\n")
}


for (comp in components) {
  cat(sprintf("Comparing component '%s':\n", comp))
  comp_cpp_1 <- result_cpp[[comp]]
  comp_cpp_2 <- result_cpp_2[[comp]]
  cmp <- all.equal(comp_cpp_1, comp_cpp_2)
  print(cmp)
  cat("\n")
}

benchmark_results = microbenchmark(
  simulation_og = simulation_og(grid_size, params, c = 2, type = "trapezoid"),
  simulation_cpp = simulation_cpp(grid_size, params, c = 2, type = "trapezoid"),
  times = 10
)

print(benchmark_results)

#------------------------------------------------------------------------------
sigma = 1
alpha = 1
beta = 0
grid_size = 50
n1 = grid_size
n2 = grid_size
nvec = c(grid_size, grid_size)
g = length(nvec)
N = prod(nvec)
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
  seed = 3)


X = spatial_process$X
M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)

# Compute all autocovariances at once using vectorized form
# Use sapply to map through each row of the m_ij matrix
C_ij_vector = compute_autocovariance_vector_cpp(X, M_ij_matrix)

# Now reshape it back to a matrix form (filled by columns)
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute Autocovariance lag 0,0
C_00 = SpatialAutoCov_cpp(X, c(0,0))

# Compute regular GammaEst matrix normalized with C_00
GammaEst = C_ij_matrix

RhoEst = C_ij_matrix / C_00

Lv1 = bandwidth_selection_spatial_cpp_v1(RhoEst, n1=n1, n2=n2)
Lv2 = bandwidth_selection_spatial_cpp_v2(X, n1=n1, n2=n2)

all.equal(Lv1, Lv2)

benchmark_results <- microbenchmark(
  v1 = bandwidth_selection_spatial_cpp_v1(RhoEst, n1=n1, n2=n2),
  v2 = bandwidth_selection_spatial_cpp_v2(X, n1=n1, n2=n2),
  times =50
)
print(benchmark_results)
