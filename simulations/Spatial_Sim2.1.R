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
### EXAMPLE TRUE COVARIANCE MATRIX STRUCTURE

# MODIFIED EXPONENTIAL

# Define the dimensions of the 2D grid
n1 = 30 # Dimension 1 size
n2 = 30 # Dimension 2 size
nvec = c(n1, n2)
N = prod(nvec)
# Generate spatial locations/coordinates of integers
grid = expand.grid(t1 = 1:n1, t2 = 1:n2)

# Smoothness parameter:
#   high makes the spatial process smooth
#   low makes the spatial process rougher
alpha1 = 1
alpha2 = 1

# Range parameter:
#   high makes the covariance decay more slowly with distance
#   low makes the covariance decay more quickly and with distance

lambda1 = 6
lambda2 = 6

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
#is_positive_semi_definite(true_cov)
title_cov = paste("True Covariance Modifed Exp; beta = ", beta,
                  "alpha =", alpha1, "lambda =", lambda1)
file_path <- file.path(wd, paste0("plots/Cova_",
                                  "beta", beta, "_",
                                  "alpha", alpha1, "_",
                                  "lambda", lambda1, ".png"))
#png(file_path, width = 800, height = 500)
plot_matrix(true_cov, main = title_cov,  labels = F)
#dev.off()
#-------------------------------------------------------------------------------
# SIMULATION EXAMPLE
# params = list(sigma = 1,
#               alpha1 = 1,
#               alpha2 = 1,
#               lambda1 = 2,
#               lambda2 = 2,
#               beta = 0,
#               test_sep= F)
# spatial_process = simulate_spatial_process(
#                       covariance_function = ModifiedExponentialCovariance,
#                       grid_size = 50,
#                       params = params,
#                       seed = 42)
# X = spatial_process$X
# true_cov = spatial_process$covariance
# title_cov = paste("true covariance modifed exp; beta = ", beta,
#                    "alpha =", alpha1, "lambda =", lambda1)
# plot_matrix(true_cov, main = title_cov,  labels = F)
# 
# title = paste("spatial process simulated with exp; beta = ", beta,
#                "alpha =", alpha1, "lambda =", lambda1)
# plot_matrix(X, labels = T, main = title)
#-------------------------------------------------------------------------------
# DEMONSTRATION M FUNCTION
# 
# Xvec = as.vector(X) # column wise vectorization
# i = 1 # Position i in the vectorized X
# j = 6 # Position j in the vectorized X
# c(m_fun(i, s=1, d=1, nvec) - m_fun(j, s=1, d=1, nvec),
#   m_fun(i, s=2, d=1, nvec) - m_fun(j, s=2, d=1, nvec))
# # OPINION: the correct way is to keep positive lags for movements that go
# # to the right and below. This presumes origin at point 1,1 being top left
# # corner.
# # In addition, consider that we need to flip the position of the generated lags
# # such that:
# #   first number = horizontal shift (x-axis)
# #   second number = vertical shift (y-axis)
# compute_m_values(i,j,nvec, d=1, flipsign = T, flipposition = T)
#-------------------------------------------------------------------------------
# DEMONSTRATION FLAT TOP TAPERS

# Xvec = as.vector(X) # column wise vectorization
# i = 1 # Position i in the vectorized X
# j = 900 # Position j in the vectorized X
# 
# # Taper parameters
# c = 1
# l = 1
# type = "rectangular"
# 
# flat_top_taper(m_ij, c=c, l=l,type=type) 

#-------------------------------------------------------------------------------
# LOOP SIMULATION FOR INCREASING GRID / SAMPLE SIZE

#alphas = seq(from = 0.5, to = 1.5, by = 0.25)
sigma = 1
alpha = 1
lambdas = seq(from = 2, to = 10, by = 2)
#lambda = 5
beta = 0
# Taper parameters
type = "rectangular"
c = 1

#l = 5 # should increase with increasing grid size
bandwidths = seq(from = 13, to = 21, by = 2)
final_df = data.frame()


tic("Whole simulations loop")
for(lambda in lambdas){
  tic("Time for lambdas")
  cat("\nPARAMETER lambda: ", lambda, "\n") 
  for(l in bandwidths){
    tic("Time for bandwidth")
    cat("\nPARAMETER l: ", l, "\n") # change name of parameter accordingly
    params = list(sigma = sigma,
                  alpha1 = alpha,
                  alpha2 = alpha,
                  lambda1 = lambda,
                  lambda2 = lambda,
                  beta = beta,
                  test_sep = F)
  
    
    tic("Loop for different grid sizes")
    Ns = c(10, 20, 30, 40, 50, 60)  # Grid sizes
    # Initialize vectors and lists to store results
    truecov_matrices = list()
    cov_matrices = list() 
    taper_matrices = list()
    taper_covariances = list()
    spectral_norm = numeric(length(Ns))
    spectral_norm_tap = numeric(length(Ns))
    spectral_norm_sepkron = numeric(length(Ns))
    # frob_norm = numeric(length(Ns)) 
    # frob_norm_tap = numeric(length(Ns)) 
    
    # Looping a simulation for each grid size
    for (iter in 1:length(Ns)){
      tic("Time for interation")
      cat("\nLoop iteration Grid size:", Ns[iter], "\n")
      grid_size = Ns[iter]
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
      GammaEst = C_ij_matrix / C_00
      
      # Compute Taper vector
      kappa_ij_vector = sapply(1:nrow(M_ij_matrix), 
                           function(idx) flat_top_taper(M_ij_matrix[idx, ],
                                                        c=c, l=l, type = type))
      # Compute taper matrix
      kappa_ij_matrix = matrix(kappa_ij_vector, nrow=N, ncol=N)
      
      # Compute Tapered Covariance Matrix
      GammaEstTaper = kappa_ij_matrix * GammaEst
      
      # Compute Separable Taper Estimator
      SepResults= Tapered_Sep_Autocovariance_Kron(X, l=l, c=c, type=type)
      GammaEstTaperSep = SepResults$KronTaperCov
      
     # Compute the Spectral norm between true_cov and estimated GammaEst_2
      truecov_matrices[iter] = true_cov
      cov_matrices[iter] = GammaEst
      taper_matrices[iter] = kappa_ij_matrix
      taper_covariances[iter] = GammaEstTaper
      spectral_norm[iter] = norm(true_cov - GammaEst, type = "2")
      spectral_norm_tap[iter] = norm(true_cov - GammaEstTaper, type = "2")
      spectral_norm_sepkron[iter] = norm(true_cov - GammaEstTaperSep, type = "2")
      toc()
    }
  toc()
  
  
  # Preparing data frame for Excel export
  results_df = data.frame(
    GridSize = Ns,
    alpha = alpha,
    lambda = lambda,
    beta = beta,
    type = type,
    c = c,
    l = l,
    SpectralNorm = spectral_norm,
    SpectralNormTaper = spectral_norm_tap,
    SpectralNormSepTaper = spectral_norm_sepkron
  )
  
  final_df = rbind(final_df, results_df)
  toc()
  }
toc()
}
toc()


# Print results
write.csv(final_df, paste0(wd, "/results/simulations_20250118.csv"))
#write.xlsx(final_df, file = paste0(wd, "/results/simulations_week40_v2.xlsx"))

for (iter in 1:length(Ns)) {
  cat("\nGrid size:", Ns[iter])
  cat("\nSpectral norm distance between true and estimated covariance:",
      spectral_norm[iter], "\n")
  cat("\nSpectral norm distance between true and taperd estimated covariance:",
      spectral_norm_tap[iter], "\n")
  cat("\nSpectral norm distance between true and separable tapered estimated
      covariance:", spectral_norm_tap[iter], "\n")

}


file_name <- paste0("sim_res_c_", c, "_l_", l,
                    "_alpha_", alpha, "_lambda_", lambda,
                    ".xlsx")
write_xlsx(results_df, path = paste0("./results/", file_name))