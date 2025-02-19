#-------------------------------------------------------------------------------
### SPATIAL PROCESS SIMULATION - PARALLEL MULTIPLE SAMPLES###
#-------------------------------------------------------------------------------
options(scipen = 5) # some decimals, set to 999 for max decimals
#options(scipen = 0) # scientific display
if (!require(sp)) install.packages("sp"); library(sp)
if (!require(tidyr)) install.packages("tidyr"); library(tidyr)
if (!require(gstat)) install.packages("gstat"); library(gstat)
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(fields)) install.packages("fields"); library(fields)
if (!require(rstudioapi)) install.packages("rstudioapi"); library(rstudioapi)
if (!require(pracma)) install.packages("pracma"); library(pracma)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(ncf)) install.packages("ncf"); library(ncf)
if (!require(openxlsx)) install.packages("openxlsx"); library(openxlsx)
if (!require(parallel)) install.packages("parallel"); library(parallel)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(parallelDist)) install.packages("parallelDist"); library(parallelDist)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),".."))
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions
#-------------------------------------------------------------------------------
# ONE SIMULATION
#-------------------------------------------------------------------------------

run_one_simulation = function(grid_size = 10,
                              set_seed = 42,
                              sigma = 1,
                              beta = 0,
                              alpha = 1,
                              lambda = 2,
                              type = "rectangular",
                              c = 1) {
  
  # Param list for ModifiedExponentialCovariance
  params = list(
    sigma = sigma,
    alpha1 = alpha,
    alpha2 = alpha,
    lambda1 = lambda,
    lambda2 = lambda,
    beta = beta,
    test_sep = FALSE
  )
  
  # 1) Simulate the spatial process
  spatial_process = simulate_spatial_process(
    covariance_function = ModifiedExponentialCovariance,
    grid_size = grid_size,
    params = params,
    seed = set_seed)
  
  X = spatial_process$X
  true_cov = spatial_process$covariance
  
  # 2) Build the matrix of lags M_ij
  nvec = c(grid_size, grid_size)
  N = prod(nvec)
  
  M_ij_list = vector("list", N*N)
  for(i in 1:N){
    for(j in 1:N){
      M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec,
                                                  flipsign = TRUE,
                                                  flipposition = TRUE)
    }
  }
  M_ij_matrix = do.call(rbind, M_ij_list)
  
  # 3) Compute all sample autocovariances
  C_ij_vector = sapply(seq_len(nrow(M_ij_matrix)), 
                       function(idx) SpatialAutoCov(X, M_ij_matrix[idx, ]))
  
  # Construct the autocovariance matrix
  C_ij_matrix = matrix(C_ij_vector, nrow = N, ncol = N)
  
  # 4) Normalize by variance C_00
  C_00    = SpatialAutoCov(X, c(0,0))
  GammaEst = C_ij_matrix / C_00 # ?????
  
  # 5) Select empirical bandwidth for each dimension
  Lvals = bandwidth_selection_spatial(
    corrMat = GammaEst,
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
  
  # 7) Build Tapered Covariance
  GammaEstTaper = kappa_ij_matrix * GammaEst
  
  # 8) Compute Separable Taper Estimator
  SepResults = Tapered_Sep_Autocovariance_Kron_multi(
    X,
    L   = L,
    c   = c,
    type= type)
  
  GammaEstTaperSep = SepResults$KronTaperCov
  
  # 9) Compute matrix norms:
  spectral_norm = norm(true_cov - GammaEst, type = "2")
  frob_norm = norm(true_cov - GammaEst, type = "F")
  
  spectral_norm_tap = norm(true_cov - GammaEstTaper, type = "2")
  frob_norm_tap = norm(true_cov - GammaEstTaper, type = "F")
  
  spectral_norm_sepkron = norm(true_cov - GammaEstTaperSep, type = "2")
  frob_norm_sepkron = norm(true_cov - GammaEstTaperSep, type = "F")
  
  selected_l1  = L[1]
  selected_l2  = L[2]
  
  # Build data frame for the current lambda and all grid sizes
  results_df = data.frame(
    GridSize           = grid_size,
    seed               = set_seed,
    alpha              = alpha,
    lambda             = lambda,
    beta               = beta,
    type               = type,
    SpectralNorm       = spectral_norm,
    SpectralNormTaper  = spectral_norm_tap,
    FrobNorm           = frob_norm,
    FrobNormTaper      = frob_norm_tap,
    SpectralNormSepKron = spectral_norm_sepkron,
    FrobNormSepKron    = frob_norm_sepkron,
    L1Selected         = selected_l1,
    L2Selected         = selected_l2
  )
  
  return(results_df)
}


#-------------------------------------------------------------------------------
# RUNNING THE SIMULATION FOR MULTIPLE SAMPLES

nSim = 100  # Number of simulations
Ns = c(10,20,30,40,50,60) # Number of grid sizes
sigma = 1
alpha = 1
lambdas = seq(from = 2, to = 6, by = 2)
beta = 0
# Taper parameters
type = "rectangular"
c = 1
final_df = data.frame()

# Detect the number of cores
num_cores = detectCores()

# Create and register a cluster
cl = makeCluster(num_cores - 2)
clusterEvalQ(cl, {
  source("R/utils.R")           
  source("R/covariance_funs.R") 
  source("R/estimators.R")      
  source("R/plotting.R") 
  library(MASS)
  library(tidyr)
  library(pracma)
  library(sp)
  library(ncf)
  library(fields)
  library(reshape2)
  library(doSNOW)
  library(parallelDist)
  })
registerDoParallel(cl)
clusterExport(cl, "run_one_simulation")

tic("Start Simulation: ")
for(grid_size in Ns){
  tic("Start Grid Size: ")
  # Loop for grid size
  cat("Grid Size: ", grid_size, "\n")
  
  for(lambda in lambdas){
    tic("Start Lambda: ")
    # Loop for Lambdas
    cat("Parameter Lambda: ", lambda, "\n")
  
  
    RESULTS = foreach(i = 1:nSim, 
                         .combine = rbind) %dopar% {

      out_df = run_one_simulation(
        grid_size = grid_size,
        set_seed = i,
        sigma = 1,
        beta = 0,
        lambda = lambda,
        alpha = 1,
        type = "rectangular",
        c = 1)
      
      return(out_df)
    }                   

    final_df = rbind(final_df, RESULTS)
    toc()
  }
  toc()
}
stopCluster(cl)
toc()


