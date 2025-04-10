# PARALLEL with MULTIPLE SAMPLES and C++ FUNCTIONS ###
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
if (!require(Rcpp)) install.packages("Rcpp"); library(Rcpp)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(wd,".."))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions
sourceCpp("src/estimators.cpp")
#-----------------------------------

one_simulation_gridsearch = function(params,
                                     l,
                                     set_seed,
                                     grid_size,
                                     type,
                                     c){
  
  # Create vector with the grid size (square grid)
  nvec = c(grid_size, grid_size)
  g = length(nvec) # number of spatial dimension (for now g = 2)
  N = prod(nvec)   # Total size of the grid n1 x n2
  
  # 1) Simulate the spatial process given the parameters and 
  #    grid sizes
  spatial_process = simulate_spatial_process(
    covariance_function = ModifiedExponentialCovariance,
    grid_size = grid_size,
    params = params,
    seed = set_seed)
  
  X = spatial_process$X                 # Spatial Process
  true_cov = spatial_process$covariance # True Covariance
  
  # 2) Computing M Matrix for lags corrected with the m functions
  M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)
  
  # 3) Compute all autocovariances
  C_ij_vector = compute_autocovariance_vector_cpp(X, M_ij_matrix)
  
  # Reshape back to a Matrix
  C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)
  
  # Compute Autocovariance lag 0,0
  C_00 = SpatialAutoCov_cpp(X, c(0,0))
  
  # Naive Autocovariance Estimator
  GammaEst = C_ij_matrix
  # Naive Autocorrelation Estimator
  RhoEst = C_ij_matrix / C_00
  
  # Assign same l for both spatial dimensions
  L = c(l, l)
  # 4) Compute Taper matrix
  kappa_ij_vector = compute_multi_taper_vector_cpp(M_ij_matrix,
                                                   L = L,
                                                   c = c,
                                                   type = type)
  
  # Turn into matrix
  kappa_ij_matrix = matrix(kappa_ij_vector, nrow = N, ncol = N)
  
  # Compute Tapered Covariance Matrix
  GammaEstTaper = kappa_ij_matrix * GammaEst
  
  # 6) Compute Separable Taper Estimator
  SepResults = Tapered_Sep_Autocovariance_Kron_multi_cpp(
    X,
    c   = c,
    L   = L,
    type= type)
  
  
  GammaEstTaperSep = SepResults$KronTaperCov
  
  truecov_matrices = true_cov
  cov_matrices = GammaEst
  taper_matrices = kappa_ij_matrix
  taper_covariances = GammaEstTaper
  RhoEst = RhoEst
  spectral_norm = norm(true_cov - GammaEst, type = "2")
  spectral_norm_tap = norm(true_cov - GammaEstTaper, type = "2")
  spectral_norm_sepkron = norm(true_cov - GammaEstTaperSep, type = "2")
  
  # Build data frame for the current lambda and all grid sizes
  results_df = data.frame(
    GridSize           = grid_size,
    seed               = set_seed,
    l                  = l,
    alpha              = params$alpha1,
    lambda             = params$lambda1,
    beta               = params$beta,
    type               = type,
    c                  = c,
    SNorm              = spectral_norm,
    SNormTapered       = spectral_norm_tap,
    SNormSeparTapered  = spectral_norm_sepkron
  )
  
  return(results_df)
}
# RUNNING THE SIMULATION FOR MULTIPLE SAMPLES

nSim = 100  # Number of simulations
Ns = c(10,20,30,40,50) # Number of grid sizes
sigma = 1
alpha = 1
lambdas = seq(from = 2, to = 6, by = 2)
Ls = seq(from = 1, to = 18, by = 1)
beta = 0

first_params = list(sigma = sigma,
                    alpha1 = alpha,
                    alpha2 = alpha,
                    beta = beta,
                    test_sep = F)

# Taper parameters
type = "trapezoid"
c = 2
final_df = data.frame()

# Detect the number of cores
num_cores = detectCores()

# Create and register a cluster
cl = makeCluster(num_cores - 2)
clusterEvalQ(cl, {
  library(Rcpp)
  library(MASS)
  library(tidyr)
  library(pracma)
  library(sp)
  library(ncf)
  library(fields)
  library(reshape2)
  library(doSNOW)
  library(parallelDist)
  source("R/utils.R")           
  source("R/covariance_funs.R") 
  source("R/estimators.R")      
  sourceCpp("src/estimators.cpp")
})

clusterExport(cl, "one_simulation_gridsearch")
#clusterCall(cl, sourceCpp("src/estimators.cpp"))
registerDoParallel(cl)
fun_list = my_cpp_source_funs("src/estimators.cpp")
cfun = fun_list$functions


tic("Start Simulation: ")
for(l in Ls){
tic("Start bandwidth: ")
cat("Bandwidth: ", l, "\n")
  for(grid_size in Ns){
    tic("Start Grid Size: ")
    # Loop for grid size
    cat("Grid Size: ", grid_size, "\n")
    
    for(lambda in lambdas){
      tic("Start Lambda: ")
      # Loop for Lambdas
      cat("Parameter Lambda: ", lambda, "\n")
      params = append(first_params, list(lambda1 = lambda, lambda2 = lambda))
      
      RESULTS = foreach(i = 1:nSim, 
                        .noexport = cfun,
                        .combine = rbind) %dopar% {
                          
                          out_df = one_simulation_gridsearch(
                            params = params,
                            l = l,
                            set_seed = i,
                            grid_size = grid_size,
                            type = type,
                            c = c)
                          
                          return(out_df)
                        }                   
      
      final_df = rbind(final_df, RESULTS)
      toc()
    }
    toc()
  }
toc()
}
toc()
stopCluster(cl)


#-------------------------------------------------------------------------------
# Save the results
resdir = file.path(wd,"results")

timestamp = format(Sys.time(), "%Y-%m-%d_%H%M")
file_name = paste0("results_gridsearch_", timestamp, ".csv")


write.csv(final_df, file = file.path(resdir, file_name),
          row.names = F)
