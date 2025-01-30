#-------------------------------------------------------------------------------
### SPATIAL PROCESS SIMULATION ###
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
if (!require(parallel)) install.packages("parallel"); library(parallel)
if (!require(reshape2)) install.packages("reshape2"); library(reshape2)
if (!require(ncf)) install.packages("ncf"); library(ncf)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
setwd(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),".."))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions
set.seed(42)
#-------------------------------------------------------------------------------
# GENERATING GRID

# Define the dimensions of the 2D grid
n1 = 5 # Dimension 1 size 
n2 = 5 # Dimension 2 size
nvec = c(n1, n2)
N = prod(nvec) 
# Generate spatial locations/coordinates of integers
grid = expand.grid(t1 = 1:n1, t2 = 1:n2)
#-------------------------------------------------------------------------------
# GENERATING ANALYTICAL TRUE COVARINACE
#-------------------------------------------------------------------------------
# Exponential

# decay = 0.5 # Use higher than 1 decay otherwise small spatial correlation!
# sigma = 1
# true_cov = ExponentialCovariance(grid, sigma, decay, method = "abs_diff")
# is_positive_semi_definite(true_cov)
# plot_matrix(true_cov, main = "Exponential True Covariance", labels = F)
#-------------------------------------------------------------------------------
# Gneiting

# decay1 = 0.2 # Use lower decay values for strong spatial correlation
# decay2 = 0.2
# sigma = 1
# beta = 1
# 
# true_cov = GneitingSCovariance(grid,sigma,
#                                decay1, decay2,
#                                beta)
# is_positive_semi_definite(true_cov)
# #plot_matrix(true_cov, main = "True Covariance", labels =F)
#-------------------------------------------------------------------------------
# MODIFIED EXPONENTIAL

# Smoothness parameter:
#   high makes the spatial process smooth
#   low makes the spatial process rougher
alpha1 = 1
alpha2 = 1

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

is_positive_semi_definite(true_cov)
title_cov = paste("True Covariance Modifed Exp; beta = ", beta,
                  "alpha =", alpha1, "lambda =", lambda1)
plot_matrix(true_cov, main = title_cov,  labels = F)
#-------------------------------------------------------------------------------
# SIMULATION
# Simulate Spatial Process with given true covariance
set.seed(42)
mvn_sim = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_cov) 
# Checking that the simulated spatial data "falls" in the correct coordinates
#grid$sim = mvn_sim 

# NOTE: the process is structured from top left instead of from bottom left
X = matrix(mvn_sim, nrow = n1, ncol = n2, byrow = T)
title = paste("Spatial Process Simulated with Exp; beta = ", beta,
             "alpha =", alpha1, "lambda =", lambda1)
#plot_matrix(X, labels = T, main = title)
#-------------------------------------------------------------------------------
# Taper parameters
c = 1
l = 1
type = "rectangular" # rectangular or trapezoid

GammaEst_loop = matrix(NA, nrow = N, ncol = N)
GammaEstTap_loop = matrix(NA, nrow = N, ncol = N)
C_00 = SpatialAutoCov(X, c(0,0))
M_ij_list = vector("list", N*N)
k_l_mat = matrix(NA, nrow = N, ncol = N)
for(i in 1:N){
  for(j in 1:N){
    m_ij = compute_m_values(i, j, nvec,
                           flipsign = TRUE,
                           flipposition = TRUE)
    M_ij_list[[(i-1)*N + j]] = m_ij
    C_ij = SpatialAutoCov(X, m_ij)
    GammaEst_loop[i, j] =  C_ij / C_00
    k_l_mat[i, j] = flat_top_taper(m_ij, c=c, l=l, type = type)
    GammaEstTap_loop[i, j] = k_l_mat[i, j] * (C_ij / C_00)
  }
}

norm(true_cov - GammaEst_loop, type = "2")
norm(true_cov - GammaEstTap_loop, type = "2")
norm(true_cov - GammaEst_loop, type = "F")
norm(true_cov - GammaEstTap_loop, type = "F")
# plot_matrix(GammaEst_loop, main = "Est Covariance")
# plot_matrix(k_l_mat, main = "Tapers")
# plot_matrix(GammaEstTap_loop, main = "Tapered Covariance")
#-------------------------------------------------------------------------------


tic()
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

# Now reshape it back to a matrix form
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

C_00 = SpatialAutoCov(X, c(0,0))
# Compute regular GammaEst matrix
GammaEst = C_ij_matrix / C_00

k_ij_vector = sapply(1:nrow(M_ij_matrix), 
                     function(idx) flat_top_taper(M_ij_matrix[idx, ],
                                                  c=c, l=l, type = type))
k_ij_matrix = matrix(k_ij_vector, nrow=N, ncol=N)

GammaEstTaper = k_ij_matrix * GammaEst
toc()


norm(GammaEst_loop - GammaEst, type = "2")
norm(GammaEstTap_loop - GammaEstTaper, type = "2")
norm(k_l_mat - k_ij_matrix, type = "2")
#-------------------------------------------------------------------------------
# ESTIMATION OF THE COVARIANCE MATRIX
# Taper parameters
c = 2
l = 1
type = "trapezoid" # rectangular or trapezoid

tic("Estimation of the covariance matrix")
M_ij_list = vector("list", N*N)
for(i in 1:N){
  for(j in 1:N){
    M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec,
                                                flipsign = TRUE,
                                                flipposition = TRUE)
  }
}
# Autocovariance at lag (0, 0)  
C_00 = SpatialAutoCov(X, c(0, 0))  

# Vectorize m_ij computation by converting list to matrix
M_ij_matrix = do.call(rbind, M_ij_list)

# Compute all autocovariances at once using vectorized form

C_ij_vector = sapply(1:nrow(M_ij_matrix), function(idx) 
                  SpatialAutoCov(X, M_ij_matrix[idx, ]))
                                                                         

# Now reshape back to a matrix form
C_ij_matrix = matrix(C_ij_vector, nrow=N, ncol=N)

# Compute regular GammaEst matrix
GammaEst = C_ij_matrix / C_00
toc()

norm(true_cov - GammaEst, type = "F")
norm(true_cov - GammaEst, type = "2")


#-------------------------------------------------------------------------------
C_00 = SpatialAutoCov(X, c(0, 0)) 

# Verione 1
GammaEst = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
for(i in 1:nrow(grid)){
  for(j in 1:nrow(grid)){
    m_ij = compute_m_values(i, j, nvec,
                            flipsign = TRUE,
                            flipposition = TRUE)
    GammaEst[i, j] = SpatialAutoCov(X, m_ij) 
  }
}

# Verione 2
GammaEst2 = matrix(NA, nrow = nrow(grid), ncol = nrow(grid))
for(i in 1:nrow(grid)){
  for(j in 1:nrow(grid)){
    m_ij = compute_m_values(i, j, nvec,
                            flipsign = F,
                            flipposition = F)
    GammaEst2[i, j] = SpatialAutoCov(X, m_ij) 
  }
}
View(GammaEst / C_00)
View(GammaEst2 / C_00)


norm(true_cov - (GammaEst / C_00), type = "2")
norm(true_cov - (GammaEst2 / C_00), type  = "2")

hvec = compute_m_values(1, 2, nvec,
                 flipsign = T,
                 flipposition = T)


h1 = hvec[1] # lag1
h2 = hvec[2] # lag2

X_mean = mean(as.vector(X))
n1 = dim(X)[1]
n2 = dim(X)[2]
N = n1*n2

#Finding the regions on which to shift the lags for the autocorrelation
# Horizontal lag
T11 = max(1, 1-h1)
T12 = min(n1, n1-h1)

# Vertical lag
T21 = max(1, 1-h2)
T22 = min(n2, n2-h2)

sub_matrix1 = X[T21:T22, T11:T12]
sub_matrix2 = X[(T21+h2):(T22+h2), (T11+h1):(T12+h1)]



for(K_h1 in 0:n1){
  for(J_h2 in 0:n2){
    
  }
}



