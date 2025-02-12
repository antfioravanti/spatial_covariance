#-------------------------------------------------------------------------------
# COVARIANCE FUNCTIONS & SIMULATIONS
#-------------------------------------------------------------------------------


if (!require(pracma)) install.packages("pracma"); library(pracma)
if (!require(MASS)) install.packages("MASS"); library(MASS)
source("R/estimators.R") # Estimators
source("R/utils.R") # Estimators
#-------------------------------------------------------------------------------
# GRIDS
#-------------------------------------------------------------------------------
generate_square_grid = function(size){
  # Define the dimensions of the 2D grid
  n1 = size # Dimension 1 size 
  n2 = size # Dimension 2 size
  nvec = c(n1, n2)
  N = prod(nvec) 
  # Generate spatial locations/coordinates of integers
  grid = expand.grid(t1 = 1:n1, t2 = 1:n2)
  return(grid)
}

#-------------------------------------------------------------------------------
# SPATIAL
#-------------------------------------------------------------------------------
# EXPONENTIAL
ExponentialCovariance = function(grid, sigma, decay, method = "abs_diff"){
  # Check that grid has the correct column names "t1" and "t2"
  expected_names <- c("t1", "t2")
  
  if (!all(colnames(grid) == expected_names)) {
    stop("Error: The grid must have column names 't1' and 't2'.")
  }
  
  # Check that the method is among the allowed ones
  if(!sum(method == c("diff", "euc", "abs_diff"))){
    stop("Specify either euc, abs_diff or diff method")
  }
  # Absolute Difference as Lag
  
  if(method == "abs_diff"){
    lag1 = abs(outer(grid$t1, grid$t1, "-"))
    lag2 = abs(outer(grid$t2, grid$t2, "-"))
    covariance = sigma^2 * exp(-lag1/decay) * exp(-lag2/decay)
    
  }else if (method == "euc"){
    # Euclidean Distance as lag
    lag = as.matrix(dist(grid, method = "euclidean"))
    covariance = sigma^2 * exp(-lag/decay)
  }
  return(covariance)
}

GneitingSCovariance = function(grid, sigma, decay1, decay2, beta){
  # Check that grid has the correct column names "t1" and "t2"
  expected_names <- c("t1", "t2")
  
  if (!all(colnames(grid) == expected_names)) {
    stop("Error: The grid must have column names 't1' and 't2'.")
  }
  
  h1 = outer(grid$t1, grid$t1, "-")
  h2 = outer(grid$t2, grid$t2, "-")
  
  covariance = (sigma^2 / (decay2 * abs(h2) + 1)) * 
    exp(-decay1 * abs(h1) / (decay2 * abs(h2) + 1)^(beta / 2))
  return(covariance)
}

# MODIFIED EXPONENTIAL
ModifiedExponentialCovariance = function(grid,
                                         sigma = 1,
                                         alpha1 = 1,
                                         alpha2 = 1,
                                         lambda1 = 1,
                                         lambda2 = 1,
                                         beta = 0,
                                         test_sep = F) {
  # Args
  #   alpha1 and alpha 2 control the smoothness in each dimension
  #   lambda1 and lamda 2 control the range in each dimension
  #   beta controls the separability
  #     beta = 0  separable
  #     beta = 1 non-separable
  
  
  # Check that grid has the correct column names "t1" and "t2"
  expected_names <- c("t1", "t2")
  
  if (!all(colnames(grid) == expected_names)) {
    stop("Error: The grid must have column names 't1' and 't2'.")
  }
  
  # Compute pairwise differences for each dimension
  h1_diff = outer(grid$t1, grid$t1, "-")  # Difference in first dimension
  h2_diff = outer(grid$t2, grid$t2, "-")  # Difference in second dimension
  
  cov = sigma^2 * exp(-(abs(h1_diff)^alpha1 / lambda1 +
                          abs(h2_diff)^alpha2 / lambda2 + 
                          beta * abs(h1_diff - h2_diff)))
  if(isTRUE(test_sep)){
    cov_sep = sigma^2*exp(-(abs(h1_diff)^alpha1 / lambda1)-beta*abs(h1_diff))* 
      exp(-(abs(h2_diff)^alpha2 / lambda2) -beta*abs(h2_diff))
    
    norm_diff = norm(cov_sep - cov, type = "2")
    return(list(covariance = cov, norm_diff = norm_diff))
  }else{
    return(cov)
  }
}

#-------------------------------------------------------------------------------
# GENERATE SIMULATION
#-------------------------------------------------------------------------------
simulate_spatial_process = function(covariance_function, 
                                    grid_size, 
                                    params = list(),
                                    seed = 42) {
  
  # Generate a square grid
  grid = generate_square_grid(grid_size)
  
  # Unpack parameters and ensure they are passed correctly
  true_cov = do.call(covariance_function, c(list(grid = grid), params))
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Simulate the spatial process with the generated covariance matrix
  mvn_sim = mvrnorm(n = 1, mu = rep(0, nrow(grid)), Sigma = true_cov)
  
  # Reshape the simulated data into a matrix corresponding to the grid
  X = matrix(mvn_sim, nrow = grid_size, ncol = grid_size, byrow = TRUE)
  
  # Return the grid and the simulated process matrix
  return(list(grid = grid, X = X, covariance = true_cov))
}