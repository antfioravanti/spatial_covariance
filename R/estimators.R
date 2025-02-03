### ESTIMATORS KERNELS AND NONPARAMETRIC FUNCTIONS ###
if (!require(pracma)) install.packages("pracma"); library(pracma)
if (!require(MASS)) install.packages("MASS"); library(MASS)
if (!require(foreach)) install.packages("foreach"); library(foreach)
if (!require(parallel)) install.packages("parallel"); library(parallel)
if (!require(doParallel)) install.packages("doParallel"); library(doParallel)

#-------------------------------------------------------------------------------
# M FUNCTIONS
#-------------------------------------------------------------------------------

# M0 function
m0_fun = function(k, d=1) {
  # M0 search function
  return((k - 1) %% d + 1)
}

# M function
m_fun = function(k, s, d=1, nvec) {
  # M search function. 
  # k is the value
  # s is the counter for the spatial dimensions 1,..., g
  # d is the dimension of the multivariate process (assume univariate d=1)
  # nvec is the vector containing the sizes of the spatial dimensions
  
  if (s == 1) {
    m = ((ceiling(k / d) - 1) %% (d * nvec[1])) + 1
  } else {
    prod_prev_n = prod(nvec[1:(s-1)])
    m = ((ceiling(k / (d * prod_prev_n)) - 1) %% 
           (d * prod_prev_n * nvec[s])) + 1
  }
  return(m)
}

compute_m0_values = function(i, j, d=1) {
  # Function to compute the m0 values
  # We are vectorizing ROW WISE
  m0_i <- m0_fun(i, d)
  m0_j <- m0_fun(j, d)
  m0_ij <- c(m0_i, m0_j)
  return(m0_ij)
}

compute_m_values = function(i,j, nvec, d=1,
                            flipsign = TRUE,
                            flipposition = FALSE) {
  # Function to get the vector M with the right lags
  g = length(nvec)
  m_ij = numeric(g) # dimension of the spatial process
  for(s in 1:g){ # s is the counter of the spatial process
    m_ij[s] =  m_fun(i, s, d, nvec) - m_fun(j, s, d, nvec)
  }
  if(flipsign==TRUE){
    # Changing the direction of the lag
    m_ij = m_ij * -1
  }
  if(flipposition==TRUE){
    # Reversing the lag ordering
    m_ij = rev(m_ij)
  }
  return(m_ij)
}

compute_M_matrix = function(N, nvec, flipsign = TRUE, flipposition = FALSE){
    # Computing M Matrix for correct lags
    M_ij_list = vector("list", N*N)
    for(i in 1:N){
      for(j in 1:N){
        M_ij_list[[(i-1)*N + j]] = compute_m_values(i, j, nvec,
                                                    flipsign = flipsign,
                                                    flipposition = flipposition)
        }
      }
    
    # Vectorize m_ij computation by converting list to matrix
    M_ij_matrix = do.call(rbind, M_ij_list)
    return(M_ij_matrix)
}



# Vectorized computation of m_ij values with correct order
compute_m_values_vectorized_corrected = function(grid, nvec, d=1,
                                                 flipsign=TRUE,
                                                 flipposition=FALSE) {
  
  # Number of spatial dimensions
  g <- length(nvec)
  
  # Function to compute m_fun for multiple values of k
  m_fun_vectorized <- function(k, s, d, nvec) {
    if (s == 1) {
      m <- ((ceiling(k / d) - 1) %% (d * nvec[1])) + 1
    } else {
      prod_prev_n <- prod(nvec[1:(s-1)])
      m <- ((ceiling(k / (d * prod_prev_n)) - 1) %% (d * prod_prev_n * nvec[s])) + 1
    }
    return(m)
  }
  
  # Initialize a matrix to store m_ij values
  m_ij_matrix <- matrix(NA, nrow = nrow(grid), ncol = g)
  
  # Loop over the spatial dimensions to compute the m values for all (i, j) pairs
  for (s in 1:g) {
    # Calculate m_i and m_j for dimension s
    m_i_s <- m_fun_vectorized(grid$t1, s, d, nvec)
    m_j_s <- m_fun_vectorized(grid$t2, s, d, nvec)
    
    # Store the difference m_ij for each pair (i, j) and dimension s
    m_ij_matrix[, s] <- m_i_s - m_j_s
  }
  
  # Apply the flipsign and flipposition if needed
  if (flipsign) {
    m_ij_matrix <- m_ij_matrix * -1
  }
  if (flipposition) {
    m_ij_matrix <- t(apply(m_ij_matrix, 1, rev))
  }
  
  # Reshape the matrix to match the structure of the original for loop
  m_ij_matrix <- matrix(m_ij_matrix, nrow = N * N, ncol = g, byrow = TRUE)
  
  return(m_ij_matrix)
}


#-------------------------------------------------------------------------------
# TAPERS
#-------------------------------------------------------------------------------

# flat_top_taper_v0 = function(x, c=1, l=1, type = "rectangular") {
#   # Rectangular kernel
#   rectangular_kernel = function(x, c) {
#     if (all(abs(x) <= c)) {
#       return(1)
#     } else {
#       return(0)
#     }
#   }
#   
#   # Trapezoid kernel
#   trapezoid_kernel = function(xi, c) {
#     return(max(min(c - abs(xi), 1), 0))
#   }
#   
#   # Handle l as a number or matrix
#   if (is.matrix(l)) {
#     # Assume third column is then the correct bandwidth
#     xl <- x[, 1:2] / x[, 3]
#   } else {
#     xl <- x / l
#   }
#   
#   # Select the taper type
#   if (type == "rectangular") {
#     kappa <- rectangular_kernel(xl, c)
#   } else if (type == "trapezoid") {
#     kappa <- prod(sapply(xl, trapezoid_kernel, c=c))
#   } else if (type == "product") {
#     # Product kernel implementation can be added here
#     stop("Product kernel is not implemented yet.")
#   } else {
#     stop("Invalid type. Choose from 'rectangular', 'trapezoid', or 'product'.")
#   }

#   } else if (type == "product") {
#     # # Product kernel
#     # kernel <- function(x, c_kappa = rep(2, g)) {
#     #   f <- function(xi) {
#     #     if (abs(xi) <= 1) {
#     #       return(1)
#     #     } else if (abs(xi) > c_kappa[which(x == xi)]) {
#     #       return(0)
#     #     } else {
#     #       return(min(2 - abs(xi), 1))
#     #     }
#     #   }
#     #   prod(sapply(x, f))
#     # }
#   } else {
#     stop("Invalid type. Choose from 'rectangular', 'trapezoid', or 'product'.")
#   }

#   return(kappa)
# }


flat_top_taper = function(xvec, c=1, l=1, type = "rectangular") {
  # Define the rectangular kernel
  rectangular_kernel = function(xvec, c) {
    if (all(abs(xvec) <= c)) {
      return(1)
    } else {
      return(0)
    }
  }

  # Define the trapezoid kernel
  trapezoid_kernel =  function(xvec, c=2) {
    # x_vec is (x_1, ..., x_g)
    vals = sapply(x_vec, function(xi) max(min(c - abs(xi), 1), 0))
    return(prod(vals))
  }

  xvecl = xvec / l

  # Select the taper type
  if (type == "rectangular") {
    kappa <- rectangular_kernel(xvecl, c)
  } else if (type == "trapezoid") {
    kappa <- prod(sapply(xvecl, trapezoid_kernel, c=c))
  } else if (type == "product") {
    # Product kernel implementation can be added here
    stop("Product kernel is not implemented yet.")
  } else {
    stop("Invalid type. Choose from 'rectangular', 'trapezoid', or 'product'.")
  }
  return(kappa)
}

flat_top_taper_1d = function(x_scalar, c=1, l=1, type="rectangular") {
  # x_scalar is a single numeric
  # l is the bandwidth (scalar)
  # c is the cutoff defining "flat" portion
  # type is "rectangular" or "trapezoid"
  
  # We interpret x_scalar / l as the scaled lag
  x_scaled = x_scalar / l
  
  rectangular_kernel = function(z, c) {
    if(abs(z) <= c) 1 else 0
  }
  
  trapezoid_kernel = function(z, c) {
    #  'max(min(c - abs(z),1),0)' 
    val = c - abs(z)
    return(max(min(val,1),0))
  }
  
  if(type == "rectangular") {
    kappa = rectangular_kernel(x_scaled, c)
  } else if(type == "trapezoid") {
    kappa = trapezoid_kernel(x_scaled, c)
  } else {
    stop("type must be 'rectangular' or 'trapezoid'.")
  }
  return(kappa)
}

flat_top_taper_multi = function(x_vec, c=1, L, type="rectangular") {
  # x_vec: a numeric vector of dimension d, e.g. (h1, h2, ..., hd)
  # L:     a numeric vector of the same length as x_vec, e.g. (l1, l2, ..., ld)
  # c:     same meaning as before
  # type:  "rectangular" or "trapezoid"
  
  d = length(x_vec)
  if(length(L) != d) {
    stop("Length of L must match the length of x_vec.")
  }
  
  # We'll multiply the tapers across all dimensions
  product_kappa = 1
  for(i in seq_len(d)) {
    # For dimension i, get taper on x_vec[i] with bandwidth L[i]
    kappa_i = flat_top_taper_1d(
      x_scalar = x_vec[i],
      c        = c,
      l        = L[i],
      type     = type
    )
    product_kappa = product_kappa * kappa_i
  }
  return(product_kappa)
}



# 
# # Compute Taper vector
# compute_taper_vector = function(M_ij_matrix,
#                                 c, l, type) {
#   if (is.matrix(l)) {
#     l_vector <- as.vector(l)
#     M_ij_with_l <- cbind(M_ij_matrix, l_vector)
#     kappa_ij_vector <- sapply(1:nrow(M_ij_with_l), function(idx) {
#       flat_top_taper(M_ij_with_l[idx, 1:2], c=c, l=M_ij_with_l[idx, 3],
#                      type=type)
#     })
#   } else {
#     kappa_ij_vector <- sapply(1:nrow(M_ij_matrix), function(idx) {
#       flat_top_taper(M_ij_matrix[idx, ], c=c, l=l, type=type)
#     })
#   }
#   return(kappa_ij_vector)
# }
# 
#-------------------------------------------------------------------------------
# SPATIAL AUTOCOVARIANCE ESTIMATORS
#-------------------------------------------------------------------------------

SpatialAutoCov_v0 = function(X, hvec){
  # Estimator for the Autocovariance function (page 2 of the manuscript)
  # Input:
  #   X     matrix of spatial / spatiotemporal observations
  #   hvec  vector containing the dimensional lags
  
  h1 = hvec[1] # lag1
  h2 = hvec[2] # lag2
  
  X_mean = mean(as.vector(X))
  n1 = dim(X)[1]
  n2 = dim(X)[2]
  N = n1*n2
  
  #Finding the regions on which to shift the lags for the autocorrelation
  T11 = max(1, 1-h1)
  T12 = min(n1, n1-h1)
  T21 = max(1, 1-h2)
  T22 = min(n2, n2-h2)
  
  #Original way but we get vectors here
  # cov_lags = 1/N* sum((X[T11:T12, T21:T22] - X_mean) %*%
  #                       t(X[(T11+h1):(T12+h1), (T21+h2):(T22+h2)] - X_mean))
  
  sub_matrix1 = X[(T11+h1):(T12+h1), (T21+h2):(T22+h2)] - X_mean
  sub_matrix2 = X[T11:T12, T21:T22] - X_mean
  
  sub_vector1 = as.vector(sub_matrix1) # column wise vectorization
  sub_vector2 = as.vector(sub_matrix2) # column wise vectorization
  
  C_h_hat = 1/N* sum(t((sub_vector1) * t((sub_vector2))))
  
  return(C_h_hat)
}


SpatialAutoCov = function(X, hvec){
  # Estimator for the Autocovariance function (page 2 of the manuscript)
  # Input:
  #   X     matrix of spatial / spatiotemporal observations
  #   hvec  vector containing the dimensional lags
  
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
  
  sub_vector1 = as.vector(sub_matrix1) # column wise vectorization
  sub_vector2 = as.vector(sub_matrix2) # column wise vectorization
  
  C_h_hat = 1/N* sum(t((sub_vector1) * t((sub_vector2))))
  
  return(C_h_hat)
}


SpatialAutoCov_loop = function(X, hvec){
  # Estimator for the Autocovariance function (page 2 of the manuscript)
  # but here we do not divide by N, but rather by the number of elements 
  # in the sum.
  # NOTE: actually NOT positive semi-definite (tested on 03.09.2024)
  
  # Input:
  #   X     matrix of spatial / spatiotemporal observations
  #   hvec  vector containing the dimensional lags
  h1 = hvec[1]
  h2 = hvec[2]
  X_mean = mean(as.vector(X))
  n1 = dim(X)[1]
  n2 = dim(X)[2]
  N = n1*n2
  T11 = max(1, 1-h1)
  T12 = min(n1, n1-h1)
  T21 = max(1, 1-h2)
  T22 = min(n2, n2-h2)
  C_h_hat=0
  N_h = 0
  for(i in T21:T22){
    for (j in T11:T12){
      C_h_hat = C_h_hat + (X[i+h2, j+h1] - X_mean) * (X[i, j] - X_mean)
      N_h = N_h + 1
    }
  }
  if(N_h >0){
    nC_h_hat = 1/N_h * C_h_hat
  }else{
    nC_h_hat = 1/N * C_h_hat
  }
  return(nC_h_hat)
}
#-------------------------------------------------------------------------------
# SEPARABLE AUTOCOVARIANCE ESTIMATOR
#-------------------------------------------------------------------------------

Tapered_Sep_Autocovariance_Kron = function(X, c=1, l=1, type = "rectangular"){
  
      n1 = ncol(X)
      n2 = nrow(X)
      g = length(dim(X))
      kappas_sep = list()
      cov_sep = list()
      tapered_cov_sep  = list()
      
      # C(0,0)
      C_00 =  SpatialAutoCov(X, rep(0, g))
      
      kappa_sep_mat = matrix(NA, nrow = n2, ncol = n1)
      C_sep_mat  = matrix(NA, nrow = n2, ncol = n1)

      for(r in 1:g){
        for(i in 1:n1){
          for(j in 1:n2){
            hvec = rep(0, g)
            hvec[r] = i - j
            kappa_sep_mat[i, j] = flat_top_taper(hvec, c=c, l=l, type = type)
            C_sep_mat[i, j] = SpatialAutoCov(X, hvec)
          }
        }
        kappas_sep[[r]] = kappa_sep_mat
        cov_sep[[r]] = C_00^(-(g-1))*C_sep_mat
        tapered_cov_sep[[r]] = kappa_sep_mat * C_sep_mat
      }
      
      KronTaperCov = tapered_cov_sep[[g]]
      for(r in (g-1):1){
        KronTaperCov = kronecker(KronTaperCov, tapered_cov_sep[[r]])
      }
  return(list(KronTaperCov = KronTaperCov, kappas_sep = kappas_sep,
              cov_sep = cov_sep, tapered_cov_sep = tapered_cov_sep))
}



Tapered_Sep_Autocovariance_Kron_v2 = function(X, c=1, L=c(1,1),
                                           type = "rectangular") {
  # X is assumed to be a g-dimensional array or matrix
  #   n1 = number of columns, n2 = number of rows
  #   and g = length(dim(X)) is presumably 2
  
  if(is.vector(L) == F){
    stop("L must be a vector containing different bandwidths per spatial
         dimension.")
  }
    
    
  n1 = ncol(X)
  n2 = nrow(X)
  g  = length(dim(X))   # number of dimensions

  # placeholders
  kappas_sep       = vector("list", g)
  cov_sep          = vector("list", g)
  tapered_cov_sep  = vector("list", g)
  
  #   C_00 is the variance at lag=0,

  C_00 = SpatialAutoCov(X, rep(0, g))
  
  # We'll create matrices kappa_sep_mat, C_sep_mat to store dimension-r slices
  kappa_sep_mat = matrix(NA, nrow = n2, ncol = n1)
  C_sep_mat     = matrix(NA, nrow = n2, ncol = n1)
  
  # Loop over each dimension r in 1..g
  for(r in seq_len(g)) {
    for(i in seq_len(n1)) {
      for(j in seq_len(n2)) {
        # hvec: a g-dimensional vector of zeros
        hvec = rep(0, g)
        # only dimension r changes: hvec[r] = i-j
        hvec[r] = i - j
        
        # Apply dimension-wise product taper
        kappa_sep_mat[i, j] = flat_top_taper(xvec = hvec,
                                             c     = c,
                                             l     = L[r],    
                                             type  = type)
        
        C_sep_mat[i, j] = SpatialAutoCov(X, hvec)
      }
    }
    
    kappas_sep[[r]]       = kappa_sep_mat
    cov_sep[[r]]          = C_00^(-(g-1)) * C_sep_mat
    tapered_cov_sep[[r]]  = kappa_sep_mat * C_sep_mat
  }
  
  # Then you do your Kronecker product over all dimensions:
  KronTaperCov = tapered_cov_sep[[g]]
  for(r in seq.int(g-1, 1, by=-1)) {
    KronTaperCov = kronecker(KronTaperCov, tapered_cov_sep[[r]])
  }
  
  return(list(KronTaperCov     = KronTaperCov,
              kappas_sep      = kappas_sep,
              cov_sep         = cov_sep,
              tapered_cov_sep = tapered_cov_sep))
}



#-------------------------------------------------------------------------------
# ADAPTIVE BANDWIDTH SELECTION
#-------------------------------------------------------------------------------

# # EMPIRICAL RULE
# select_bandwidth_empirical <- function(rho_hat,
#                                        N,
#                                        C0 = 2,
#                                        K_T = max(5, sqrt(log10(T))),
#                                        c_eff = 0.5) {
#   # Bandwidth selection based on empirical rule for flat-top kernel
#   # rho_hat: estimated correlation values (matrix)
#   # N: sample (grid) size
#   # C0: constant threshold
#   # K_T: function of sample size T
#   # c_eff: effective constant for scaling
#   
#   n <- nrow(rho_hat)
#   S <- matrix(0, n, n)
#   
#   for (j in 1:n) {
#     for (k in j:n) {
#       if (j == k) {
#         # Case j = k
#         q_hat <- which.max(sapply(0:K_T, function(m) abs(rho_hat[j, j + m]) < C0 * sqrt(log10(T) / T))) - 1
#         S[j, k] <- max(ceiling(q_hat / c_eff), 1)
#       } else {
#         # Case j != k
#         q_hat_jk <- which.max(sapply(0:K_T, function(m) abs(rho_hat[j, k + m]) < C0 * sqrt(log10(T) / T))) - 1
#         q_hat_kj <- which.max(sapply(0:K_T, function(m) abs(rho_hat[k, j + m]) < C0 * sqrt(log10(T) / T))) - 1
#         q_hat <- max(q_hat_jk, q_hat_kj)
#         S[j, k] <- S[k, j] <- max(ceiling(q_hat / c_eff), 1)
#       }
#     }
#   }
#   
#   return(S)
# }

# EMPIRICAL RULE BANDWIDTH FROM POLITIS (2011)
# select_bandwidth_empirical <- function(GammaEst,
#                                        N,
#                                        C0 = 2,
#                                        K_N = max(4, sqrt(log10(N))),
#                                        c_eff = 0.5) {
#   # Bandwidth selection based on empirical rule for flat-top kernel
#   # GammaEst: estimated correlation values (matrix)
#   # N: sample size
#   # C0: constant threshold
#   # K_T: function of sample size T
#   # c_eff: effective constant for scaling
#   
#   n <- nrow(GammaEst)
#   S <- matrix(0, n, n)
#   
#   for (j in 1:n) {
#     for (k in j:n) {
#       if (j == k) {
#         # Case j = k
#         q_hat <- which.max(sapply(1:K_N, function(m) {
#           if ((j + m) > n) {
#             return(TRUE)
#           }
#           abs(GammaEst[j, j + m]) < C0 * sqrt(log10(N) / N)
#         })) - 1
#         S[j, k] <- max(ceiling(q_hat / c_eff), 1)
#       } else {
#         # Case j != k
#         q_hat_jk <- which.max(sapply(0:K_N, function(m) {
#           if ((k + m) > n) {
#             return(TRUE)
#           }
#           abs(GammaEst[j, k + m]) < C0 * sqrt(log10(N) / N)
#         })) - 1
#         
#         q_hat_kj <- which.max(sapply(0:K_N, function(m) {
#           if ((j + m) > n) {
#             return(TRUE)
#           }
#           abs(GammaEst[k, j + m]) < C0 * sqrt(log10(N) / N)
#         })) - 1
#         
#         q_hat <- max(q_hat_jk, q_hat_kj)
#         S[j, k] <- S[k, j] <- max(ceiling(q_hat / c_eff), 1)
#       }
#     }
#   }
#   return(S)
# }


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
