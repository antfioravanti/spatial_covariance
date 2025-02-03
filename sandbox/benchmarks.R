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




# Benchmark both versions
benchmark_results <- microbenchmark(
  R_version = compute_m_values(5, 3, nvec, d = 1, flipsign = TRUE, flipposition = FALSE),
  Cpp_version = compute_m_values_cpp(5, 3, nvec, d = 1, flipsign = TRUE, flipposition = FALSE),
  times = 10000
)

print(benchmark_results)

# Test input
i <- 5
j <- 3
nvec <- c(20,20)
d <- 1
flipsign <- TRUE
flipposition <- FALSE

# Run the R function
result_r <- compute_m_values(i, j, nvec, d, flipsign, flipposition)

# Run the C++ function
result_cpp <- compute_m_values_cpp(i, j, nvec, d, flipsign, flipposition)


i <- 5
j <- 3

d <- 1
flipsign <- TRUE
flipposition <- FALSE


k = -1
s = 2
d = 1
m_fun(k = k, s=s, d=1, nvec=nvec)
m_fun_cpp(k = k, s=s, d=1, nvec=nvec)



prod_prev_n = prod(nvec[1:(s-1)])
m = ((ceiling(k / (d * prod_prev_n)) - 1) %% 
       (d * prod_prev_n * nvec[s])) + 1


