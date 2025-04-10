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
if (!require(wesanderson)) install.packages("wesanderson"); library(wesanderson)
if (!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)
if (!require(pals)) install.packages("pals"); library(pals)
if (!require(microbenchmark)) install.packages("microbenchmark"); library(microbenchmark)
if (!require(wesanderson)) install.packages("wesanderson"); library(wesanderson)
if (!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)
if (!require(pals)) install.packages("pals"); library(pals)
#if (!require(here)) install.packages("here"); library(here)
# Set the working directory to where the current file is saved
wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(wd,".."))
source("R/utils.R") # General functions
source("R/covariance_funs.R") # Analytical covariances
source("R/estimators.R") # Estimators
source("R/plotting.R")  # Plotting functions
#-------------------------------------------------------------------------------

wes_pal = wes_palette("Zissou1", n=21, type = c( "continuous"))
cols = brewer.pal(3, "BuGn")
pal = colorRampPalette(cols)
coolwarm_pal = coolwarm(10)
#-------------------------------------------------------------------------------
load(file.path(wd, "satellite-data.RData"))

# Check the names of the list elements
names(data)

# Create a layout: 2 rows x 5 columns, where the 10th panel is for the legend
layout_matrix <- matrix(c(1,2,3,4,9, 5,6,7,8,9), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# Plot the 8 NDVI matrices in panels 1-8
for(i in 1:8) {
  date_name <- names(data)[i]
  ndvi_matrix <- data[[date_name]]
  
  plot_matrix(ndvi_matrix,
              main = date_name, 
              col = rev(wes_pal),
              labels = TRUE)
}

# Now add a common legend in panel 9.
# Determine the range of NDVI values across your matrices
ndvi_range <- range(unlist(data[1:8]), na.rm = TRUE)

# Adjust margins if needed
par(mar = c(5, 4, 2, 2))
image.plot(legend.only = TRUE, 
           zlim = ndvi_range, 
           col = rev(wes_pal), 
           legend.width = 1, 
           legend.mar = 4)
#-------------------------------------------------------------------------------
dates_to_plot <- c("2013-08-09", "2019-07-09")

# Set up a layout: 1 row with 3 columns.
layout_matrix <- matrix(c(1, 2, 3), nrow = 1, byrow = TRUE)
layout(layout_matrix, widths = c(2, 2, 0.5))

# Plot the NDVI matrices for each date in panels 1 and 2.
for(date_name in dates_to_plot) {
  par(mar = c(4, 4, 2, 0.2))
  ndvi_matrix <- data[[date_name]]
  plot_matrix(ndvi_matrix,
              main = date_name, 
              col = rev(wes_pal),
              labels = TRUE)
}

# Calculate the overall NDVI range from the selected matrices.
ndvi_range <- range(unlist(data[dates_to_plot]), na.rm = TRUE)

# Switch to the legend panel (panel 3).
# Increase the right margin to push the legend further to the left.
par(mar = c(4, 0, 2, 5))
plot.new()  # Start a fresh plot in the legend panel

library(fields)
image.plot(legend.only = TRUE, 
           zlim = ndvi_range, 
           col = rev(wes_pal), 
           horizontal = FALSE,   # vertical legend
           legend.width = 10,     # internal legend width
           legend.mar = 5,        # increased internal margin for the legend
           axis.args = list(cex.axis = 1))



#-------------------------------------------------------------------------------

# Define the functions for each taper

# Rectangular taper: 1 for |x| <= 1, 0 otherwise.
rect_taper <- function(x) {
  ifelse(abs(x) <= 1, 1, 0)
}

# Trapezoidal taper: 1 for |x| <= c, 0 for |x| > c.
# With c_val = 2, this becomes 1 for |x| <= 2 and 0 for |x| > 2.
trap_taper <- function(x, c_val) {
  ifelse(abs(x) <= 1, 1,
         ifelse((abs(x)> 1 & (abs(x)<= c_val)), c_val - abs(x), 0)
  )
}

# Differentiable taper as defined:
# lambda(s) = 1 if |s| <= c
#           = exp( -b * exp(-b/ (|s|-c)^2) / (|s|-1)^2 ) if c < |s| < 1
#           = 0 if |s| >= 1
diff_taper <- function(x, b, c_val) {
  ifelse(abs(x) <= c_val, 1,
         ifelse(abs(x) >= 1, 0,
                exp(-b * exp(-b / (abs(x) - c_val)^2) / (abs(x) - 1)^2)
         ))
}

# Set parameter values:
b <- 2
c_rect <- 1      # for rectangular taper (though not needed in function)
c_trap <- 2      # for trapezoidal taper, here it becomes a "wider" rectangle
c_diff <- 0.05   # for differentiable taper

# Define the evaluation grid (common domain for all plots)
x <- seq(-2.5, 2.5, length.out = 1000)

# Compute the taper values:
y_rect   <- rect_taper(x)
y_trap   <- trap_taper(x, c_trap)
y_diff   <- diff_taper(x, b, c_diff)

# Set up a 1x3 plotting area:
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

# Plot Rectangular Taper
plot(x, y_rect, type = "l", col = "blue", lwd = 2, ylim = c(0, 1.2),
     xlab = "x", ylab = "Taper Value", main = "Rectangular Taper")

# Plot Trapezoidal Taper
plot(x, y_trap, type = "l", col = "red", lwd = 2, ylim = c(0, 1.2),
     xlab = "x", ylab = "Taper Value", main = "Trapezoidal Taper")

# Plot Differentiable Taper
plot(x, y_diff, type = "l", col = "darkgreen", lwd = 2, ylim = c(0, 1.2),
     xlab = "x", ylab = "Taper Value", main = "Differentiable Taper")

# Reset plotting layout to default
par(mfrow = c(1, 1))
