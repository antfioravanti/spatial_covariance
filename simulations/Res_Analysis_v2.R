#--------------------------------------------------
# 1) Load necessary packages
#--------------------------------------------------
if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
if (!require(tidyr)) install.packages("tidyr"); library(tidyr)
if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
if (!require(readr)) install.packages("readr"); library(readr)
if (!require(knitr)) install.packages("knitr"); library(knitr)
if (!require(kableExtra)) install.packages("kableExtra"); library(kableExtra)


wd = file.path(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(wd,".."))


#-------------------------------------------------------------------------------
# NORM v BANDWIDTH - with two levels: beta
#-------------------------------------------------------------------------------

# Read data
resdir <- file.path(wd, "results")
df <- read.csv(file.path(resdir, "res_grid_search_all_01-04-2025.csv"), sep = ";")

# Set filtering parameters:
alpha_filter  <- 1
type_filter   <- "rectangular"
lambda_sep = 6
lambda_nonsep = 10

# We want to include both beta = 0 with lambda = 4 and beta = 1 with lambda = 8
df_filtered <- df %>%
  filter(
    alpha    == alpha_filter,
    type     == type_filter,
    GridSize %in% c(20, 30, 50),
    ((beta == 0) & (lambda == lambda_sep)) | ((beta == 1) & (lambda == lambda_nonsep))
  )

# Group by GridSize, beta, and l to compute mean performance for each norm
df_summarized <- df_filtered %>%
  group_by(GridSize, beta, l) %>%
  summarise(
    meanSNorm             = mean(SNorm, na.rm = TRUE),
    meanSNormTapered      = mean(SNormTapered, na.rm = TRUE),
    meanSNormSeparTapered = mean(SNormSeparTapered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(
    Naive = meanSNorm,
    Tapered = meanSNormTapered,
    SeparableTapered = meanSNormSeparTapered
  )

# Reshape the data to long format
df_long <- df_summarized %>%
  pivot_longer(
    cols = c("Naive", "Tapered", "SeparableTapered"),
    names_to = "Estimator",
    values_to = "Value"
  )

# Identify the minimum Value for each combination of GridSize, beta, and Estimator
df_min <- df_long %>%
  group_by(GridSize, beta, Estimator) %>%
  slice_min(Value, n = 1, with_ties = FALSE) %>%
  ungroup()

# Define colors for the estimators
estimator_colors <- c("Naive" = "black", "Tapered" = "blue", "SeparableTapered" = "red")

# Set up a 2-row by 3-column plotting layout and reserve outer bottom margin space for the legend
par(mfrow = c(2, 3), oma = c(8, 0, 0, 0), mar = c(4, 4, 3, 1))

# Get unique beta values and grid sizes (ordered)
beta_values <- sort(unique(df_long$beta))   # e.g., 0 and 1
grid_sizes  <- sort(unique(df_long$GridSize))  # e.g., 20, 30, 50

# Loop over beta (rows) and GridSize (columns) to create the plots
for(b in beta_values) {
  for(gs in grid_sizes) {
    # Subset data for this beta and grid size
    df_sub <- subset(df_long, beta == b & GridSize == gs)
    
    # Set x and y ranges
    x_range <- range(df_sub$l, na.rm = TRUE)
    y_range <- range(df_sub$Value, na.rm = TRUE)
    
    # Create an empty plot
    plot(x_range, y_range, type = "n",
         main = paste("β =", b, "\nGridSize =", gs, "x", gs),
         xlab = "l", ylab = "Mean Norm Diff")
    
    # For each estimator, plot its line and points
    for(est in names(estimator_colors)) {
      df_est <- subset(df_sub, Estimator == est)
      if(nrow(df_est) > 0) {
        ord <- order(df_est$l)
        lines(df_est$l[ord], df_est$Value[ord], col = estimator_colors[est], lwd = 2)
        points(df_est$l, df_est$Value, col = estimator_colors[est], pch = 16)
        
        # Add label for the minimum point
        df_est_min <- subset(df_min, beta == b & GridSize == gs & Estimator == est)
        
        if(nrow(df_est_min) > 0) {
          pos_val <- ifelse(est == "Naive", 1, 3)  # Below for Naive, above for others
          x_coord <- if(est == "Naive") df_est_min$l + 5 else df_est_min$l
          text(x_coord, df_est_min$Value,
               labels = round(df_est_min$Value, 2),
               pos = pos_val, cex = 0.9, col = estimator_colors[est])
        }
      }
    }
  }
}

# Create a new plotting region at the bottom for the legend
par(xpd = NA)  # Allow drawing outside the plot region
par(fig = c(0, 1, 0, 0.1), new = TRUE, mar = c(0,0,0,0))
plot.new()
legend(x = 0.5, y = -0.5, legend = names(estimator_colors),
       col = estimator_colors, lwd = 2, pch = 16, horiz = TRUE,
       bty = "n", cex = 1.2, xjust = 0.5, yjust = 0.5)
par(xpd = FALSE)

#-------------------------------------------------------------------------------
# GRID SIZE V NORM
#-------------------------------------------------------------------------------

# Read data
resdir <- file.path(wd, "results")
df <- read.csv(file.path(resdir, "res_grid_search_all_01-04-2025.csv"), sep = ";")

# Set filtering parameters:
alpha_filter   <- 1
type_filter    <- "trapezoid"
lambda_sep     <- 2    # for beta = 0 (separable)
lambda_nonsep  <- 6   # for beta = 1 (nonseparable)

# Choose one bandwidth value for each case:
l_filter_beta0 <- 4   # chosen bandwidth for separable case (β = 0)
l_filter_beta1 <- 4   # chosen bandwidth for nonseparable case (β = 1)

# Filter data: we select rows corresponding to the chosen l values
df_filtered <- df %>%
  filter(
    alpha    == alpha_filter,
    type     == type_filter,
    GridSize != 10,
    (
      (beta == 0 & lambda == lambda_sep & l == l_filter_beta0) |
        (beta == 1 & lambda == lambda_nonsep & l == l_filter_beta1)
    )
  )

# Summarize: group by GridSize and beta, then compute means and SDs
df_summarized <- df_filtered %>%
  group_by(GridSize, beta) %>%
  summarise(
    Naive_mean = mean(SNorm, na.rm = TRUE),
    Naive_sd   = sd(SNorm, na.rm = TRUE),
    Tapered_mean = mean(SNormTapered, na.rm = TRUE),
    Tapered_sd   = sd(SNormTapered, na.rm = TRUE),
    SeparableTapered_mean = mean(SNormSeparTapered, na.rm = TRUE),
    SeparableTapered_sd   = sd(SNormSeparTapered, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  # Compute standard errors and 95% confidence intervals using t quantiles:
  mutate(
    Naive_se = Naive_sd / sqrt(n),
    Tapered_se = Tapered_sd / sqrt(n),
    SeparableTapered_se = SeparableTapered_sd / sqrt(n),
    Naive_lower = Naive_mean - qt(0.975, df = n - 1) * Naive_se,
    Naive_upper = Naive_mean + qt(0.975, df = n - 1) * Naive_se,
    Tapered_lower = Tapered_mean - qt(0.975, df = n - 1) * Tapered_se,
    Tapered_upper = Tapered_mean + qt(0.975, df = n - 1) * Tapered_se,
    SeparableTapered_lower = SeparableTapered_mean - qt(0.975, df = n - 1) * SeparableTapered_se,
    SeparableTapered_upper = SeparableTapered_mean + qt(0.975, df = n - 1) * SeparableTapered_se
  )

# Reshape the data: pivot_longer so that we have one row per GridSize, beta, and Estimator.
# We use a names pattern to split the column names into Estimator and Statistic.
df_long <- df_summarized %>%
  pivot_longer(
    cols = c(Naive_mean, Naive_lower, Naive_upper,
             Tapered_mean, Tapered_lower, Tapered_upper,
             SeparableTapered_mean, SeparableTapered_lower, SeparableTapered_upper),
    names_to = c("Estimator", ".value"),
    names_pattern = "(.*)_(.*)"
  )
# Now df_long has columns: GridSize, beta, Estimator, mean, lower, upper.

# Define colors for the estimators
estimator_colors <- c("Naive" = "black", "Tapered" = "blue", "SeparableTapered" = "red")

# Set up a 1-row by 2-column layout for the two β values.
par(mfrow = c(1, 2), oma = c(4, 0, 2, 0), mar = c(5, 5, 4, 2))

# Define translucent fill colors for CIs
ci_fill <- c("Naive" = adjustcolor("black", alpha.f = 0.15),
             "Tapered" = adjustcolor("blue", alpha.f = 0.15),
             "SeparableTapered" = adjustcolor("red", alpha.f = 0.15))

# Get unique beta values (expected: 0 and 1)
beta_values <- sort(unique(df_long$beta))

# Loop over each beta value to create a panel
for (b in beta_values) {
  df_sub <- subset(df_long, beta == b)
  
  # Set x and y axis limits
  x_range <- range(df_sub$GridSize, na.rm = TRUE)
  y_range <- range(df_sub$lower, df_sub$upper, na.rm = TRUE)
  
  # Empty plot
  plot(x_range, y_range, type = "n",
       main = paste("β =", b),
       xlab = "GridSize", ylab = "Mean Spectral Norm Diff",
       xaxt = "n")
  axis(1, at = unique(df_sub$GridSize))
  
  for (est in names(estimator_colors)) {
    df_est <- subset(df_sub, Estimator == est)
    if (nrow(df_est) > 0) {
      ord <- order(df_est$GridSize)
      x_vals <- df_est$GridSize[ord]
      y_mean <- df_est$mean[ord]
      y_lower <- df_est$lower[ord]
      y_upper <- df_est$upper[ord]
      
      # Draw shaded CI ribbon using polygon
      polygon(c(x_vals, rev(x_vals)),
              c(y_lower, rev(y_upper)),
              col = ci_fill[est], border = NA)
      
      # Draw mean line
      lines(x_vals, y_mean, col = estimator_colors[est], lwd = 2)
      points(x_vals, y_mean, col = estimator_colors[est], pch = 16)
      
      # (Optional annotation removed)
    }
  }
}

# Add legend manually at specific coordinates
par(xpd = NA)
# Use a smaller top boundary for the legend region (e.g., 0.1 becomes 0.08)
par(fig = c(0, 1, 0, 0.08), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()
usr <- par("usr")
# Lower the legend by subtracting a small offset (e.g., 0.05 times the vertical extent)
legend(x = (usr[1] + usr[2]) / 2, 
       y = usr[3] + 0.5 * (usr[4] - usr[3]) - 0.05 * (usr[4] - usr[3]),
       legend = names(estimator_colors),
       col = estimator_colors, lwd = 2, pch = 16, horiz = TRUE,
       bty = "n", cex = 1.2, xjust = 0.5)
par(xpd = FALSE)
