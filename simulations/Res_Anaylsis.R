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
#------------------------------------------------------------------------------
# GridSize, show how the mean norms vary with l
#------------------------------------------------------------------------------
resdir = file.path(wd, "results")
df = read.csv(file.path(resdir, "res_grid_search_all_01-04-2025.csv"),
              sep = ";")

# 3) Set the parameters you want to filter
#    Adjust these values to match your desired subset
alpha_filter  = 1
beta_filter   = 1
lambda_filter = 8
type_filter   = "rectangular"

# 4) Filter data for chosen alpha, beta, lambda, type
#    Exclude GridSize = 10
df_filtered <- df %>%
  filter(
    alpha    == alpha_filter,
    beta     == beta_filter,
    lambda   == lambda_filter,
    type     == type_filter,
    GridSize == c(20,30,50)
  )

# 5) Group by GridSize and l, then compute mean of each norm
df_summarized <- df_filtered %>%
  group_by(GridSize, l) %>%
  summarise(
    meanSNorm             = mean(SNorm, na.rm = TRUE),
    meanSNormTapered      = mean(SNormTapered, na.rm = TRUE),
    meanSNormSeparTapered = mean(SNormSeparTapered, na.rm = TRUE),
    .groups = "drop"
  )


# 6) Reshape (pivot) the data to long format for ggplot
df_summarized = df_summarized %>% 
            rename(
              Naive = meanSNorm,
              Tapered = meanSNormTapered,
              SeparableTapered = meanSNormSeparTapered
  )

df_long <- df_summarized %>%
  pivot_longer(
    cols      = c("Naive", "Tapered", "SeparableTapered"),
    names_to  = "Estimator",
    values_to = "Value"
  )

# 7) Identify the unique lowest Value for each Estimator within each GridSize
df_min <- df_long %>%
  group_by(GridSize, Estimator) %>%
  slice_min(Value, n = 1, with_ties = FALSE) %>%
  ungroup()

par(mfrow=c(3,3))
# 8) Plot: For each GridSize, show how the mean norms vary with l
#    and label the lowest point for each estimator
ggplot(df_long, aes(x = l, y = Value, color = Estimator)) +
  geom_line(size = 1.2) +   # Thicker lines
  geom_point(size = 2) +    # Points for clarity
  # Label the lowest points for non-Naive (placed above)
  geom_text(
    data = df_min %>% filter(Estimator != "Naive"),
    aes(label = round(Value, 2)),
    vjust = -0.5,
    size = 4,
    show.legend = FALSE
  ) +
  # Label the lowest point for Naive (placed below)
  geom_text(
    data = df_min %>% filter(Estimator == "Naive"),
    aes(label = round(Value, 2)),
    vjust = 1.5,
    size = 4,
    show.legend = FALSE
  ) +
  facet_wrap(~ GridSize, scales = "free_y",
             labeller = labeller(GridSize = function(x) paste("N =", x, "x", x))) +
  scale_x_continuous(
    breaks = unique(df_long$l),  # Use all unique l values as breaks
    labels = function(x) format(x, nsmall = 0),  # Ensure integer labels
    minor_breaks = NULL  # Remove minor breaks for clarity
  ) +
  labs(
    title = "",
    x     = "l",
    y     = "Mean Spectral Norm Difference"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),  # Larger x-axis text
    axis.text.y = element_text(size = 12),  # Larger y-axis text
    axis.title.x = element_text(size = 14),  # Larger x-axis title
    axis.title.y = element_text(size = 14),  # Larger y-axis title
    legend.position = "bottom",  # Move legend below the plot
    legend.text = element_text(size = 12),  # Larger legend text
    legend.title = element_text(size = 14),  # Larger legend title (if used)
    strip.text = element_text(size = 14, face = "bold")  # Increase facet titles size
  )


#-------------------------------------------------------------------------------
# PLOT OF GridSize vs. Spectral Norms with Confidence Intervals
#-------------------------------------------------------------------------------

# Filter the data for the selected parameters

l_filter      <- 6  # Example value, adjust as needed
alpha_filter  <- 1
beta_filter   <- 1
lambda_filter <- 10
type_filter   <- "rectangular"

df_filtered <- df %>%
  filter(
    l       == l_filter,
    alpha   == alpha_filter,
    beta    == beta_filter,
    lambda  == lambda_filter,
    type    == type_filter
  )
# Compute mean, standard deviation, and confidence intervals

df_summarized <- df_filtered %>%
  group_by(GridSize) %>%
  summarise(
    meanSNorm             = mean(SNorm, na.rm = TRUE),
    sdSNorm               = sd(SNorm, na.rm = TRUE),
    meanSNormTapered      = mean(SNormTapered, na.rm = TRUE),
    sdSNormTapered        = sd(SNormTapered, na.rm = TRUE),
    meanSNormSeparTapered = mean(SNormSeparTapered, na.rm = TRUE),
    sdSNormSeparTapered   = sd(SNormSeparTapered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Compute confidence intervals (assuming 100 samples per configuration)
  mutate(
    ciSNorm_lower         = meanSNorm - 1.96 * (sdSNorm / sqrt(100)),
    ciSNorm_upper         = meanSNorm + 1.96 * (sdSNorm / sqrt(100)),
    ciTapered_lower       = meanSNormTapered - 1.96 * (sdSNormTapered / sqrt(100)),
    ciTapered_upper       = meanSNormTapered + 1.96 * (sdSNormTapered / sqrt(100)),
    ciSeparTapered_lower  = meanSNormSeparTapered - 1.96 * (sdSNormSeparTapered / sqrt(100)),
    ciSeparTapered_upper  = meanSNormSeparTapered + 1.96 * (sdSNormSeparTapered / sqrt(100))
  )


# Reshape (pivot) the data to long format for ggplot
df_long <- df_summarized %>%
  pivot_longer(
    cols = c(meanSNorm, meanSNormTapered, meanSNormSeparTapered),
    names_to = "Estimator",
    values_to = "Value"
  ) %>%
  mutate(
    Estimator = recode(Estimator,
                       "meanSNorm" = "Classical",
                       "meanSNormTapered" = "Tapered",
                       "meanSNormSeparTapered" = "SeparableTapered"),
    CI_lower = case_when(
      Estimator == "Classical" ~ ciSNorm_lower,
      Estimator == "Tapered" ~ ciTapered_lower,
      Estimator == "SeparableTapered" ~ ciSeparTapered_lower
    ),
    CI_upper = case_when(
      Estimator == "Classical" ~ ciSNorm_upper,
      Estimator == "Tapered" ~ ciTapered_upper,
      Estimator == "SeparableTapered" ~ ciSeparTapered_upper
    )
  )


# 6) Plot: GridSize vs. Spectral Norms with Confidence Intervals
ggplot(df_long, aes(x = GridSize, y = Value, color = Estimator)) +
  # Add confidence interval as a shaded region
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Estimator), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +   # Thicker lines for mean curves
  geom_point(size = 2) +    # Points for clarity
  geom_linerange(aes(ymin = CI_lower, ymax = CI_upper), linetype = "dotted", size = 0.6) + # Dotted CIs
  labs(
    # title = paste0("Spectral Norms vs. GridSize (l = ", l_filter, 
    #                ", α = ", alpha_filter, ", β = ", beta_filter, 
    #                ", λ = ", lambda_filter, ")"),
    x = "Grid Size",
    y = "Mean Spectral Norm Difference"
  ) +
  scale_x_continuous(
    breaks = unique(df_long$GridSize)  # Ensure only actual GridSize values are used as breaks
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Larger x-axis text
    axis.text.y = element_text(size = 12),  # Larger y-axis text
    axis.title.x = element_text(size = 14),  # Larger x-axis title
    axis.title.y = element_text(size = 14),  # Larger y-axis title
    legend.position = c(0.5, -0.2),  # Move legend below the plot
    legend.text = element_text(size = 12),  # Larger legend text
    legend.title = element_text(size = 14),  # Larger legend title (if used)
    strip.text = element_text(size = 14, face = "bold")  # Increase facet titles size
  )
#-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
# SUMMARY TABLES
#-------------------------------------------------------------------------------

alpha_filter  <- 1
beta_filter   <- 0
lambda_filter <- 4
type_filter   <- "rectangular"

# 4) Filter data for chosen alpha, beta, lambda, type; exclude GridSize = 10
df_filtered <- df %>%
  filter(
    alpha    == alpha_filter,
    beta     == beta_filter,
    lambda   == lambda_filter,
    type     == type_filter
  )

# 5) Group by GridSize and l, then compute mean and standard deviation for each norm
df_summarized <- df_filtered %>%
  group_by(GridSize, l) %>%
  summarise(
    meanSNorm             = mean(SNorm, na.rm = TRUE),
    sdSNorm               = sd(SNorm, na.rm = TRUE),
    meanSNormTapered      = mean(SNormTapered, na.rm = TRUE),
    sdSNormTapered        = sd(SNormTapered, na.rm = TRUE),
    meanSNormSeparTapered = mean(SNormSeparTapered, na.rm = TRUE),
    sdSNormSeparTapered   = sd(SNormSeparTapered, na.rm = TRUE),
    .groups = "drop"
  )

# 6) Reshape the table into a wide format

# For the naïve (Classical) estimator, average over l for each GridSize
df_naive <- df_summarized %>%
  group_by(GridSize) %>%
  summarise(
    Classical_mean = mean(meanSNorm, na.rm = TRUE),
    Classical_sd   = mean(sdSNorm, na.rm = TRUE)
  ) %>%
  ungroup()

# For the Tapered estimator, select the row with the minimum meanSNormTapered per GridSize
df_tapered <- df_summarized %>%
  group_by(GridSize) %>%
  slice_min(meanSNormTapered, with_ties = FALSE) %>%
  ungroup() %>%
  select(GridSize, optimal_l_taper = l, Tapered_mean = meanSNormTapered, Tapered_sd = sdSNormTapered)

# For the Separable & Tapered estimator, select the row with the minimum meanSNormSeparTapered per GridSize
df_sep <- df_summarized %>%
  group_by(GridSize) %>%
  slice_min(meanSNormSeparTapered, with_ties = FALSE) %>%
  ungroup() %>%
  select(GridSize, optimal_l_sep = l, SeparableTapered_mean = meanSNormSeparTapered, SeparableTapered_sd = sdSNormSeparTapered)

# Join the three summaries by GridSize to get one wide table
df_table <- df_naive %>%
  left_join(df_tapered, by = "GridSize") %>%
  left_join(df_sep, by = "GridSize")

# Optionally, view the reshaped table
print(df_table)

latex_table <- kable(
  df_table,
  format = "latex",
  booktabs = TRUE,
  escape= FALSE,
  col.names = c("Grid Size", "Mean", "SD",
                "$l_{\\text{tap}}$", "Mean", "SD",
                "$l_{\\text{sep}}$", "Mean", "SD"),
  caption = "Simulation Results for Different Estimators and GridSizes"
) %>%
  add_header_above(c(" " = 1,
                     "Naïve" = 2,
                     "Tapered " = 3,
                     "Separable \\& Tapered" = 3)) %>%
  kable_styling(latex_options =  c("hold_position"), full_width = FALSE)

cat(latex_table)
