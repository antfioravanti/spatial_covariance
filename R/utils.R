#-------------------------------------------------------------------------------
### UTILS FUNCTIONS ###
#-------------------------------------------------------------------------------
# Function to check that a matrix is positive semi definite
is_positive_semi_definite = function(matrix) {
  eigenvalues <- eigen(matrix)$values
  return(all(eigenvalues > 0))
}
# 
# # Function to check if covariance matrix is spatially separable
# check_separability_svd <- function(cov_matrix, tolerance = 1e-6) {
#   # Perform singular value decomposition (SVD)
#   svd_result <- svd(cov_matrix)
#   
#   # Singular values from the decomposition
#   singular_values <- svd_result$d
#   
#   # Calculate the ratio of the largest singular value to the sum of all singular values
#   largest_singular_value_ratio <- singular_values[1] / sum(singular_values)
#   
#   # Check if the second largest singular value is below a tolerance
#   is_separable <- singular_values[2] < tolerance
#   
#   # Output result
#   if (is_separable) {
#     message("The covariance matrix is approximately separable.")
#   } else {
#     message("The covariance matrix is not separable.")
#   }
#   
#   # Return details for further inspection
#   return(list(
#     singular_values = singular_values,
#     largest_singular_value_ratio = largest_singular_value_ratio,
#     is_separable = is_separable
#   ))
# }

my_cpp_source_funs = function(file) {
  tmp = new.env(parent=parent.frame())
  sourceCpp(file, env = tmp)
  funs <- names(tmp)[unlist(eapply(tmp, is.function))]
  for(x in names(tmp)) {
    assign(x, tmp[[x]], envir = parent.frame())
  }
  return(list(functions=funs))
}