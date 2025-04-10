if (!require(Rcpp)) install.packages("Rcpp"); library(Rcpp)
if (!require(wesanderson)) install.packages("wesanderson"); library(wesanderson)
if (!require(RColorBrewer)) install.packages("RColorBrewer"); library(RColorBrewer)

sourceCpp("src/estimators.cpp")

#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

plot_matrix = function(X, nvec,
                       main = "Matrix", 
                       labels = FALSE,
                       show_cell_num = FALSE,
                       show_index = FALSE, 
                       show_M_num = FALSE, 
                       show_values = FALSE,
                       col_palette = NA,
                       show_legend = FALSE){
  
  if((length(col_palette) == 1 && is.na(col_palette))){
    # Plot Matrix by keeping fixed the cell positions
    image(t(X)[, nrow(X):1], xaxt = "n", yaxt = "n", main = main)
    
  }else{
    # Plot Matrix by keeping fixed the cell positions
    # with custom color palette
    image(t(X)[, nrow(X):1], col = col_palette,
          xaxt = "n", yaxt = "n", main = main)
  }

  if(labels){
    # Add x-axis labels
    axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
    
    # Add y-axis labels
    axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
  }
  
  if (show_cell_num) {
    # Add cell numbers inside each cell (row-wise numbering)
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        cell_num = (i - 1) * ncol(X) + j  
        text(
          x = (j - 1) / (ncol(X) - 1), 
          y = 1 - (i - 1) / (nrow(X) - 1), 
          labels = cell_num, 
          col = "black", 
          cex = 0.8)
      }
    }
  }
  
  if (show_index){
    # Add index labels (i,j)
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        text(x = x_coord, y = y_coord, labels = paste0("(", i, ",", j, ")"), 
             col = "black", cex = 0.6)
      }
    }
  }
  
  if (show_M_num){
    if(is.null(nvec)){
      stop("Provide an nvec")
    }
    
    N = prod(nvec)
    # compute_M_matrix_cpp and compute_m_values_cpp are assumed to be available functions
    M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        # Compute the text placement
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        m_values = compute_m_values_cpp(i, j, nvec = nvec)
        M_values = paste(m_values, collapse = ", ")
        text(x = x_coord, y = y_coord, labels = M_values, col = "black",
             cex = 0.6)
      }
    }
  }
  
  if(show_values){
    # Print matrix values (up to 4 decimal places) in each cell
    for (i in seq_len(nrow(X))) {
      for (j in seq_len(ncol(X))) {
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        val_str = sprintf("%.4f", X[i, j])
        text(
          x = x_coord,
          y = y_coord,
          labels = val_str,
          col = "black",
          cex = 0.7)
      }
    }
  }
  
  if(show_legend){
    # Add a color legend (color bar) using the fields package.
    # The legend will reflect the range of values in X.
    if (!requireNamespace("fields", quietly = TRUE)) {
      stop("Please install the 'fields' package to display a legend (install.packages('fields')).")
    }
    # Adjust the legend margins and width as needed.
    fields::image.plot(legend.only = TRUE, zlim = range(X, na.rm = TRUE),
                       col = col_palette, legend.width = 1, legend.mar = 4.1)
  }
}



plot_matrix_grayscale = function(X){
  nvec = dim(X)
  return(filled.contour(x=1:nvec[1], y=1:nvec[2],
                        X, color.palette=gray.colors))
}
plot_3D_matrix = function(X){
  nvec = dim(X)
  return(persp(x=(1:nvec[1]), y=(1:nvec[2]), X,
               theta=45, phi=35, r=5, expand=0.6, axes=T,
               ticktype="detailed", xlab="t1 ", ylab="t2", zlab="X_t1,t2"))
}
