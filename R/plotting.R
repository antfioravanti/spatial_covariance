if (!require(Rcpp)) install.packages("Rcpp"); library(Rcpp)
sourceCpp("src/estimators.cpp")

#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

plot_matrix = function(X, nvec,
                       main = "Matrix", 
                       labels = F,
                       show_cell_num = F,
                       show_index = F, 
                       show_M_num = F, 
                       show_values = F){
  # Plot Matrix by keeping fixed the cell positions
  
  image(t(X)[, nrow(X):1], xaxt = "n", yaxt = "n", main = main)
  
  if(labels == T){
    # Add x-axis labels
    axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
    
    # Add y-axis labelscell_num
    axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
  }
  
  if (show_cell_num) {
    # Add numbers inside each cell
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        cell_num = (i - 1) * ncol(X) + j  # Row-wise numbering
        text(
          x = (j - 1) / (ncol(X) - 1), 
          y = 1 - (i - 1) / (nrow(X) - 1), 
          labels = cell_num, 
          col = "black", 
          cex = 0.8)  }}
  }
  
  if (show_index){
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        text(x = x_coord, y = y_coord, labels = paste0("(", i, ",", j, ")"), 
             col = "black", cex = 0.6)  }}
  }
  
  if (show_M_num){
    if(is.null(nvec)){
      stop("Provide an nvec")
    }
    
    N = prod(nvec)
    M_ij_matrix = compute_M_matrix_cpp(N, nvec = nvec)
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        # Compute the text placement
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        # Compute row-wise numbering
        m_values = compute_m_values_cpp(i, j, nvec = nvec)
        M_values = paste(m_values, collapse = ", ")
        text(x = x_coord, y = y_coord, labels = M_values, col = "black",
             cex = 0.6)   }}
  }
  
  
  if(show_values){
    # Print matrix values (up to 4 decimal places) in each cell
    for (i in seq_len(nrow(X))) {
      for (j in seq_len(ncol(X))) {
        x_coord = (j - 1) / (ncol(X) - 1)
        y_coord = 1 - (i - 1) / (nrow(X) - 1)
        val_str = sprintf("%.4f", X[i, j])  # up to 4 decimal places
        text(
          x = x_coord,
          y = y_coord,
          labels = val_str,
          col = "black",
          cex = 0.7)  }}
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
