
#-------------------------------------------------------------------------------
# PLOTTING
#-------------------------------------------------------------------------------

plot_matrix = function(X, main = "Matrix", labels = F){
  # Plot Matrix by keeping fixed the cell positions
  
  image(t(X)[, nrow(X):1], xaxt = "n", yaxt = "n", main = main)
  
  if(labels == T){
  # Add x-axis labels
  axis(1, at = seq(0, 1, length.out = ncol(X)), labels = 1:ncol(X))
  
  # Add y-axis labels
  axis(2, at = seq(0, 1, length.out = nrow(X)), labels = nrow(X):1)
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
